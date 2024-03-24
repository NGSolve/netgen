#include "debugging.hpp"

#include <libusockets.h>

#include <core/python_ngcore.hpp>

#include "pybind11/gil.h"

#ifdef NETGEN_DEBUGGING_GUI
#include "websockets/App.h"
#endif  // NETGEN_DEBUGGING_GUI

namespace netgen {
unique_ptr<Mesh> GetOpenElements(const Mesh& m, int dom) {
  static Timer t("GetOpenElements");
  RegionTimer rt(t);
  auto mesh = make_unique<Mesh>();
  *mesh = m;

  Array<bool, PointIndex> interesting_points(mesh->GetNP());
  interesting_points = false;

  mesh->FindOpenElements(dom);
  NgArray<Element2d> openelements;
  openelements = mesh->OpenElements();

  for (auto& el : openelements)
    for (auto i : el.PNums()) interesting_points[i] = true;

  for (auto& el : mesh->VolumeElements()) {
    int num_interesting_points = 0;

    for (auto pi : el.PNums())
      if (interesting_points[pi]) num_interesting_points++;

    if (num_interesting_points == 0) el.Delete();
    el.SetIndex(num_interesting_points);
  }

  mesh->SetMaterial(1, "1_point");
  mesh->SetMaterial(2, "2_points");
  mesh->SetMaterial(3, "3_points");
  mesh->SetMaterial(4, "4_points");
  mesh->Compress();

  mesh->ClearSurfaceElements();

  for (auto& el : openelements) mesh->AddSurfaceElement(el);

  return mesh;
}

unique_ptr<Mesh> FilterMesh(const Mesh& m, FlatArray<PointIndex> points,
                            FlatArray<SurfaceElementIndex> sels,
                            FlatArray<ElementIndex> els) {
  static Timer t("GetOpenElements");
  RegionTimer rt(t);
  auto mesh_ptr = make_unique<Mesh>();
  auto& mesh = *mesh_ptr;
  mesh = m;

  Array<bool, PointIndex> keep_point(mesh.GetNP());
  Array<bool, SurfaceElementIndex> keep_sel(mesh.GetNSE());
  Array<bool, ElementIndex> keep_el(mesh.GetNE());
  mesh.LineSegments().DeleteAll();

  keep_point = false;
  for (auto pi : points) keep_point[pi] = true;

  auto set_keep = [&](auto& input, auto& keep_array, auto& els) {
    keep_array = false;
    for (auto ind : input) keep_array[ind] = true;

    for (auto ind : Range(els)) {
      bool& keep = keep_array[ind];
      if (keep) continue;

      for (auto pi : mesh[ind].PNums()) keep |= keep_point[pi];

      if (!keep) mesh[ind].Delete();
    }

    for (auto i = 0; i < els.Size(); i++)
      if (els[i].IsDeleted()) {
        els.DeleteElement(i);
        i--;
      }
  };

  set_keep(sels, keep_sel, mesh.SurfaceElements());
  set_keep(els, keep_el, mesh.VolumeElements());
  // mesh.Compress();

  return mesh_ptr;
}

#ifdef NETGEN_DEBUGGING_GUI
DebuggingGUI::DebuggingGUI() { Start(); }
DebuggingGUI::~DebuggingGUI() { Stop(); }

struct PerSocketData {};
using WebSocket = uWS::WebSocket<false, true, PerSocketData>;

void DebuggingGUI::Start() {
  gui_thread = thread([&]() {
    auto a = uWS::App();
    loop = a.getLoop();
    app = &a;

    auto read_file = [](const string& path) {
      std::ifstream in(path);
      std::stringstream buffer;
      buffer << in.rdbuf();
      return buffer.str();
    };

    const string webgui_file = "/home/matthias/src/webgui/dist/webgui.js";
    const string main_file = "/home/matthias/a.html";

    a.get("/", [&](auto* res, auto* req) { res->end(read_file(main_file)); })
        // .any("*", [&](auto* res, auto* req) { cout << "any " << endl; })
        .get("/index.html",
             [&](auto* res, auto* req) { res->end(read_file(main_file)); })
        .get("/webgui.js",
             [&](auto* res, auto* req) { res->end(read_file(webgui_file)); })
        .connect("/*",
                 [&](auto* res, auto* req) {
                   cout << "connect " << req->getUrl() << endl;
                 })

        .listen(7871, [&](auto* token_) {
          token = token_;
          if (token) {
            cout << "Listening on port " << 7871 << endl;
          }
        });

    a.ws<PerSocketData>(
        "/ws",
        {.compression = uWS::CompressOptions(uWS::DEDICATED_COMPRESSOR_4KB |
                                             uWS::DEDICATED_DECOMPRESSOR),
         .maxPayloadLength = 100 * 1024 * 1024,
         .idleTimeout = 16,
         .maxBackpressure = 100 * 1024 * 1024,
         .closeOnBackpressureLimit = false,
         .resetIdleTimeoutOnSend = false,
         .sendPingsAutomatically = true,
         .upgrade = nullptr,
         .open =
             [&](auto* ws) {
               /* Open event here, you may access ws->getUserData() which points
                * to a PerSocketData struct */
               // cout << "open websocket " << endl;
               websockets.insert(ws);
               Send(data);
             },
         .message =
             [](auto* ws, std::string_view message, uWS::OpCode opCode) {
               cout << "got message " << message << endl;
               ws->send(message, opCode, message.length() > 16 * 1024);
             },
         .dropped =
             [](auto* /*ws*/, std::string_view /*message*/,
                uWS::OpCode /*opCode*/) {
               /* A message was dropped due to set maxBackpressure and
                * closeOnBackpressureLimit limit */
               cout << "ws: dropped message" << endl;
             },
         .drain =
             [](auto* /*ws*/) {
               /* Check ws->getBufferedAmount() here */
               cout << "ws: drain" << endl;
             },
         .close =
             [&](auto* ws, int /*code*/, std::string_view /*message*/) {
               /* You may access ws->getUserData() here */
               // cout << "close" << endl;
               websockets.erase(ws);
             }});
    a.run();
    cout << "Exiting GUI thread" << endl;
    loop = nullptr;
    app = nullptr;
  });
}
void DebuggingGUI::Stop() {
  cout << "close socket" << endl;
  for (auto* ws : websockets) ((WebSocket*)ws)->close();
  us_listen_socket_close(0, (struct us_listen_socket_t*)token);
  cout << "join thread" << endl;
  gui_thread.join();
  cout << "joined thread" << endl;
}

void DebuggingGUI::DrawMesh(const Mesh& m) {
  py::gil_scoped_acquire acquire;
  const auto webgui = py::module::import("netgen.webgui");
  const auto dumps = py::module::import("json").attr("dumps");
  mesh = make_shared<Mesh>();
  *mesh = m;
  const auto py_data =
      webgui.attr("Draw")(mesh, py::arg("show") = false).attr("GetData")();
  const string d = py::cast<string>(dumps(py_data));
  data = json::parse(d);
  Send(d);
}

void DebuggingGUI::DrawPoints(FlatArray<Point<3>> position, string name,
                              string color) {
  DrawObject(position, "points", name, color);
}

void DebuggingGUI::DrawLines(FlatArray<Point<3>> position, string name,
                              string color) {
  DrawObject(position, "lines", name, color);
}
void DebuggingGUI::DrawTrigs(FlatArray<Point<3>> position, string name,
                              string color) {
  // Array<Point<3>> pnts;
  // cout << "trigs " << position << endl;
  // for (auto i : Range(position.Size()/3)) {
  //   pnts.Append(position[3*i+0]);
  //   pnts.Append(position[3*i+1]);
  //   pnts.Append(position[3*i+1]);
  //   pnts.Append(position[3*i+2]);
  //   pnts.Append(position[3*i+2]);
  //   pnts.Append(position[3*i+0]);

  // }
  // cout << "pnts size " << pnts.Size() << endl;
  DrawObject(position, "trigs", name, color);
}

void DebuggingGUI::DrawObject(FlatArray<Point<3>> position, string type, string name, string color) {
  json p = json::array_t();
  for (auto pnt : position) {
    p.push_back(pnt[0]);
    p.push_back(pnt[1]);
    p.push_back(pnt[2]);
  }
  json d;
  d["type"] = type;
  d["name"] = name;
  d["color"] = color;
  d["position"] = p;
  data["objects"].push_back(d);
  Send(data);
}

void DebuggingGUI::Send(const string& message) {
  if (loop) {
    ((uWS::Loop*)loop)->defer([=]() {
      for (auto* ws : websockets) {
        ((WebSocket*)ws)->send(message, uWS::OpCode::TEXT);
      }
    });
  }
}

DebuggingGUI debug_gui = DebuggingGUI();

#endif  // NETGEN_DEBUGGING_GUI
}  // namespace netgen
