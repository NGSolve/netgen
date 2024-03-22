#include "meshing.hpp"

#define NETGEN_DEBUGGING_GUI

#ifdef NETGEN_DEBUGGING_GUI
#include "json.hpp"
using json = nlohmann::json;
#endif  // NETGEN_DEBUGGING_GUI

namespace netgen {
unique_ptr<Mesh> GetOpenElements(const Mesh& m, int dom = 0);

unique_ptr<Mesh> FilterMesh(
    const Mesh& m, FlatArray<PointIndex> points,
    FlatArray<SurfaceElementIndex> sels = Array<SurfaceElementIndex>{},
    FlatArray<ElementIndex> els = Array<ElementIndex>{});

#ifdef NETGEN_DEBUGGING_GUI
class DebuggingGUI {
 public:
  DebuggingGUI();
  ~DebuggingGUI();

  void Start();
  void Stop();

  void DrawMesh(const string& name, const Mesh& m);
  void DrawPoints(const string& name, const Mesh& m,
                  FlatArray<PointIndex> points);
  void DrawLines(const string& name, const Mesh& m,
                 FlatArray<SegmentIndex> lines);
  void DrawTrigs(const string& name, const Mesh& m,
                 FlatArray<SurfaceElementIndex> trigs);
  void DrawTets(const string& name, const Mesh& m,
                FlatArray<ElementIndex> tets);
  void AddComponent(const Mesh& m);

 private:
  thread gui_thread;
  void *app, *loop, *token;
  std::set<void*> websockets;
  json data;

  void Send(const json& data) { Send(data.dump()); }
  void Send(const string& message);
};

extern DebuggingGUI debug_gui;
#endif  // NETGEN_DEBUGGING_GUI

}  // namespace netgen
