#include <regex>

#include <meshing.hpp>
#include "rw_medit.hpp"

namespace netgen
{
void ReadMeditFormat (Mesh & mesh, const filesystem::path & filename, map<tuple<int,int>, int> & index_map)
{
  static Timer tall("ReadMeditMesh"); RegionTimer rtall(tall);
  if(!filesystem::exists(filename))
    throw Exception("File does not exist: " + filename.string());
  auto fin = ifstream(filename);
  string token;
  int version, dim;
  mesh.ClearFaceDescriptors();

  int index_cnt[4] = {0,0,0,0};
  auto getIndex = [&](int eldim, int index) {
    if(index_map.count({eldim,index})==0) {
      auto n = ++index_cnt[eldim];
      index_map[{eldim, index}] = n;
      if(eldim==2) {
        auto fd = FaceDescriptor(n-1,1,0,0);
        fd.SetBCProperty(n);
        mesh.AddFaceDescriptor (fd);
      }
    }
    return index_map[{eldim, index}];
  };

  while(true) {
    fin >> token;
    int index;
    // cout << "token: " << token << endl;
    if(token == "End") {
      break;
    }
    else if(token == "" || std::regex_match(token, std::regex("^[\\s]*$"))) {
      continue;
    }
    else if(token == "MeshVersionFormatted") {
      fin >> version;
    }
    else if(token == "Dimension") {
      fin >> dim;
      mesh.SetDimension(dim);
    }
    else if(token == "Vertices") {
      int nvert;
      fin >> nvert;
      Point<3> p{0.,0.,0.};
      for([[maybe_unused]] auto k : Range(nvert)) {
        for(auto i : Range(dim))
          fin >> p[i];
          fin >> index;
          mesh.AddPoint(p);
      }
    }
    else if(token == "Edges") {
      int nedge;
      fin >> nedge;
      Segment seg;
      for([[maybe_unused]] auto k : Range(nedge)) {
        for(auto i : Range(2))
          fin >> seg[i];
          fin >> seg.edgenr;
          seg.edgenr = getIndex(1, seg.edgenr);
          seg.si = seg.edgenr;
          mesh.AddSegment(seg);
      }
    }
    else if(token == "Triangles") {
      int ntrig, index;
      fin >> ntrig;
      Element2d sel;
      for([[maybe_unused]] auto k : Range(ntrig)) {
        for(auto i : Range(3))
          fin >> sel[i];
          fin >> index;
          sel.SetIndex(getIndex(2, index));
          mesh.AddSurfaceElement(sel);
      }
    }
    else if(token == "Tetrahedra") {
      int ntet;
      fin >> ntet;
      Element el(4);
      for([[maybe_unused]] auto k : Range(ntet)) {
        for(auto i : Range(4))
          fin >> el[i];
          fin >> index;
          el.SetIndex(getIndex(3, index));
          el.Invert();
          mesh.AddVolumeElement(el);
      }
    }
    else if(token == "Corners") {
      int ncorners;
      fin >> ncorners;
      Element0d el;
      for([[maybe_unused]] auto k : Range(ncorners)) {
          fin >> el.pnum;
      }
    }
    else if(token == "RequiredVertices") {
      int nverts;
      fin >> nverts;
      int vert;
      for([[maybe_unused]] auto k : Range(nverts)) {
          fin >> vert;
      }
    }
    else if(token == "Normals") {
      int nnormals;
      fin >> nnormals;
      Vec<3> normal;
      for([[maybe_unused]] auto k : Range(nnormals)) {
          fin >> normal[0];
          fin >> normal[1];
          fin >> normal[2];
      }
    }
    else if(token == "NormalAtVertices") {
      int nnormals;
      fin >> nnormals;
      int vert;
      int normal;
      for([[maybe_unused]] auto k : Range(nnormals)) {
        fin >> normal;
        fin >> vert;
      }
    }
    else if(token == "Tangents") {
      int ntangents;
      fin >> ntangents;
      Vec<3> tangent;
      for([[maybe_unused]] auto k : Range(ntangents)) {
        fin >> tangent[0];
        fin >> tangent[1];
        fin >> tangent[2];
      }
    }
    else if(token == "TangentAtVertices") {
      int ntangents;
      fin >> ntangents;
      int vert;
      int tangent;
      for([[maybe_unused]] auto k : Range(ntangents)) {
        fin >> tangent;
        fin >> vert;
      }
    }
    else if(token == "Ridges") {
      int nridges;
      fin >> nridges;
      int ridge;
      for([[maybe_unused]] auto k : Range(nridges)) {
        fin >> ridge;
      }
    }
    else {
      cout << "unknown token " << token << endl;
      int nitems;
      fin >> nitems;
      string s;
      for([[maybe_unused]] auto i : Range(nitems))
        fin >> s; // read one line
    }
  }
}

void ReadMeditFormat (Mesh & mesh, const filesystem::path & filename)
{
  map<tuple<int, int>, int> index_map;
  ReadMeditFormat(mesh, filename, index_map);
}


void WriteMeditFormat (const Mesh & mesh, const filesystem::path & filename, map<tuple<int,int>, int> & index_map)
{
  static Timer tall("WriteMeditFormat"); RegionTimer rtall(tall);
  auto fout = ofstream(filename);
  fout << "MeshVersionFormatted 2\n";
  fout << "Dimension\n" << mesh.GetDimension() << endl;
  fout << "Vertices\n" << mesh.GetNP() << endl;
  int base_index = 0;
  int max_index = 0;
  auto getIndex = [&](int i, int dim) {
    max_index = max(max_index, i+base_index);
    auto index = base_index+i;
    index_map[{dim,i}] = index;
    return index;
  };
  fout << setprecision(16);

  for(const auto & p : mesh.Points())
  {
    for(auto i : Range(mesh.GetDimension()))
      fout << p[i] << ' ';
    fout << getIndex(1, 0) << endl;
  }

  base_index = max_index;
  fout << "Edges\n" << mesh.GetNSeg() << endl;
  for(const auto & seg : mesh.LineSegments())
    fout << seg[0] << ' ' << seg[1] << ' ' << getIndex(seg.edgenr, 1) << endl;

  base_index = max_index;
  fout << "Triangles\n" << mesh.GetNSE() << endl;
  for(const auto & sel : mesh.SurfaceElements())
    fout << sel[0] << ' ' << sel[1] << ' ' << sel[2] << ' ' << getIndex(sel.GetIndex(), 2) << endl;

  base_index = max_index;
  fout << "Tetrahedra\n" << mesh.GetNE() << endl;
  for(const auto & el : mesh.VolumeElements())
    fout << el[0] << ' ' << el[1] << ' ' << el[2] << ' ' << el[3] << '\t' << getIndex(el.GetIndex(), 3) << endl;

  fout << "End" << endl;
}

void WriteMeditFormat (const Mesh & mesh, const filesystem::path & filename)
{
  map<tuple<int,int>, int> index_map;
  WriteMeditFormat(mesh, filename, index_map);
}

static RegisterUserFormat reg_medit ("Medit Format", {".mesh"},
                                     static_cast<void(*)(Mesh &, const filesystem::path&)>(ReadMeditFormat),
                                     static_cast<void(*)(const Mesh &, const filesystem::path&)>(WriteMeditFormat));
} // namespace netgen
