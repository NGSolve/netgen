#include <meshing.hpp>
#include "writeuser.hpp"

#include <variant>

namespace netgen
{
void ReadMeditFormat (Mesh & mesh, const filesystem::path & filename)
{
  static Timer tall("ReadMeditMesh"); RegionTimer rtall(tall);
  auto fin = ifstream(filename);
  string token;
  int version, dim;
  mesh.ClearFaceDescriptors();

  int index_cnt[4] = {0,0,0,0};
  std::map<int, int>index_map[4];
  auto getIndex = [&](int dim, int index) {
    if(index_map[dim].count(index) == 0) {
      index_map[dim][index] = ++index_cnt[dim];
      if(dim==2) {
        auto fd = FaceDescriptor(index_map[dim][index]-1,1,0,0);
        fd.SetBCProperty(index_map[dim][index]);
        mesh.AddFaceDescriptor (fd);
      }
    }
    return index_map[dim][index];
  };

  while(true) {
    fin >> token;
    int index;
    if(token == "End")
      break;

    if(token == "MeshVersionFormatted") {
      fin >> version;
    }
    if(token == "Dimension") {
      fin >> dim;
      mesh.SetDimension(dim);
    }
    if(token == "Vertices") {
      int nvert;
      fin >> nvert;
      Point<3> p;
      for(auto k : Range(nvert)) {
        for(auto i : Range(dim))
        fin >> p[i];
        fin >> index;
        mesh.AddPoint(p);
      }
    }
    if(token == "Edges") {
      int nedge;
      fin >> nedge;
      Segment seg;
      for(auto k : Range(nedge)) {
        for(auto i : Range(2))
        fin >> seg[i];
        fin >> seg.edgenr;
        seg.edgenr = getIndex(1, seg.edgenr);
        seg.si = seg.edgenr;
        mesh.AddSegment(seg);
      }
    }
    if(token == "Triangles") {
      int ntrig, index;
      fin >> ntrig;
      Element2d sel;
      for(auto k : Range(ntrig)) {
        for(auto i : Range(3))
        fin >> sel[i];
        fin >> index;
        sel.SetIndex(getIndex(2, index));
        mesh.AddSurfaceElement(sel);
      }
    }
    if(token == "Tetrahedra") {
      int ntet;
      fin >> ntet;
      Element el(4);
      for(auto k : Range(ntet)) {
        for(auto i : Range(4))
        fin >> el[i];
        fin >> index;
        el.SetIndex(getIndex(3, index));
        el.Invert();
        mesh.AddVolumeElement(el);
      }
    }
  }
}

void WriteMeditFormat (const Mesh & mesh, const filesystem::path & filename)
{
  static Timer tall("WriteMeditFormat"); RegionTimer rtall(tall);
  auto fout = ofstream(filename);
  fout << "MeshVersionFormatted 2\n";
  fout << "Dimension\n" << mesh.GetDimension() << endl;
  fout << "Vertices\n" << mesh.GetNP() << endl;
  int base_index = 0;
  int max_index = 0;
  auto getIndex = [&](int i) {
    max_index = max(max_index, i+base_index);
    return base_index+i;
  };
  fout << setprecision(16);

  for(const auto & p : mesh.Points())
  {
    for(auto i : Range(mesh.GetDimension()))
      fout << setw(20) << p[i];
    fout << setw(6) << getIndex(1) << endl;
  }

  base_index = max_index;
  fout << "Edges\n" << mesh.GetNSeg() << endl;
  for(const auto & seg : mesh.LineSegments())
  fout << seg[0] << ' ' << seg[1] << ' ' << getIndex(seg.edgenr) << endl;

  base_index = max_index;
  fout << "Triangles\n" << mesh.GetNSE() << endl;
  for(const auto & sel : mesh.SurfaceElements())
  fout << sel[0] << ' ' << sel[1] << ' ' << sel[2] << ' ' << getIndex(sel.GetIndex()) << endl;

  base_index = max_index;
  fout << "Tetrahedra\n" << mesh.GetNE() << endl;
  for(const auto & el : mesh.VolumeElements())
  fout << el[0] << ' ' << el[1] << ' ' << el[2] << ' ' << el[3] << '\t' << getIndex(el.GetIndex()) << endl;

  fout << "End" << endl;
}

static RegisterUserFormat reg_medit ("Medit Format", {".mesh"}, ReadMeditFormat, WriteMeditFormat);
} // namespace netgen
