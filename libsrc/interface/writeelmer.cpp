
//
//  Write Elmer file
//
//

#include <mystdlib.h>

#include <myadt.hpp>
#include <linalg.hpp>
#include <csg.hpp>
#include <meshing.hpp>
#include <sys/stat.h>

#include "writeuser.hpp"


namespace netgen
{
  extern MeshingParameters mparam;


void WriteElmerFormat (const Mesh &mesh,
                       const filesystem::path &dirname)
{
  cout << "write elmer mesh files" << endl;

  std::map<ELEMENT_TYPE, int> tmap;
  tmap[TRIG] = 303;
  tmap[TRIG6] = 306;
  tmap[QUAD] = 404;
  tmap[QUAD8] = 408;
  tmap[TET] = 504;
  tmap[TET10] = 510;
  tmap[PYRAMID] = 605;
  tmap[PYRAMID13] = 613;
  tmap[PRISM] = 706;
  tmap[PRISM15] = 715;
  tmap[HEX] = 808;
  tmap[HEX20] = 820;

  std::map<int, Array<int,int>> pmap;
  pmap[TRIG]  = {1,2,3};
  pmap[TRIG6] = {1,2,3, 6,4,5};
  pmap[QUAD]  = {1,2,3,4};
  pmap[QUAD8] = {1,2,3,4, 5,8,6,7};
  pmap[TET]   = {1,2,3,4};
  pmap[TET10] = {1,2,3,4, 5,8,6,7,9,10};
  pmap[PYRAMID]={1,2,3,4,5};
  pmap[PYRAMID13]= {1,2,3,4,5,6,7,8,9,10,11,12,13};
  pmap[PRISM] = {1,2,3,4,5,6};
  pmap[PRISM15] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
  pmap[HEX]   = {1,2,3,4,5,6,7,8};
  pmap[HEX20] = {1,2,3,4,5,8,6,7,8, 9,12,10,11, 17,20,19,18, 13,16,14,15};

  int np = mesh.GetNP();
  int ne = mesh.GetNE();
  int nse = mesh.GetNSE();
  int i, j;
  // char str[200];
  
  int inverttets = mparam.inverttets;
  int invertsurf = mparam.inverttrigs;

  filesystem::create_directories(dirname);

  auto get_name = [&dirname]( string s ) {
      return filesystem::path(dirname).append(s);
  };

  ofstream outfile_h(get_name("mesh.header"));
  ofstream outfile_n(get_name("mesh.nodes"));
  ofstream outfile_e(get_name("mesh.elements"));
  ofstream outfile_b(get_name("mesh.boundary"));
  ofstream outfile_names(get_name("mesh.names"));

  for( auto codim : IntRange(0, mesh.GetDimension()-1) )
  {
    auto & names = const_cast<Mesh&>(mesh).GetRegionNamesCD(codim);

    for (auto i0 : Range(names) )
    {
      if(names[i0] == nullptr)
        continue;
      string name = *names[i0];
      if(name == "" || name == "default")
        continue;
      outfile_names << "$" << name << "=" << i0+1 << "\n";
    }
  }

  auto get3FacePoints = [](const Element2d & el)
  {
      INDEX_3 i3;
      INDEX_4 i4;
      auto eltype = el.GetType();
      switch (eltype)
      {
          case TRIG:
          case TRIG6:
              i3 = {el[0], el[1], el[2]};
              i3.Sort();
              break;
          case QUAD:
          case QUAD8:
              i4 = {el[0], el[1], el[2], el[3]};
              i4.Sort();
              i3 = {i4[0], i4[1], i4[2]};
              break;
          default:
             throw Exception("Got invalid type (no face)");
      }
      return i3;
  };

  // fill hashtable

  // use lowest three point numbers of lowest-order face to index faces
  INDEX_3_HASHTABLE<int> face2volelement(ne);

  for (int i = 1; i <= ne; i++)
    {
      const Element & el = mesh.VolumeElement(i);

      // getface not working for second order elements -> reconstruct linear element here
      Element linear_el = el;
      linear_el.SetNP(el.GetNV()); // GetNV returns 8 for HEX20 for instance

      for (auto j : Range(1,el.GetNFaces()+1))
	{
          Element2d face;
          linear_el.GetFace(j, face);
	  face2volelement.Set (get3FacePoints(face), i);
          cout << "set " << get3FacePoints(face) << "\tto " << i << endl;
	}
    }

//  outfile.precision(6);
//  outfile.setf (ios::fixed, ios::floatfield);
//  outfile.setf (ios::showpoint);
  
  std::map<ELEMENT_TYPE, size_t> elcount;
  
  for (i = 1; i <= np; i++)
    {
      const Point3d & p = mesh.Point(i);
      
      outfile_n << i << " -1 ";
      outfile_n << p.X() << " ";
      outfile_n << p.Y() << " ";
      outfile_n << p.Z() << "\n";
    }

  for (i = 1; i <= ne; i++)
    {
      Element el = mesh.VolumeElement(i);
      if (inverttets) el.Invert();
      auto eltype = el.GetType();
      elcount[eltype]++;
      outfile_e << i << " " << el.GetIndex() << " " << tmap[eltype] <<  "  ";

      auto & map = pmap[eltype];
      for (j = 1; j <= el.GetNP(); j++)
	{
	  outfile_e << " ";
	  outfile_e << el.PNum(map[j-1]);
	}
      outfile_e << "\n";
    }

  for (i = 1; i <= nse; i++)
    {
      Element2d el = mesh.SurfaceElement(i);
      if (invertsurf) el.Invert();
      auto eltype = el.GetType();
      elcount[eltype]++;
            
      int elind = face2volelement.Get(get3FacePoints(el));
      cout << "get " << get3FacePoints(el) << "\t " << elind << endl;

      outfile_b << i << " " << mesh.GetFaceDescriptor(el.GetIndex()).BCProperty() << 
         " " << elind << " 0 "  << tmap[eltype] << "    ";

      auto & map = pmap[el.GetType()];
      for (j = 1; j <= el.GetNP(); j++)
	{
	  outfile_b << " ";
	  outfile_b << el.PNum(map[j-1]);
	}
      outfile_b << "\n";
    }

  outfile_h << np << " " << ne << " " << nse << "\n";
  outfile_h << "2"     << "\n";

  for( auto & [eltype,count] : elcount )
      outfile_h << tmap[eltype] << " " << count << "\n";
}

}
