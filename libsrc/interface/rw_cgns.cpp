#include <meshing.hpp>
#include "writeuser.hpp"

#ifdef NG_CGNS
#include <variant>

#include <cgnslib.h>

namespace netgen::cg
{
  typedef ngcore::ClosedHashTable<ngcore::INT<3,size_t>, size_t> PointTable;

  int getDim(ElementType_t type)
    {
      switch(type)
        {
          case BAR_2:
          case BAR_3:
              return 1;
          case TRI_3:
          case TRI_6:
          case QUAD_4:
          case QUAD_8:
              return 2;
          case TETRA_4:
          case TETRA_10:
          case PYRA_5:
          case PYRA_13:
          case HEXA_8:
          case HEXA_20:
          case PENTA_6:
          case PENTA_15:
              return 3;
          default:
              throw Exception("Read CGNS: unknown element type " + string(cg_ElementTypeName(type)));
        }
    }

  Segment ReadCGNSElement1D( ElementType_t type, FlatArray<cgsize_t> verts )
    {
      int np;
      cg_npe(type, &np);

      Segment s;
      for (auto i : Range(np))
          s[i] = verts[i];
      return s;
    }

  Element2d ReadCGNSElement2D( ElementType_t type, FlatArray<cgsize_t> verts )
    {
//       static constexpr int map_tri3[]  = {0,2,1};
      static constexpr int map_tri6[]  = {0,2,1,3,5,4}; // untested
//       static constexpr int map_quad4[] = {0,3,2,1};
      static constexpr int map_quad8[] = {0,3,2,1,4,7,6,5}; // untested

      const int * map = nullptr;
      switch(type)
        {
          case TRI_3:
//               map = map_tri3;
              break;
          case QUAD_4:
//               map = map_quad4;
              break;
          case TRI_6:
              map = map_tri6;
              break;
          case QUAD_8:
              map = map_quad8;
              break;
          default:
              throw Exception("Read CGNS: unknown element type " + string(cg_ElementTypeName(type)));
        }

      int np;
      cg_npe(type, &np);

      Element2d el(np);
      for (auto i : Range(np))
          el[i] = verts[i];
      return el;
    }

  Element ReadCGNSElement3D( ElementType_t type, FlatArray<cgsize_t> verts )
    {
      static constexpr int map_tet4[]   = {0,2,1,3};
      static constexpr int map_prism6[] = {0,2,1,3,5,4};
      static constexpr int map_pyra5[]  = {0,3,2,1,4};
      static constexpr int map_hexa8[]  = {0,3,2,1,4,7,6,5};
      int np;
      cg_npe(type, &np);

      const int * map = nullptr;
      switch(type)
        {
          case TETRA_4:
              map = map_tet4; break;
          case PYRA_5:
              map = map_pyra5; break;
          case PENTA_6:
              map = map_prism6; break;
          case HEXA_8:
              map = map_hexa8; break;
              // TODO: Second order elements
          case TETRA_10:
          case PYRA_13:
          case HEXA_20:
          case PENTA_15:
          default:
              throw Exception("Read CGNS: unknown element type " + string(cg_ElementTypeName(type)));
        }

      Element el(np);
      for (auto i : Range(np))
          el[i] = verts[map[i]];
      return el;
    }

  void WriteCGNSElement( const Segment & el, Array<cgsize_t> & verts )
    {
      verts.Append(BAR_2);
      verts.Append(el[0]);
      verts.Append(el[1]);
    }

  void WriteCGNSElement( const Element2d & el, Array<cgsize_t> & verts )
    {
      static constexpr int map_tri6[]  = {0,2,1,3,5,4}; // untested
      static constexpr int map_quad8[] = {0,3,2,1,4,7,6,5}; // untested

      ElementType_t type;

      const int * map = nullptr;
      switch(el.GetType())
        {
          case TRIG:
              type = TRI_3;
              break;
          case QUAD:
              type = QUAD_4;
              break;
          case TRIG6:
              type = TRI_6;
              map = map_tri6;
              break;
          case QUAD8:
              type = QUAD_8;
              map = map_quad8;
              break;
              // TODO: Second order elements
          default:
              throw Exception("Write CGNS: unknown element type " + ToString(el.GetType()));
        }

      verts.Append(type);

      for (auto i : Range(el.GetNP()))
          verts.Append(el[i]);
    }

  void WriteCGNSElement( const Element & el, Array<cgsize_t> & verts )
    {
      static constexpr int map_tet4[]   = {0,2,1,3};
      static constexpr int map_prism6[] = {0,2,1,3,5,4};
      static constexpr int map_pyra5[]  = {0,3,2,1,4};
      static constexpr int map_hexa8[]  = {0,3,2,1,4,7,6,5};

      ElementType_t type;

      const int * map = nullptr;
      switch(el.GetType())
        {
          case TET:
              map = map_tet4;
              type = TETRA_4;
              break;
          case PYRAMID:
              type = PYRA_5;
              map = map_pyra5;
              break;
          case PRISM:
              type = PENTA_6;
              map = map_prism6;
              break;
          case HEX:
              type = HEXA_8;
              map = map_hexa8;
              break;
              // TODO: Second order elements
          default:
              throw Exception("Write CGNS: unknown element type " + ToString(el.GetType()));
        }

      verts.Append(type);

      for (auto i : Range(el.GetNP()))
          verts.Append(el[map[i]]);
    }

  int WriteCGNSRegion( const Mesh & mesh, int dim, int index, int fn, int base, int zone, int ne_before )
  {
    int meshdim = mesh.GetDimension();
    int codim = meshdim-dim;

    if(codim < 0 || codim > 2)
      return 0;

    // make sure that each material/boundary name is unique
    string prefix[] = { "dom_", "bnd_", "bbnd_" };
    string name = prefix[meshdim-dim] + ToString(index) + "_";

    if(codim==0) name += mesh.GetMaterial(index+1);
    if(codim==1) name += *mesh.GetBCNamePtr(index);
    if(codim==2) name += mesh.GetCD2Name(index);

    int ne = 0;
    Array<int> data;

    if(dim==3)
      for(const auto el : mesh.VolumeElements())
        if(el.GetIndex()==index)
        {
          ne++;
          WriteCGNSElement(el, data);
        }

    if(dim==2)
      for(const auto el : mesh.SurfaceElements())
        if(el.GetIndex()==index)
        {
          ne++;
          WriteCGNSElement(el, data);
        }

    if(dim==1)
      for(const auto el : mesh.LineSegments())
        if(el.si==index)
        {
          ne++;
          WriteCGNSElement(el, data);
        }

    if(ne==0)
      return 0;

    int section;
    int start = 1;
    int end = ne;
#if CGNS_VERSION < 3400
    cg_section_write(fn,base,zone, name.c_str(), MIXED, ne_before+1, ne_before+ne, 0, &data[0], &section);
#else
    cg_poly_section_write(fn,base,zone, name.c_str(), MIXED, ne_before+1, ne_before+ne, 0, &data[0], nullptr, &section);
#endif

    return ne;
  }

  // maps cgns node type to ngsolve node type
  // enum NODE_TYPE { NT_VERTEX = 0, NT_EDGE = 1, NT_FACE = 2, NT_CELL = 3, NT_ELEMENT = 4, NT_FACET = 5 };
  int getNodeType( GridLocation_t location )
    {
      switch(location)
        {
          case Vertex:
              return 0;
          case CellCenter:
              return 3;
          case FaceCenter:
              return 2;
          case EdgeCenter:
              return 1;
          default:
              throw Exception("Read CGNS: unknown grid location " + string(cg_GridLocationName(location)));
        }
    }

  GridLocation_t getCGNodeType( int node_type )
  {
    switch(node_type)
    {
      case 0:
        return Vertex;
      case 1:
        return EdgeCenter;
      case 2:
        return FaceCenter;
      case 3:
        return CellCenter;
      default:
        throw Exception("Write CGNS: unknown node type " + ToString(node_type));
    }
  }


  struct Solution
    {
      int fn, base, zone, solution;
      string name;
      GridLocation_t location; // solution is defined on either cells, faces, edges or vertices
      PointSetType_t  point_type;
      cgsize_t n_points;

      Array<string> field_names;
      Array<DataType_t> field_datatypes;

      Solution() = default;

      Solution(int fn_, int base_, int zone_, int solution_)
          : fn(fn_), base(base_), zone(zone_), solution(solution_)
        {
          char solname[100];
          cg_sol_info(fn, base, zone, solution, solname, &location);
          name = solname;
          cg_sol_ptset_info(fn, base, zone, solution, &point_type, &n_points);

          int n_fields = 0;
          cg_nfields(fn, base, zone, solution, &n_fields);

          field_names.SetSize(n_fields);
          field_datatypes.SetSize(n_fields);
          for(auto fi : Range(n_fields))
            {
              char buf[100];
              cg_field_info(fn, base, zone, solution, fi+1, &field_datatypes[fi], buf);
              field_names[fi] = buf;
            }
        }
    };

  struct Zone
    {
      ZoneType_t zone_type;
      int fn, base, zone;
      int first_index_1d, first_index_2d, first_index_3d;
      int nv=0, ne_1d=0, ne_2d=0, ne_3d=0;

      Array<string> names_1d, names_2d, names_3d;

      string name;
      cgsize_t size[3];

      Array<Solution> solutions;

      Zone(int fn_, int base_, int zone_)
          : fn(fn_), base(base_), zone(zone_)
        {
          cg_zone_type(fn, base, zone, &zone_type);
          char zone_name[100];
          cg_zone_read(fn,base,zone, zone_name, size);
          nv = size[0];

          int n_solutions;
          cg_nsols(fn, base, zone, &n_solutions);

          solutions.SetSize(n_solutions);
          for(auto si : Range(n_solutions))
              solutions[si] = Solution{fn, base, zone, si+1};
        }

      void ReadSolutions( int meshdim, std::vector<string> & sol_names, std::vector<Array<double>> & sol_values, std::vector<int> & sol_locations )
        {
          static Timer tall("CGNS::ReadSolutions"); RegionTimer rtall(tall);
          for (auto & sol : solutions)
            {
              for (auto fi : Range(sol.field_names.Size()))
                {
                  cgsize_t size = sol.n_points;
                  size=0; // TODO: check if sol.point_type is a list or range, and handle appropriately
                  if(size==0)
                    {
                      switch(sol.location)
                        {
                          case Vertex:
                              size = nv;
                              break;
                          case CellCenter:
                              size = (meshdim == 3 ? ne_3d : ne_2d);
                              break;
                          case FaceCenter:
                          case IFaceCenter:
                          case JFaceCenter:
                          case KFaceCenter:
                          case EdgeCenter:
                          default:
                              throw Exception("Read CGNS: unknown grid location " + string(cg_GridLocationName(sol.location)));
                        }
                    }

                  auto values = Array<double>(size);

                  cgsize_t imin = 1UL;
                  cg_field_read(fn, base, zone, sol.solution, sol.field_names[fi].c_str(), RealDouble, &imin, &size, &values[0]);
                  sol_names.push_back(sol.field_names[fi]);
                  sol_values.emplace_back(std::move(values));
                  sol_locations.push_back(getNodeType(sol.location));
                }
            }
        }

      void ReadMesh( Mesh & mesh, PointTable & point_table )
        {
          static Timer tall("CGNS::ReadMesh-Zone"); RegionTimer rtall(tall);
          static Timer tsection("CGNS::ReadMesh-Section");
          first_index_1d = mesh.GetRegionNamesCD(2).Size();
          first_index_2d = mesh.GetRegionNamesCD(1).Size();
          first_index_3d = mesh.GetRegionNamesCD(0).Size();

          Array<double> x(nv), y(nv), z(nv);
          cgsize_t imin=1;
          cg_coord_read(fn,base,zone, "CoordinateX", RealDouble, &imin, &nv, &x[0]);
          cg_coord_read(fn,base,zone, "CoordinateY", RealDouble, &imin, &nv, &y[0]);
          cg_coord_read(fn,base,zone, "CoordinateZ", RealDouble, &imin, &nv, &z[0]);

          Array<cgsize_t> point_map(nv);

          for(auto i : Range(nv))
          {
            ngcore::INT<3,size_t> hash = {*reinterpret_cast<size_t*>(&x[i]), *reinterpret_cast<size_t*>(&y[i]), *reinterpret_cast<size_t*>(&z[i])};
            size_t pi_ng;
            size_t pos;
            // check if this point is new
            if( point_table.PositionCreate (hash, pos) )
            {
              pi_ng = mesh.AddPoint( {x[i], y[i], z[i]} );
              point_table.SetData(pos, pi_ng);
            }
            else
              point_table.GetData(pos, pi_ng);

            point_map[i] = pi_ng;
          }

          int nsections;
          cg_nsections(fn, base, zone, &nsections);

          int index_1d = first_index_1d;
          int index_2d = first_index_2d;
          int index_3d = first_index_3d;

          for (auto section : Range(1,nsections+1))
            {
              RegionTimer rtsection(tsection);
              char sec_name[100];
              ElementType_t type;
              cgsize_t start, end;
              int nbndry, parent_flag;

              cg_section_read(fn, base, zone, section, sec_name,  &type, &start, &end, &nbndry, &parent_flag);
              PrintMessage(4, "Read section ", section, " with name ", sec_name, " and element type ", cg_ElementTypeName(type));

              string ngname{sec_name};

              for (char & c : ngname)
                  if(c==' ')
                      c = '_';


              if(type==MIXED)
                {
                  bool have_1d_elements = false;
                  bool have_2d_elements = false;
                  bool have_3d_elements = false;

                  cgsize_t nv;
                  cg_ElementDataSize(fn, base, zone, section, &nv);

                  Array<cgsize_t> vertices(nv);
#if CGNS_VERSION < 3400
                  cg_elements_read(fn, base, zone, section, &vertices[0], nullptr);
#else
                  cg_poly_elements_read(fn, base, zone, section, &vertices[0], nullptr, nullptr);
#endif

                  size_t vi = 0;
                  while(vi<nv)
                    {
                      auto type = static_cast<ElementType_t>(vertices[vi++]);
                      int dim = getDim(type);

                      int np;
                      cg_npe(type, &np);

                      for (auto & v : vertices.Range(vi, vi+np))
                        v = point_map[v-1];

                      if(dim==1)
                        {
                          if(!have_1d_elements)
                          {
                            index_1d++;
                            have_1d_elements = true;
                            mesh.AddEdgeDescriptor(EdgeDescriptor{});
                            names_1d.Append(ngname);
                          }
                          auto el = ReadCGNSElement1D(type, vertices.Range(vi, vertices.Size()));
                          el.si = index_1d;
                          mesh.AddSegment(el);
                          vi += el.GetNP();
                          ne_1d++;
                        }

                      if(dim==2)
                        {
                          if(!have_2d_elements)
                          {
                            index_2d++;
                            have_2d_elements = true;
                            mesh.AddFaceDescriptor(FaceDescriptor(index_2d, 1, 0, 1));
                            names_2d.Append(ngname);
                          }
                          auto el = ReadCGNSElement2D(type, vertices.Range(vi, vertices.Size()));
                          el.SetIndex(index_2d);
                          mesh.AddSurfaceElement(el);
                          vi += el.GetNP();
                          ne_2d++;
                        }

                      if(dim==3)
                        {
                          if(!have_3d_elements)
                          {
                            index_3d++;
                            have_3d_elements = true;
                            names_3d.Append(ngname);
                          }

                          auto el = ReadCGNSElement3D(type, vertices.Range(vi, vertices.Size()));
                          el.SetIndex(index_3d);
                          mesh.AddVolumeElement(el);
                          vi += el.GetNP();
                          ne_3d++;
                        }
                    }
                }
              else
                {
                  int dim = getDim(type);

                  cgsize_t nv;
                  cg_ElementDataSize(fn, base, zone, section, &nv);
                  int np=0;
                  cg_npe(type, &np);

                  Array<cgsize_t> vertices(nv);
                  cg_elements_read(fn, base, zone, section, &vertices[0], nullptr);
                  for (auto & v : vertices)
                    v = point_map[v-1];
                  int ne_section = nv/np;

                  if(dim==1)
                    {
                      index_1d++;
                      mesh.AddEdgeDescriptor(EdgeDescriptor{});
                      names_1d.Append(ngname);
                      for(auto i : Range(ne_section))
                        {
                          auto el = ReadCGNSElement1D(type, vertices.Range(np*i, np*(i+1)));
                          el.si = index_1d;
                          mesh.AddSegment(el);
                        }
                      ne_1d += ne_section;
                    }

                  if(dim==2)
                    {
                      index_2d++;
                      mesh.AddFaceDescriptor(FaceDescriptor(index_2d, 1, 0, 1));
                      names_2d.Append(ngname);
                      for(auto i : Range(ne_section))
                        {
                          auto el = ReadCGNSElement2D(type, vertices.Range(np*i, np*(i+1)));
                          el.SetIndex(index_2d);
                          mesh.AddSurfaceElement(el);
                        }
                      ne_2d += ne_section;
                    }

                  if(dim==3)
                    {
                      index_3d++;
                      names_3d.Append(ngname);
                      for(auto i : Range(ne_section))
                        {
                          auto el = ReadCGNSElement3D(type, vertices.Range(np*i, np*(i+1)));
                          el.SetIndex(index_3d);
                          mesh.AddVolumeElement(el);
                        }
                      ne_3d += ne_section;
                    }
                }
            }

          mesh.GetRegionNamesCD(2).SetSize(index_1d);
          mesh.GetRegionNamesCD(1).SetSize(index_2d);
          mesh.GetRegionNamesCD(0).SetSize(index_3d);
          mesh.GetRegionNamesCD(2) = nullptr;
          mesh.GetRegionNamesCD(1) = nullptr;
          mesh.GetRegionNamesCD(0) = nullptr;
        }

      void SetNames( Mesh & mesh )
        {
          if(mesh.GetDimension() == 2)
          {
            for (auto i : Range(names_1d.Size()))
              mesh.SetBCName(first_index_1d + i, names_1d[i]);

            for (auto i : Range(names_2d.Size()))
              mesh.SetMaterial(first_index_2d + i +1, names_2d[i]);
          }
          else
          {
            for (auto i : Range(names_1d.Size()))
              mesh.SetCD2Name(first_index_1d + i +1, names_1d[i]);

            for (auto i : Range(names_2d.Size()))
            {
              mesh.SetBCName(first_index_2d + i, names_2d[i]);
              mesh.GetFaceDescriptor(first_index_2d + i +1).SetDomainIn(first_index_3d+1);
            }

            for (auto i : Range(names_3d.Size()))
              mesh.SetMaterial(first_index_3d + i +1, names_3d[i]);
          }
        }
    };
}

namespace netgen
{
  int ReadCGNSMesh (Mesh & mesh, const string & filename, Array<unique_ptr<cg::Zone>> & zones)
    {
      mesh.SetDimension(3);
      static Timer tall("CGNS::ReadMesh"); RegionTimer rtall(tall);
      int fn;
      cg_open(filename.c_str(),CG_MODE_READ,&fn);

      int base = 1;
      int nzones;
      cg_nzones(fn, base, &nzones);

      int n_vertices = 0;
      for (auto zi : Range(1, nzones+1))
      {
        int size[3];
        char name[100];
        cg_zone_read(fn,base,zi, name, size);
        n_vertices += size[0];
      }

      cg::PointTable points(2*n_vertices);

      for (auto zi : Range(1, nzones+1))
        {
          ZoneType_t zone_type;
          cg_zone_type(fn, base, zi, &zone_type);
          if(zone_type != Unstructured )
            {
              PrintMessage(2, "skipping zone with type ", cg_ZoneTypeName(zone_type) );
              continue;
            }
          auto zone = make_unique<cg::Zone>(fn, base, zi);
          zone->ReadMesh( mesh, points );
          zones.Append(std::move(zone));
        }

      if(mesh.GetNE() == 0)
        mesh.SetDimension(2);

      for (auto & zone : zones)
        zone->SetNames(mesh);

      mesh.UpdateTopology();
      const auto & topo = mesh.GetTopology();

      for (auto sei : Range(mesh.SurfaceElements()))
      {
        int ei0, ei1;
        topo.GetSurface2VolumeElement (sei+1, ei0, ei1);
        auto si = mesh.SurfaceElement(sei).GetIndex();
        auto & fd = mesh.GetFaceDescriptor(si);

        if(ei0>0)
        {
          int i0 = mesh.VolumeElement(ei0).GetIndex();
          if(fd.DomainIn()!=i0)
            fd.SetDomainOut(i0);
        }

        if(ei1>0)
        {
          int i1 = mesh.VolumeElement(ei1).GetIndex();
          if(fd.DomainIn()!=i1)
            fd.SetDomainOut(i1);
        }
      }
      return fn;
    }

  void ReadCGNSMesh (Mesh & mesh, const string & filename)
    {
      Array<unique_ptr<cg::Zone>> zones;
      int fn = ReadCGNSMesh(mesh, filename, zones);
      cg_close(fn);
    }

  // Reads mesh and solutions of .csns file
  tuple<shared_ptr<Mesh>, vector<string>, vector<Array<double>>, vector<int>> ReadCGNSFile(string filename, int base)
    {
      static Timer tall("CGNS::ReadFile"); RegionTimer rtall(tall);

      auto mesh = make_shared<Mesh>();
      Array<unique_ptr<cg::Zone>> zones;
      int fn = ReadCGNSMesh(*mesh, filename, zones);

      std::vector<string> names;
      std::vector<Array<double>> values;
      std::vector<int> locations;

      for (auto & zone : zones)
        zone->ReadSolutions( mesh->GetDimension(), names, values, locations );

      cg_close(fn);
      return std::make_tuple(mesh, names, values, locations);
    }

  void WriteCGNSMesh (const Mesh & mesh, int fn, int & base, int & zone)
    {
      int dim = mesh.GetDimension();
      cg_base_write(fn, "mesh", dim, dim, &base);

      int nv = static_cast<int>(mesh.GetNV());
      int ne = mesh.GetNE();

      Array<double> x, y, z;
      for(auto & p : mesh.Points())
      {
        x.Append(p[0]);
        y.Append(p[1]);
        z.Append(p[2]);
      }

      cgsize_t isize[3] = { nv, ne, 0 };
      cg_zone_write(fn,base, "mesh", isize, Unstructured, &zone);

      int coord;
      cg_coord_write(fn,base,zone, RealDouble, "CoordinateX", &x[0], &coord);
      cg_coord_write(fn,base,zone, RealDouble, "CoordinateY", &y[0], &coord);
      cg_coord_write(fn,base,zone, RealDouble, "CoordinateZ", &z[0], &coord);

      int imax3 = 0;
      for(const auto & el : mesh.VolumeElements())
        imax3 = max(imax3, el.GetIndex());

      int imax2 = 0;
      for(const auto & el : mesh.SurfaceElements())
        imax2 = max(imax2, el.GetIndex());

      int imax1 = 0;
      for(const auto & el : mesh.LineSegments())
        imax1 = max(imax1, el.si);

      int ne_written = 0;
      int meshdim = mesh.GetDimension();

      for(const auto i : IntRange(imax3))
        ne_written += cg::WriteCGNSRegion(mesh, 3, i+1, fn, base, zone, ne_written);

      for(const auto i : IntRange(imax2))
        ne_written += cg::WriteCGNSRegion(mesh, 2, i+1, fn, base, zone, ne_written);

      for(const auto i : IntRange(imax1))
        ne_written += cg::WriteCGNSRegion(mesh, 1, i+1, fn, base, zone, ne_written);

    }

  void WriteCGNSMesh (const Mesh & mesh, const string & filename)
    {
      static Timer tall("CGNS::WriteMesh"); RegionTimer rtall(tall);
      int fn, base, zone;
      cg_open(filename.c_str(),CG_MODE_WRITE,&fn);

      WriteCGNSMesh(mesh, fn, base, zone);

      cg_close(fn);
    }


  void WriteCGNSFile(shared_ptr<Mesh> mesh, string filename, vector<string> fields, vector<Array<double>> values, vector<int> locations)
  {
      static Timer tall("CGNS::WriteFile"); RegionTimer rtall(tall);
      int fn, base, zone;
      cg_open(filename.c_str(),CG_MODE_WRITE,&fn);

      WriteCGNSMesh(*mesh, fn, base, zone);

      for(auto i : IntRange(fields.size()))
      {
        int section, field;
        string name = "solution_" + ToString(i);


        cg_sol_write(fn, base, zone, name.c_str(), cg::getCGNodeType(locations[i]), &section);
        cg_field_write(fn, base, zone, section, RealDouble, fields[i].c_str(), &values[i][0], &field);
      }

      cg_close(fn);
  }

}

#else // NG_CGNS

namespace netgen
{
  void ReadCGNSMesh (Mesh & mesh, const string & filename)
    {
      PrintMessage(1, "Could not import CGNS mesh: Netgen was built without CGNS support");
    }

  tuple<shared_ptr<Mesh>, vector<string>, vector<Array<double>>, vector<int>> ReadCGNSFile(string filename, int base)
    {
      throw Exception("Netgen was built without CGNS support");
    }

  void WriteCGNSMesh (const Mesh & mesh, const string & filename)
    {
      PrintMessage(1, "Could not write CGNS mesh: Netgen was built without CGNS support");
    }

  void WriteCGNSFile(shared_ptr<Mesh> mesh, string filename, vector<string> fields, vector<Array<double>> values, vector<int> locations)
    {
      throw Exception("Netgen was built without CGNS support");
    }

}

#endif // NG_CGNS
