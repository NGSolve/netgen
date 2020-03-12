#include <meshing.hpp>
#include "writeuser.hpp"

#ifdef NG_CGNS
#include <variant>

#include <cgnslib.h>

namespace netgen::cg
{
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

  Segment ReadCGNSElement1D( ElementType_t type, FlatArray<cgsize_t> verts, int vert_offset=0 )
    {
      int np;
      cg_npe(type, &np);

      Segment s;
      for (auto i : Range(np))
          s[i] = vert_offset+verts[i];
      return s;
    }

  Element2d ReadCGNSElement2D( ElementType_t type, FlatArray<cgsize_t> verts, int vert_offset=0 )
    {
      static constexpr int map_tri3[]  = {0,2,1};
      static constexpr int map_tri6[]  = {0,2,1,3,5,4}; // untested
      static constexpr int map_quad4[] = {0,3,2,1};
      static constexpr int map_quad8[] = {0,3,2,1,4,7,6,5}; // untested

      const int * map = nullptr;
      switch(type)
        {
          case TRI_3:
              map = map_tri3;
              break;
          case QUAD_4:
              map = map_quad4;
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
          el[i] = vert_offset+verts[map[i]];
      return el;
    }

  Element ReadCGNSElement3D( ElementType_t type, FlatArray<cgsize_t> verts, int vert_offset=0 )
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
          el[i] = vert_offset+verts[map[i]];
      return el;
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
      int nv, ne, first_vertex, first_mat, first_bc;
      Array<int> materials;
      Array<int> boundaries;
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

      void ReadSolutions( std::vector<string> & sol_names, std::vector<Array<double>> & sol_values, std::vector<int> & sol_locations )
        {
          static Timer tall("CGNS::ReadSolutions"); RegionTimer rtall(tall);
          for (auto & sol : solutions)
            {
              for (auto fi : Range(sol.field_names.Size()))
                {
                  cgsize_t size = sol.n_points;
                  if(size==0)
                    {
                      switch(sol.location)
                        {
                          case Vertex:
                              size = nv;
                              break;
                          case CellCenter:
                              size = ne;
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

                  size = size==0 ? nv : size;
                  auto values = Array<double>(size);

                  cgsize_t imin = 1UL;
                  cg_field_read(fn, base, zone, sol.solution, sol.field_names[fi].c_str(), RealDouble, &imin, &size, &values[0]);
                  sol_names.push_back(sol.field_names[fi]);
                  sol_values.emplace_back(std::move(values));
                  sol_locations.push_back(sol.location);
                }
            }
        }

      void ReadMesh( Mesh & mesh )
        {
          static Timer tall("CGNS::ReadMesh-Zone"); RegionTimer rtall(tall);
          static Timer tsection("CGNS::ReadMesh-Section");
          first_vertex = mesh.GetNP();
          first_mat = mesh.GetRegionNamesCD(0).Size();
          first_bc = mesh.GetRegionNamesCD(1).Size();
          ne = 0;

          Array<float> x(nv), y(nv), z(nv);
          cgsize_t imin=1;
          cg_coord_read(fn,base,zone, "CoordinateX", RealSingle, &imin, &nv, &x[0]);
          cg_coord_read(fn,base,zone, "CoordinateY", RealSingle, &imin, &nv, &y[0]);
          cg_coord_read(fn,base,zone, "CoordinateZ", RealSingle, &imin, &nv, &z[0]);

          for(auto i : Range(nv))
              mesh.AddPoint( {x[i], y[i], z[i]} );

          int nsections;
          cg_nsections(fn, base, zone, &nsections);

          int bc = first_bc;
          int material = first_mat;
          for (auto section : Range(1,nsections+1))
            {
              RegionTimer rtsection(tsection);
              char name[100];
              ElementType_t type;
              cgsize_t start, end;
              int nbndry, parent_flag;

              cg_section_read(fn, base, zone, section, name,  &type, &start, &end, &nbndry, &parent_flag);
              PrintMessage(4, "Read section ", section, " with name ", name, " and element type ", cg_ElementTypeName(type));

              string ngname{name};

              for (char & c : ngname)
                  if(c==' ')
                      c = '_';


              if(type==MIXED)
                {
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

                      if(dim==1)
                        {
                          auto el = ReadCGNSElement1D(type, vertices.Range(vi, vertices.Size()), first_vertex);
                          mesh.AddSegment(el);
                          vi += el.GetNP();
                        }

                      if(dim==2)
                        {
                          if(!have_2d_elements)
                          {
                            bc++;
                            have_2d_elements = true;
                            mesh.AddFaceDescriptor(FaceDescriptor(bc, 1, 0, 1));
                            mesh.SetBCName(bc-1, ngname);
                          }
                          auto el = ReadCGNSElement2D(type, vertices.Range(vi, vertices.Size()), first_vertex);
                          el.SetIndex(bc);
                          mesh.AddSurfaceElement(el);
                          vi += el.GetNP();
                        }

                      if(dim==3)
                        {
                          if(!have_3d_elements)
                          {
                            material++;
                            have_3d_elements = true;
                            mesh.SetMaterial(material, ngname);
                          }

                          auto el = ReadCGNSElement3D(type, vertices.Range(vi, vertices.Size()), first_vertex);
                          el.SetIndex(material);
                          mesh.AddVolumeElement(el);
                          vi += el.GetNP();
                          ne++;
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
                  int ne_section = nv/np;

                  if(dim==1)
                    {
                      for(auto i : Range(ne_section))
                        {
                          auto el = ReadCGNSElement1D(type, vertices.Range(np*i, np*(i+1)), first_vertex);
                          mesh.AddSegment(el);
                        }
                    }

                  if(dim==2)
                    {
                      bc++;
                      mesh.AddFaceDescriptor(FaceDescriptor(bc, 1, 0, 1));
                      for(auto i : Range(ne_section))
                        {
                          auto el = ReadCGNSElement2D(type, vertices.Range(np*i, np*(i+1)), first_vertex);
                          el.SetIndex(bc);
                          mesh.AddSurfaceElement(el);
                        }
                      mesh.SetBCName(bc-1, ngname);
                    }

                  if(dim==3)
                    {
                      material++;
                      for(auto i : Range(ne_section))
                        {
                          auto el = ReadCGNSElement3D(type, vertices.Range(np*i, np*(i+1)), first_vertex);
                          el.SetIndex(material);
                          mesh.AddVolumeElement(el);
                        }
                      mesh.SetMaterial(material, ngname);
                      ne += ne_section;
                    }
                }
            }
        }
    };
}

namespace netgen
{
  void ReadCGNSMesh (Mesh & mesh, const string & filename)
    {
      static Timer tall("CGNS::ReadMesh"); RegionTimer rtall(tall);
      int fn;
      cg_open(filename.c_str(),CG_MODE_READ,&fn);

      int base = 1;
      int nzones;
      cg_nzones(fn, base, &nzones);

      int bc = 0;
      int material = 0;

      for (auto zi : Range(1, nzones+1))
        {
          ZoneType_t zone_type;
          cg_zone_type(fn, base, zi, &zone_type);
          if(zone_type != Unstructured )
            {
              PrintMessage(2, "skipping zone with type ", cg_ZoneTypeName(zone_type) );
              continue;
            }
          cg::Zone zone(fn, base, zi);
          zone.ReadMesh( mesh );
        }
    }

  // Reads mesh and solutions of .csns file
  tuple<shared_ptr<Mesh>, vector<string>, vector<Array<double>>, vector<int>> ReadCGNSFile(string filename, int base)
    {
      static Timer tall("CGNS::ReadFile"); RegionTimer rtall(tall);
      int fn;
      cg_open(filename.c_str(),CG_MODE_READ,&fn);

      int nbases;
      cg_nbases(fn, &nbases);

      int nzones;
      cg_nzones(fn, base, &nzones);

      auto mesh = make_shared<Mesh>();

      int bc = 0;
      int material = 0;


      std::vector<string> names;
      std::vector<Array<double>> values;
      std::vector<int> locations;

      for (auto zi : Range(1, nzones+1))
        {
          ZoneType_t zone_type;
          cg_zone_type(fn, base, zi, &zone_type);
          if(zone_type != Unstructured )
            {
              clog << "skipping zone with type " << cg_ZoneTypeName(zone_type) << endl;
              continue;
            }
          cg::Zone zone(fn, base, zi);
          zone.ReadMesh( *mesh );
          zone.ReadSolutions( names, values, locations );
        }

      cg_close(fn);
      return std::make_tuple(mesh, names, values, locations);
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
}

#endif // NG_CGNS
