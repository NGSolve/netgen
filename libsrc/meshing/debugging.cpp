#include <meshing.hpp>
#include "debugging.hpp"

namespace netgen
{
    unique_ptr<Mesh> GetOpenElements( const Mesh & m, int dom, bool only_quads )
    {
        static Timer t("GetOpenElements"); RegionTimer rt(t);
        auto mesh = make_unique<Mesh>();
        *mesh = m;

        Array<bool, PointIndex> interesting_points(mesh->GetNP());
        interesting_points = false;

        mesh->FindOpenElements(dom);
        NgArray<Element2d> openelements;
        openelements = mesh->OpenElements();

        for (auto & el : openelements)
            for (auto i : el.PNums())
                interesting_points[i] = true;

        for (auto & el : mesh->VolumeElements())
        {
            int num_interesting_points = 0;

            for (auto pi : el.PNums())
                if(interesting_points[pi])
                    num_interesting_points++;

            if(num_interesting_points==0)
                el.Delete();
            el.SetIndex(num_interesting_points);
        }

        mesh->SetMaterial(1, "1_point");
        mesh->SetMaterial(2, "2_points");
        mesh->SetMaterial(3, "3_points");
        mesh->SetMaterial(4, "4_points");

        mesh->ClearSurfaceElements();

        for (auto & el : openelements)
            if(!only_quads || el.GetNP() == 4)
                mesh->AddSurfaceElement( el );

        mesh->Compress();
        return mesh;
    }

    unique_ptr<Mesh> FilterMesh( const Mesh & m, FlatArray<PointIndex> points, FlatArray<SurfaceElementIndex> sels, FlatArray<ElementIndex> els )
    {
        static Timer t("GetOpenElements"); RegionTimer rt(t);
        auto mesh_ptr = make_unique<Mesh>();
        auto & mesh = *mesh_ptr;
        mesh = m;

        Array<bool, PointIndex> keep_point(mesh.GetNP());
        Array<bool, SurfaceElementIndex> keep_sel(mesh.GetNSE());
        Array<bool, ElementIndex> keep_el(mesh.GetNE());
        mesh.LineSegments().DeleteAll();

        keep_point = false;
        for(auto pi : points)
            keep_point[pi] = true;

        auto set_keep = [&] (auto & input, auto & keep_array, auto & els)
        {
            keep_array = false;
            for(auto ind : input)
                keep_array[ind] = true;

            for(auto ind : Range(els))
            {
                bool & keep = keep_array[ind];
                if(keep) continue;

                for(auto pi : mesh[ind].PNums())
                    keep |= keep_point[pi];

                if(!keep)
                    mesh[ind].Delete();
            }

            for(auto i = 0; i<els.Size(); i++)
                if(els[i].IsDeleted())
                {
                    els.DeleteElement(i);
                    i--;
                }
        };

        set_keep(sels, keep_sel, mesh.SurfaceElements());
        set_keep(els, keep_el, mesh.VolumeElements());
        //mesh.Compress();

        return mesh_ptr;
    }

    void CheckMesh (const Mesh& mesh, MESHING_STEP step, const char * filename, int line)
    {
      if (step == MESHCONST_OPTVOLUME)
        {
          bool have_error = false;
          for (auto el : mesh.VolumeElements())
            {
              if(el.GetType() != TET)
                continue;
              double volume = el.Volume(mesh.Points());
              if (volume < 0)
                {
                  if(!have_error && line != -1)
                    cerr << "Negative volume in mesh at " << filename << ":" << line << endl;
                  have_error = true;
                  cerr << "volume of element " << el << " is negative: " << volume << endl;
                }
            }
          if (have_error)
            {
              string s;
              if(line != -1)
                {
                  s += filename;
                  s += ":";
                  s += ToString(line);
                  s += "\t";
                }
              throw Exception(s + "Negative volume");
            }

        CheckElementsAroundEdges(mesh);
        }
    }

    void CheckElementsAroundEdges (const Mesh& mesh)
    {
      static Mesh last_good_mesh;

      Array<std::tuple<PointIndex, PointIndex>> edges;
      auto elementsonnode = mesh.CreatePoint2ElementTable();
      BuildEdgeList(mesh, elementsonnode, edges);
      mesh.BoundaryEdge(1, 2); // trigger build of boundary edges

      ArrayMem<ElementIndex, 20> hasbothpoints;
      for (auto [pi0, pi1] : edges)
        {
          if (mesh.BoundaryEdge(pi0, pi1))
            continue;

          hasbothpoints.SetSize(0);
          for (ElementIndex ei : elementsonnode[pi0])
            if (mesh[ei].PNums().Contains(pi1))
              hasbothpoints.Append(ei);

          bool skip = false;
          for (ElementIndex ei : hasbothpoints)
            {
              if (mesh[ei].GetType() != TET)
                {
                  skip = true;
                  break;
                }
            }
          if (skip)
            continue;

          int nsuround = hasbothpoints.Size();
          ArrayMem<PointIndex, 50> suroundpts(nsuround + 1);
          suroundpts = PointIndex::INVALID;
          ArrayMem<bool, 50> tetused(nsuround);
          tetused = false;
          tetused[0] = true;

          auto el = mesh[hasbothpoints[0]];
          PointIndex pi2 = PointIndex::INVALID;
          PointIndex pi3 = PointIndex::INVALID;
          for (auto pi : el.PNums())
            if (pi != pi0 && pi != pi1)
              {
                pi3 = pi2;
                pi2 = pi;
              }
          suroundpts[0] = pi2;
          suroundpts[1] = pi3;

          for (auto i : Range(2, nsuround + 1))
            {
              PointIndex oldpi = suroundpts[i - 1];
              PointIndex newpi = PointIndex::INVALID;

              for (int k = 0; k < nsuround && !newpi.IsValid(); k++)
                if (!tetused[k])
                  {
                    const Element& nel = mesh[hasbothpoints[k]];
                    for (int k2 = 0; k2 < 4 && !newpi.IsValid(); k2++)
                      if (nel[k2] == oldpi)
                        {
                          newpi = nel[0] - pi0 + nel[1] - pi1 + nel[2] - oldpi + nel[3];
                          tetused[k] = true;
                          suroundpts[i] = newpi;

                          ArrayMem<PointIndex, 4> nelpts{nel[0], nel[1], nel[2], nel[3]};
                          ArrayMem<PointIndex, 4> check_points{pi0, pi1, oldpi, newpi};
                          QuickSort(check_points);
                          QuickSort(nelpts);
                          if (check_points != nelpts)
                            {
                              cerr << __FILE__ << ":" << __LINE__ << "\tFound error" << endl;
                              cerr << "i = " << i << endl;
                              cerr << "oldpi = " << oldpi << endl;
                              cerr << "newpi = " << newpi << endl;
                              cerr << "Elements: " << endl;
                              cerr << "nel " << nel << endl;
                              for (auto ei : hasbothpoints)
                                cerr << mesh[ei] << endl;
                              cerr << endl;
                              cerr << "check_points: " << check_points << endl;
                              cerr << "nelpts: " << nelpts << endl;
                              cerr << "hasbothpoints: " << hasbothpoints << endl;
                              cerr << "suroundpts: " << suroundpts << endl;
                              throw Exception("Found error");
                            }
                        }
                  }
            }
          if (suroundpts.Last() != suroundpts[0])
            {
              cerr << __FILE__ << ":" << __LINE__ << "\tFound error" << endl;
              cerr << "hasbothpoints: " << hasbothpoints << endl;
              cerr << "suroundpts: " << suroundpts << endl;
              for (auto ei : hasbothpoints)
                cerr << mesh[ei] << endl;
              throw Exception("Found error");
            }
        }
    }
} // namespace netgen
