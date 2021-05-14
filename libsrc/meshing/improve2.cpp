#include <mystdlib.h>

#include "meshing.hpp"
#include <opti.hpp>


namespace netgen
{

  class trionedge
  {
  public:
    SurfaceElementIndex tnr;
    int sidenr;

    trionedge () { tnr = 0; sidenr = 0; }
    trionedge (SurfaceElementIndex atnr, int asidenr)
    { tnr = atnr; sidenr = asidenr; }
  };


  bool MeshOptimize2d :: EdgeSwapping (const int usemetric,
    Array<Neighbour> &neighbors,
    Array<bool> &swapped,
    const SurfaceElementIndex t1, const int o1,
    const int t,
    Array<int,PointIndex> &pdef,
    const bool check_only)
  {
    bool should;
    bool do_swap = false;

    SurfaceElementIndex t2 = neighbors[t1].GetNr (o1);
    int o2 = neighbors[t1].GetOrientation (o1);

    if (t2 == -1) return false;
    if (swapped[t1] || swapped[t2]) return false;

    const int faceindex = mesh[t1].GetIndex();
    const int surfnr = mesh.GetFaceDescriptor (faceindex).SurfNr();

    PointIndex pi1 = mesh[t1].PNumMod(o1+1+1);
    PointIndex pi2 = mesh[t1].PNumMod(o1+1+2);
    PointIndex pi3 = mesh[t1].PNumMod(o1+1);
    PointIndex pi4 = mesh[t2].PNumMod(o2+1);

    PointGeomInfo gi1 = mesh[t1].GeomInfoPiMod(o1+1+1);
    PointGeomInfo gi2 = mesh[t1].GeomInfoPiMod(o1+1+2);
    PointGeomInfo gi3 = mesh[t1].GeomInfoPiMod(o1+1);
    PointGeomInfo gi4 = mesh[t2].GeomInfoPiMod(o2+1);

    bool allowswap = true;

    Vec<3> auxvec1 = mesh[pi3]-mesh[pi4];
    Vec<3> auxvec2 = mesh[pi1]-mesh[pi4];

    allowswap = allowswap && fabs(1.-(auxvec1*auxvec2)/(auxvec1.Length()*auxvec2.Length())) > 1e-4;

    if(!allowswap)
        return false;

    // normal of new
    Vec<3> nv1 = Cross (auxvec1, auxvec2);

    auxvec1 = mesh.Point(pi4)-mesh.Point(pi3);
    auxvec2 = mesh.Point(pi2)-mesh.Point(pi3);
    allowswap = allowswap && fabs(1.-(auxvec1*auxvec2)/(auxvec1.Length()*auxvec2.Length())) > 1e-4;


    if(!allowswap)
        return false;

    Vec<3> nv2 = Cross (auxvec1, auxvec2);


    // normals of original
    Vec<3> nv3 = Cross (mesh[pi1]-mesh[pi4], mesh[pi2]-mesh[pi4]);
    Vec<3> nv4 = Cross (mesh[pi2]-mesh[pi3], mesh[pi1]-mesh[pi3]);

    nv3 *= -1;
    nv4 *= -1;
    nv3.Normalize();
    nv4.Normalize();

    nv1.Normalize();
    nv2.Normalize();

    auto nvp3 = geo.GetNormal (surfnr, mesh.Point(pi3), &gi3);

    nvp3.Normalize();

    auto nvp4 = geo.GetNormal (surfnr, mesh.Point(pi4), &gi4);

    nvp4.Normalize();



    double critval = cos (M_PI / 6);  // 30 degree
    allowswap = allowswap &&
        (nv1 * nvp3 > critval) &&
        (nv1 * nvp4 > critval) &&
        (nv2 * nvp3 > critval) &&
        (nv2 * nvp4 > critval) &&
        (nvp3 * nv3 > critval) &&
        (nvp4 * nv4 > critval);


    double horder = Dist (mesh[pi1], mesh[pi2]);

    if ( // nv1 * nv2 >= 0 &&
            nv1.Length() > 1e-3 * horder * horder &&
            nv2.Length() > 1e-3 * horder * horder &&
            allowswap )
      {
        if (!usemetric)
          {
            int e = pdef[pi1] + pdef[pi2] - pdef[pi3] - pdef[pi4];
            double d =
                Dist2 (mesh[pi1], mesh[pi2]) -
                Dist2 (mesh[pi3], mesh[pi4]);

            should = e >= t && (e > 2 || d > 0);
          }
        else
          {
            double loch = mesh.GetH(mesh[pi1]);
            should =
                CalcTriangleBadness (mesh[pi4], mesh[pi3], mesh[pi1], metricweight, loch) +
                CalcTriangleBadness (mesh[pi3], mesh[pi4], mesh[pi2], metricweight, loch) <
                CalcTriangleBadness (mesh[pi1], mesh[pi2], mesh[pi3], metricweight, loch) +
                CalcTriangleBadness (mesh[pi2], mesh[pi1], mesh[pi4], metricweight, loch);
          }

        if (allowswap)
          {
            Element2d sw1 (pi4, pi3, pi1);
            Element2d sw2 (pi3, pi4, pi2);

            int legal1 =
                mesh.LegalTrig (mesh[t1]) +
                mesh.LegalTrig (mesh[t2]);
            int legal2 =
                mesh.LegalTrig (sw1) + mesh.LegalTrig (sw2);

            if (legal1 < legal2) should = true;
            if (legal2 < legal1) should = false;
          }

        do_swap = should;
        if (should && !check_only)
          {
            // do swapping !

            mesh[t1] = { { pi1, gi1 }, { pi4, gi4 }, { pi3, gi3 } };
            mesh[t2] = { { pi2, gi2 }, { pi3, gi3 }, { pi4, gi4 } };

            pdef[pi1]--;
            pdef[pi2]--;
            pdef[pi3]++;
            pdef[pi4]++;

            swapped[t1] = true;
            swapped[t2] = true;
          }
      }
    return do_swap;
  }

 
  void MeshOptimize2d :: EdgeSwapping (int usemetric)
  {
    static Timer timer("EdgeSwapping (2D)"); RegionTimer reg(timer);
    static Timer timer_nb("EdgeSwapping-Find neighbors");
    if (usemetric)
      PrintMessage (3, "Edgeswapping, metric");
    else
      PrintMessage (3, "Edgeswapping, topological");

    static Timer timerstart("EdgeSwapping 2D start");
    timerstart.Start();

    Array<SurfaceElementIndex> seia;
    bool mixed = false;

    if(faceindex==0)
      {
        seia.SetSize(mesh.GetNSE());
        ParallelFor( Range(seia), [&] (auto i) NETGEN_LAMBDA_INLINE
            {
              SurfaceElementIndex sei(i);
              seia[i] = sei;
              if (mesh[sei].GetNP() != 3)
              {
                  const auto & sel = mesh[sei];
                  for(auto i : Range(sel.GetNP()))
                    if(mesh[sel[i]].Type() == INNERPOINT)
                      mixed = true;
              }
            });
      }
    else
      {
        mesh.GetSurfaceElementsOfFace (faceindex, seia);
        for (SurfaceElementIndex sei : seia)
            if (mesh[sei].GetNP() != 3)
                mixed = true;
      }

    if(mixed)
        return GenericImprove();
      
    Array<Neighbour> neighbors(mesh.GetNSE());
    auto elements_on_node = mesh.CreatePoint2SurfaceElementTable(faceindex);

    Array<bool> swapped(mesh.GetNSE());
    Array<int,PointIndex> pdef(mesh.GetNP());
    Array<double,PointIndex> pangle(mesh.GetNP());

    static const double minangle[] = { 0, 1.481, 2.565, 3.627, 4.683, 5.736, 7, 9 };


    if(faceindex == 0)
      {
        ParallelFor( Range(pangle), [&] (auto i) NETGEN_LAMBDA_INLINE
            {
              pangle[i] = 0.0;
            });
      }
    else
      {
        ParallelFor( Range(seia), [&] (auto i) NETGEN_LAMBDA_INLINE
            {
              const Element2d & sel = mesh[seia[i]];
              for (int j = 0; j < 3; j++)
                  pangle[sel[j]] = 0.0;
            });
      }

    ParallelFor( Range(seia), [&] (auto i) NETGEN_LAMBDA_INLINE
        {
          const Element2d & sel = mesh[seia[i]];
          for (int j = 0; j < 3; j++)
            {
              POINTTYPE typ = mesh[sel[j]].Type();
              if (typ == FIXEDPOINT || typ == EDGEPOINT)
                {
                  AtomicAdd(pangle[sel[j]],
                    Angle (mesh[sel[(j+1)%3]] - mesh[sel[j]],
                           mesh[sel[(j+2)%3]] - mesh[sel[j]]));
                }
            }
       });

    ParallelFor( Range(seia), [&] (auto i) NETGEN_LAMBDA_INLINE
        {
          const Element2d & sel = mesh[seia[i]];
          for (int j = 0; j < 3; j++)
            {
              PointIndex pi = sel[j];
              if (mesh[pi].Type() == INNERPOINT || mesh[pi].Type() == SURFACEPOINT)
                pdef[pi] = -6;
              else
                for (int j = 0; j < 8; j++)
                  if (pangle[pi] >= minangle[j])
                    pdef[pi] = -1-j;
            }
       });

    ParallelFor( Range(seia), [this, &pdef, &neighbors, &seia, &elements_on_node] (auto i) NETGEN_LAMBDA_INLINE
        {
          auto sei = seia[i];
          for (PointIndex pi : mesh[sei].template PNums<3>())
              AsAtomic(pdef[pi])++;
          for (int j = 0; j < 3; j++)
            {
              neighbors[sei].SetNr (j, -1);
              neighbors[sei].SetOrientation (j, 0);
            }

          const auto sel = mesh[sei];
          for (int j = 0; j < 3; j++)
            {
              PointIndex pi1 = sel.PNumMod(j+2);
              PointIndex pi2 = sel.PNumMod(j+3);

              for (auto sei_other : elements_on_node[pi1])
                {
                  if(sei_other==sei) continue;
                  const auto & other = mesh[sei_other];
                  int pi1_other = -1;
                  int pi2_other = -1;
                  bool common_edge = false;
                  for (int k = 0; k < 3; k++)
                    {
                      if(other[k] == pi1)
                          pi1_other = k;
                      if(other[k] == pi2)
                        {
                          pi2_other = k;
                          common_edge = true;
                        }
                    }

                  if(common_edge)
                    {
                      neighbors[sei].SetNr (j, sei_other);
                      neighbors[sei].SetOrientation (j, 3-pi1_other-pi2_other);
                    }
                }
            }
        });


    for (SurfaceElementIndex sei : seia)
      swapped[sei] = false;
    
    timerstart.Stop();
  


    Array<std::pair<SurfaceElementIndex,int>> improvement_candidates(3*seia.Size());
    atomic<int> cnt(0);

    int t = 4;
    bool done = false;
    while (!done && t >= 2)
      {
        cnt = 0;
        ParallelFor( Range(seia), [&] (auto i) NETGEN_LAMBDA_INLINE
          {
            SurfaceElementIndex t1 = seia[i];

            if (mesh[t1].IsDeleted())
              return;

            if (mesh[t1].GetIndex() != faceindex)
              return;

            if (multithread.terminate)
              throw NgException ("Meshing stopped");

            for (int o1 = 0; o1 < 3; o1++)
                if(EdgeSwapping(usemetric, neighbors, swapped, t1, o1, t, pdef, true))
                    improvement_candidates[cnt++]= std::make_pair(t1,o1);
          });

        auto elements_with_improvement = improvement_candidates.Range(cnt.load());
        QuickSort(elements_with_improvement);

        for (auto [t1,o1] : elements_with_improvement)
            done |= EdgeSwapping(usemetric, neighbors, swapped, t1, o1, t, pdef, false);
	t--;
      }

    mesh.SetNextTimeStamp();
  }




  double CombineImproveEdge( Mesh & mesh,
                           const Table<SurfaceElementIndex, PointIndex> & elementsonnode,
                           Array<Vec<3>, PointIndex> & normals,
                           Array<bool, PointIndex> & fixed,
                           PointIndex pi1, PointIndex pi2,
                           bool check_only = true)
  {
    Vec<3> nv;
    ArrayMem<SurfaceElementIndex, 20> hasonepi, hasbothpi;

    if (!pi1.IsValid() || !pi2.IsValid())
        return 0.0;

    bool debugflag = 0;

    if (debugflag)
      {
        (*testout) << "Combineimprove "
            << "pi1 = " << pi1 << " pi2 = " << pi2 << endl;
      }

    /*
    // save version:
    if (fixed.Get(pi1) || fixed.Get(pi2))
    return 0.0;
    if (pi2 < pi1) swap (pi1, pi2);
    */

    // more general
    if (fixed[pi2])
        Swap (pi1, pi2);

    if (fixed[pi2])
        return 0.0;

    double loch = mesh.GetH (mesh[pi1]);

    for (SurfaceElementIndex sei2 : elementsonnode[pi1])
      {
        const Element2d & el2 = mesh[sei2];

        if (el2.IsDeleted()) continue;

        if (el2[0] == pi2 || el2[1] == pi2 || el2[2] == pi2)
          {
            hasbothpi.Append (sei2);
            nv = Cross (Vec3d (mesh[el2[0]], mesh[el2[1]]),
                    Vec3d (mesh[el2[0]], mesh[el2[2]]));
          }
        else
          {
            hasonepi.Append (sei2);
          }
      }

    if(hasbothpi.Size()==0)
        return 0.0;


    nv = normals[pi1];


    for (SurfaceElementIndex sei2 :  elementsonnode[pi2])
      {
        const Element2d & el2 = mesh[sei2];
        if (el2.IsDeleted()) continue;
        if (!el2.PNums<3>().Contains (pi1))
            hasonepi.Append (sei2);
      }

    double bad1 = 0;
    int illegal1 = 0, illegal2 = 0;
    /*
       for (SurfaceElementIndex sei : hasonepi)
       {
       const Element2d & el = mesh[sei];
       bad1 += CalcTriangleBadness (mesh[el[0]], mesh[el[1]], mesh[el[2]],
       nv, -1, loch);
       illegal1 += 1-mesh.LegalTrig(el);
       }
       */
    for (const Element2d & el : mesh.SurfaceElements()[hasonepi])
      {
        bad1 += CalcTriangleBadness (mesh[el[0]], mesh[el[1]], mesh[el[2]],
                nv, -1, loch);
        illegal1 += 1-mesh.LegalTrig(el);
      }

    for (int k = 0; k < hasbothpi.Size(); k++)
      {
        const Element2d & el = mesh[hasbothpi[k]];
        bad1 += CalcTriangleBadness (mesh[el[0]], mesh[el[1]], mesh[el[2]],
                nv, -1, loch);
        illegal1 += 1-mesh.LegalTrig(el);
      }

    double bad2 = 0;
    for (int k = 0; k < hasonepi.Size(); k++)
      {
        Element2d el = mesh[hasonepi[k]];
        for (auto i : Range(3))
            if(el[i]==pi2)
                el[i] = pi1;

        double err =
            CalcTriangleBadness (mesh[el[0]], mesh[el[1]], mesh[el[2]],
                    nv, -1, loch);
        bad2 += err;

        Vec<3> hnv = Cross (Vec3d (mesh[el[0]],
                    mesh[el[1]]),
                Vec3d (mesh[el[0]],
                    mesh[el[2]]));
        if (hnv * nv < 0)
            bad2 += 1e10;

        for (int l = 0; l < 3; l++)
          {
            if ( (normals[el[l]] * nv) < 0.5)
                bad2 += 1e10;
          }

        illegal2 += 1-mesh.LegalTrig(el);
      }

    if (debugflag)
      {
        (*testout) << "bad1 = " << bad1 << ", bad2 = " << bad2 << endl;
      }

    bool should = (illegal2<=illegal1 && bad2 < bad1 && bad2 < 1e4);
    if(illegal2 < illegal1)
      {
        should = true;
        bad1 += 1e4;
      }

    double d_badness = should * (bad2-bad1);

    if(check_only)
        return d_badness;

    if (should)
      {
        /*
           (*testout) << "combine !" << endl;
           (*testout) << "bad1 = " << bad1 << ", bad2 = " << bad2 << endl;
           (*testout) << "illegal1 = " << illegal1 << ", illegal2 = " << illegal2 << endl;
           (*testout) << "loch = " << loch << endl;
           */

        PointGeomInfo gi;
        // bool gi_set(false);

        /*
           Element2d *el1p(NULL);
           int l = 0;
           while(mesh[elementsonnode[pi1][l]].IsDeleted() && l<elementsonnode[pi1].Size()) l++;
           if(l<elementsonnode[pi1].Size())
           el1p = &mesh[elementsonnode[pi1][l]];
           else
           cerr << "OOPS!" << endl;

           for (l = 0; l < el1p->GetNP(); l++)
           if ((*el1p)[l] == pi1)
           {
           gi = el1p->GeomInfoPi (l+1);
        // gi_set = true;
        }
        */
        for (SurfaceElementIndex sei : elementsonnode[pi1])
          {
            const Element2d & el1p = mesh[sei];
            if (el1p.IsDeleted()) continue;

            for (int l = 0; l < el1p.GetNP(); l++)
                if (el1p[l] == pi1)
                    // gi = el1p.GeomInfoPi (l+1);
                    gi = el1p.GeomInfo()[l];
            break;
          }


        // (*testout) << "Connect point " << pi2 << " to " << pi1 << "\n";
        // for (int k = 0; k < elementsonnode[pi2].Size(); k++)
        for (SurfaceElementIndex sei2 : elementsonnode[pi2])
          {
            Element2d & el = mesh[sei2];
            if (el.IsDeleted()) continue;
            if (el.PNums().Contains(pi1)) continue;

            for (auto l : Range(el.GetNP()))
              {
                if (el[l] == pi2)
                  {
                    el[l] = pi1;
                    el.GeomInfo()[l] = gi;
                  }

                fixed[el[l]] = true;
              }
          }

        for (auto sei : hasbothpi)
            mesh[sei].Delete();

      }
    return d_badness;
  }

  void MeshOptimize2d :: CombineImprove ()
  {
    SplitImprove();
    PrintMessage (3, "Combine improve");

    if (multithread.terminate)
        throw NgException ("Meshing stopped");

    static Timer timer ("Combineimprove 2D");
    RegionTimer reg (timer);

    static Timer timerstart ("Combineimprove 2D start");
    timerstart.Start();


    static Timer timerstart1 ("Combineimprove 2D start1");
    timerstart1.Start();

    
    Array<SurfaceElementIndex> seia;

    if(faceindex)
        mesh.GetSurfaceElementsOfFace (faceindex, seia);
    else
      {
        seia.SetSize(mesh.GetNSE());
        ParallelFor( IntRange(mesh.GetNSE()), [&seia] (auto i) NETGEN_LAMBDA_INLINE
                { seia[i] = i; });
      }

    bool mixed = false;
    ParallelFor( Range(seia), [&] (auto i) NETGEN_LAMBDA_INLINE
            {
                if (mesh[seia[i]].GetNP() != 3)
                    mixed = true;
            });

    if(mixed)
        return;

    int np = mesh.GetNP();

    auto elementsonnode = mesh.CreatePoint2SurfaceElementTable(faceindex);

    int ntasks = ngcore::TaskManager::GetMaxThreads();
    Array<std::tuple<PointIndex, PointIndex>> edges;

    BuildEdgeList( mesh, elementsonnode, edges );

    Array<bool,PointIndex> fixed(np);
    ParallelFor( fixed.Range(), [&fixed] (auto i) NETGEN_LAMBDA_INLINE
            { fixed[i] = false; });

    ParallelFor( edges.Range(), [&] (auto i) NETGEN_LAMBDA_INLINE
            {
              auto [pi0, pi1] = edges[i];
              if (mesh.IsSegment (pi0, pi1))
                {
                  fixed[pi0] = true;
                  fixed[pi1] = true;
                }
            });

    timerstart1.Stop();

    ParallelFor( mesh.LockedPoints().Range(), [&] (auto i) NETGEN_LAMBDA_INLINE
            {
              fixed[mesh.LockedPoints()[i]] = true;
            });


    Array<Vec<3>,PointIndex> normals(np);

    ParallelFor( mesh.Points().Range(), [&] (auto pi) NETGEN_LAMBDA_INLINE
        {
            if (elementsonnode[pi].Size())
              {
                Element2d & hel = mesh[elementsonnode[pi][0]];
                for (int k = 0; k < 3; k++)
                  if (hel[k] == pi)
                    {
                      const int faceindex = hel.GetIndex();
                      const int surfnr = mesh.GetFaceDescriptor (faceindex).SurfNr();
                      normals[pi] = geo.GetNormal (surfnr, mesh[pi], &hel.GeomInfoPi(k+1));
                      break;
                    }
              }
        }, TasksPerThread(4));

    timerstart.Stop();

    // Find edges with improvement
    Array<std::tuple<double, int>> candidate_edges(edges.Size());
    std::atomic<int> improvement_counter(0);

    ParallelFor( Range(edges), [&] (auto i) NETGEN_LAMBDA_INLINE
      {
        auto [pi1, pi2] = edges[i];
        double d_badness = CombineImproveEdge(mesh, elementsonnode, normals, fixed, pi1, pi2, true);
        if(d_badness < 0.0)
            candidate_edges[improvement_counter++] = make_tuple(d_badness, i);
      }, TasksPerThread(4));

    auto edges_with_improvement = candidate_edges.Part(0, improvement_counter.load());
    QuickSort(edges_with_improvement);

    for(auto [d_badness, ei] : edges_with_improvement)
      {
        auto [pi1, pi2] = edges[ei];
        CombineImproveEdge(mesh, elementsonnode, normals, fixed, pi1, pi2, false);
      }

    //  mesh.Compress();
    mesh.SetNextTimeStamp();
  }

  void MeshOptimize2d :: SplitImprove()
  {
    if (!faceindex)
      {
        PrintMessage (3, "Split improve");

        mesh.CalcSurfacesOfNode(); // TODO: needed?
        for (faceindex = 1; faceindex <= mesh.GetNFD(); faceindex++)
          {
            SplitImprove();

            if (multithread.terminate)
                throw NgException ("Meshing stopped");
          }

        faceindex = 0;
        mesh.Compress(); // TODO: needed?
        return;
      }

    Array<SurfaceElementIndex> elements;
    mesh.GetSurfaceElementsOfFace (faceindex, elements);

    // return if we have quads in this surface
    for (auto & ei : elements)
        if (mesh[ei].GetNP() != 3)
            return;

    // maps from edges to adjacent trigs
    INDEX_2_HASHTABLE<tuple<SurfaceElementIndex, SurfaceElementIndex>> els_on_edge(2*elements.Size() + 2);

    // build els_on_edge table
    for (SurfaceElementIndex sei : elements)
      {
        const Element2d & sel = mesh[sei];

        for (int j = 0; j < 3; j++)
          {
            PointIndex pi1 = sel.PNumMod(j+2);
            PointIndex pi2 = sel.PNumMod(j+3);

            if (mesh.IsSegment (pi1, pi2)) continue;

            INDEX_2 ii2 (pi1, pi2);
            ii2.Sort();
            if (els_on_edge.Used (ii2))
              {
                auto els = els_on_edge.Get(ii2);
                get<1>(els) = sei;
                els_on_edge.Set(ii2, els);
              }
            else
              {
                els_on_edge.Set (ii2, make_tuple(sei, sei));
              }
          }
      }

    // split edges of illegal trigs
    for (SurfaceElementIndex sei : elements)
      {
        Element2d & sel = mesh[sei];

        if (sel.IsDeleted()) continue;

        // TODO: split also bad trigs, nut just illegal ones
        if (mesh.LegalTrig(sel)) continue;

        // find longest edge
        INDEX_2 edge;
        double edge_len = 0;
        PointIndex pi1, pi2, pi3, pi4;
        PointGeomInfo gi1, gi2, gi3, gi4;
        for(auto j : Range(1,4))
          {
            auto test_pi1 = sel.PNumMod(j);
            auto test_pi2 = sel.PNumMod(j+1);
            if (mesh.IsSegment(test_pi1, test_pi2))
              continue;
            auto len = (mesh[test_pi2]-mesh[test_pi1]).Length();
            if(len > edge_len)
              {
                edge = {test_pi1, test_pi2};
                edge.Sort();
                edge_len = len;
                pi1 = test_pi1;
                pi2 = test_pi2;
                pi3 = sel.PNumMod(j+2);
                gi1 = sel.GeomInfoPiMod(j);
                gi2 = sel.GeomInfoPiMod(j+1);
                gi3 = sel.GeomInfoPiMod(j+2);
              }
          }
        if(!edge_len)
          throw Exception("Couldn't find edge to split, something is wrong");
        // get neighbor element
        auto els = els_on_edge.Get(edge);
        SurfaceElementIndex other_i = get<0>(els);
        if(other_i==sei) other_i = get<1>(els);
        auto & other = mesh[other_i];

        // find opposite point of neighbor element
        for (int j = 0; j < 3; j++)
          if(other[j]!=pi1 && other[j]!=pi2)
            {
              pi4 = other[j];
              gi4 = other.GeomInfoPi(j);
              break;
            }

        // split edge pi1,pi2
        Point<3> p5;
        PointIndex pi5;
        PointGeomInfo gi5;

        geo.PointBetween(mesh[pi1], mesh[pi2], 0.5,
                         faceindex,
                         gi1, gi2, p5, gi5);

        pi5 = mesh.AddPoint(p5);

        Element2d e1(3);
        e1.SetIndex(faceindex);
        e1={ {pi1,gi1}, {pi5,gi5}, {pi3,gi3} };
        mesh.AddSurfaceElement( e1 );

        Element2d e2(3);
        e2.SetIndex(faceindex);
        e2 ={ {pi5,gi5}, {pi2,gi2}, {pi3,gi3} };
        mesh.AddSurfaceElement( e2 );

        Element2d e3(3);
        e3.SetIndex(faceindex);
        e3 ={ {pi1,gi1}, {pi4,gi4}, {pi5,gi5} };
        mesh.AddSurfaceElement( e3 );

        Element2d e4(3);
        e4.SetIndex(faceindex);
        e4 ={ {pi4,gi4}, {pi2,gi2}, {pi5,gi5} };
        mesh.AddSurfaceElement( e4 );

        sel.Delete();
        other.Delete();
      }

    mesh.SetNextTimeStamp();
  }

  void MeshOptimize2d :: CheckMeshApproximation (Mesh & mesh)
  {
    // Check angles between elements and normals at corners
    /*
  
    int i, j;
    int ne = mesh.GetNSE();
    int surfnr;
  
    Vec3d n, ng;
    NgArray<Vec3d> ngs(3);

    (*mycout) << "Check Surface Approximation" << endl;
    (*testout) << "Check Surface Approximation" << endl;

    for (i = 1; i <= ne; i++)
    {
    const Element2d & el = mesh.SurfaceElement(i);
    surfnr = mesh.GetFaceDescriptor (el.GetIndex()).SurfNr();
    Vec3d n = Cross (mesh.Point (el.PNum(1)) - mesh.Point (el.PNum(2)),
    mesh.Point (el.PNum(1)) - mesh.Point (el.PNum(3)));
    n /= n.Length();

    for (j = 1; j <= el.GetNP(); j++)
    {
    SelectSurfaceOfPoint (mesh.Point(el.PNum(j)), el.GeomInfoPi(j));
    GetNormalVector (surfnr, mesh.Point(el.PNum(j)), ng);
    ng /= ng.Length();
    ngs.Elem(j) = ng;

    double angle =  (180.0 / M_PI) * Angle (n, ng);
    if (angle > 60)
    {
    (*testout) << "el " << i << " node " << el.PNum(j)
    << "has angle = " << angle << endl;
    }
    }	

    for (j = 1; j <= 3; j++)
    {
    double angle =  (180.0 / M_PI) * Angle (ngs.Get(j), ngs.Get(j%3+1));
    if (angle > 60)
    {
    (*testout) << "el " << i << " node-node " 
    << ngs.Get(j) << " - " << ngs.Get(j%3+1)
    << " has angle = " << angle << endl;
    }
    }
    }
    */
  }
}
