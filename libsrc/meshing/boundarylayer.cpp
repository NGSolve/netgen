#include <mystdlib.h>
#include "meshing.hpp"

namespace netgen
{

   void InsertVirtualBoundaryLayer (Mesh & mesh)
   {
      cout << "Insert virt. b.l." << endl;

      int surfid;

      cout << "Boundary Nr:";
      cin >> surfid;

      int i;
      int np = mesh.GetNP();

      cout << "Old NP: " << mesh.GetNP() << endl;
      cout << "Trigs: " << mesh.GetNSE() << endl;

      NgBitArray bndnodes(np);
      NgArray<int> mapto(np);

      bndnodes.Clear();
      for (i = 1; i <= mesh.GetNSeg(); i++)
      {
         int snr = mesh.LineSegment(i).edgenr;
         cout << "snr = " << snr << endl;
         if (snr == surfid)
         {
            bndnodes.Set (mesh.LineSegment(i)[0]);
            bndnodes.Set (mesh.LineSegment(i)[1]);
         }
      }
      for (i = 1; i <= mesh.GetNSeg(); i++)
      {
         int snr = mesh.LineSegment(i).edgenr;
         if (snr != surfid)
         {
            bndnodes.Clear (mesh.LineSegment(i)[0]);
            bndnodes.Clear (mesh.LineSegment(i)[1]);
         }
      }

      for (i = 1; i <= np; i++)
        {
          if (bndnodes.Test(i))
            mapto.Elem(i) = mesh.AddPoint (mesh.Point (i));
          else
            mapto.Elem(i) = 0;
        }

      for (i = 1; i <= mesh.GetNSE(); i++)
      {
         Element2d & el = mesh.SurfaceElement(i);
         for (int j = 1; j <= el.GetNP(); j++)
            if (mapto.Get(el.PNum(j)))
               el.PNum(j) = mapto.Get(el.PNum(j));
      }


      int nq = 0;
      for (i = 1; i <= mesh.GetNSeg(); i++)
      {
         int snr = mesh.LineSegment(i).edgenr;
         if (snr == surfid)
         {
            int p1 = mesh.LineSegment(i)[0];
            int p2 = mesh.LineSegment(i)[1];
            int p3 = mapto.Get (p1);
            if (!p3) p3 = p1;
            int p4 = mapto.Get (p2);
            if (!p4) p4 = p2;

            Element2d el(QUAD);
            el.PNum(1) = p1;
            el.PNum(2) = p2;
            el.PNum(3) = p3;
            el.PNum(4) = p4;
            el.SetIndex (2);
            mesh.AddSurfaceElement (el);
            nq++;
         }
      }

      cout << "New NP: " << mesh.GetNP() << endl;
      cout << "Quads: " << nq << endl;
   }





/*
   Philippose Rajan - 11 June 2009

   Function to calculate the surface normal at a given 
   vertex of a surface element, with respect to that 
   surface element.

   This function is used by the boundary layer generation 
   function, in order to calculate the effective direction 
   in which the prismatic layer should grow
*/
  inline Vec<3> GetSurfaceNormal(Mesh & mesh, const Element2d & el)
   {
     auto v0 = mesh[el[0]];
     auto v1 = mesh[el[1]];
     auto v2 = mesh[el[2]];
     Vec<3> vec1 = v1-v0;
     Vec<3> vec2 = v2-v0;
     Vec<3> normal = Cross(vec1, vec2);
     normal.Normalize();
     return normal;
   }

/*
    Philippose Rajan - 11 June 2009
    modified by Christopher Lackner Apr 2020
    
    Added an initial experimental function for 
    generating prismatic boundary layers on 
    a given set of surfaces.
    
    The number of layers, height of the first layer 
    and the growth / shrink factor can be specified 
    by the user

    Currently, the layer height is calculated using:
    height = h_first_layer * (growth_factor^(num_layers - 1))
*/
   void GenerateBoundaryLayer (Mesh & mesh, const BoundaryLayerParameters & blp)
   {
      PrintMessage(1, "Generating boundary layer...");
      PrintMessage(3, "Old NP: ", mesh.GetNP());
      PrintMessage(3, "Old NSE: ",mesh.GetNSE());

      map<tuple<int, int, int>, int> domains_to_surf_index;
      map<tuple<PointIndex, PointIndex>, int> pi_to_edgenr;

      map<int, int> last_layer_surface_index_map;
      int max_surface_index = mesh.GetNFD();

      int max_edge_nr = -1;
      for(const auto& seg : mesh.LineSegments())
        if(seg.edgenr > max_edge_nr)
          max_edge_nr = seg.edgenr;

      for(int layer = blp.heights.Size(); layer >= 1; layer--)
        {
          PrintMessage(3, "Generating layer: ", layer);

          auto map_surface_index = [&](auto si)
            {
              if(last_layer_surface_index_map.find(si) == last_layer_surface_index_map.end())
                {
                  last_layer_surface_index_map[si] = ++max_surface_index;
                  auto& old_fd = mesh.GetFaceDescriptor(si);
                  int domout = blp.outside ? old_fd.DomainOut() : blp.new_matnrs[layer-1];
                  int domin = blp.outside ? blp.new_matnrs[layer-1] : old_fd.DomainIn();
                  FaceDescriptor fd(max_surface_index-1,
                                    domin, domout, -1);
                  fd.SetBCProperty(max_surface_index);
                  mesh.AddFaceDescriptor(fd);
                  return max_surface_index;
                }
              return last_layer_surface_index_map[si];
            };

          mesh.UpdateTopology();
          auto& meshtopo = mesh.GetTopology();

          auto layerht = blp.heights[layer-1];

          PrintMessage(5, "Layer Height = ", layerht);

          // Need to store the old number of points and
          // surface elements because there are new points and
          // surface elements being added during the process
          int np = mesh.GetNP();
          int nse = mesh.GetNSE();
          int ne = mesh.GetNE();

          // Safety measure to ensure no issues with mesh
          // consistency
          int nseg = mesh.GetNSeg();

          // Indicate which points need to be remapped
          BitArray bndnodes(np+1);  // big enough for 1-based array

          // Map of the old points to the new points
          Array<PointIndex, PointIndex> mapto(np);

          // Growth vectors for the prismatic layer based on
          // the effective surface normal at a given point
          Array<Vec<3>, PointIndex> growthvectors(np);
          growthvectors = 0.;

          // Bit array to identify all the points belonging
          // to the surface of interest
          bndnodes.Clear();

          // Run through all the surface elements and mark the points
          // belonging to those where a boundary layer has to be created.
          // In addition, also calculate the effective surface normal
          // vectors at each of those points to determine the mesh motion
          // direction
          PrintMessage(3, "Marking points for remapping...");

          for(const auto& sel : mesh.SurfaceElements())
            if (blp.surfid.Contains(sel.GetIndex()))
              {
                auto n2 = GetSurfaceNormal(mesh,sel);
                if(!blp.outside)
                  n2 *= -1;
                for(auto pi : sel.PNums())
                  {
                    // Set the bitarray to indicate that the
                    // point is part of the required set
                    bndnodes.SetBit(pi);

                    // Add the surface normal to the already existent one
                    // (This gives the effective normal direction at corners
                    //  and curved areas)

                    auto& n1 = growthvectors[pi];
                    if(n1.Length() == 0) { n1 = n2; continue; }
                    auto n1n2 = n1 * n2;
                    auto n1n1 = n1 * n1;
                    auto n2n2 = n2 * n2;
                    if(n2n2 - n1n2*n1n2/n1n1 == 0) { n1 = n2; continue; }
                    n1 += (n2n2 - n1n2)/(n2n2 - n1n2*n1n2/n1n1) * (n2 - n1n2/n1n1 * n1);
                  }
              }

          if (!blp.grow_edges)
            for(const auto& sel : mesh.LineSegments())
              {
                bndnodes.Clear(sel[0]);
                bndnodes.Clear(sel[1]);
              }

          // Add additional points into the mesh structure in order to
          // clone the surface elements.
          // Also invert the growth vectors so that they point inwards,
          // and normalize them
          PrintMessage(3, "Cloning points and calculating growth vectors...");

          for (PointIndex pi = 1; pi <= np; pi++)
            {
              if (bndnodes.Test(pi))
                mapto[pi] = mesh.AddPoint(mesh[pi]);
              else
                mapto[pi].Invalidate();
            }

          // Add quad surface elements at edges for surfaces which
          // don't have boundary layers

          // Bit array to keep track of segments already processed
          BitArray segsel(nseg);

          // Set them all to "1" to initially activate all segments
          segsel.Set();

          PrintMessage(3, "Adding 2D Quad elements on required surfaces...");

          if(blp.grow_edges)
            for(SegmentIndex sei = 0; sei < nseg; sei++)
              {
                PointIndex seg_p1 = mesh[sei][0];
                PointIndex seg_p2 = mesh[sei][1];

                // Only go in if the segment is still active, and if both its
                // surface index is part of the "hit-list"
                if(segsel.Test(sei))
                  {
                    if(blp.surfid.Contains(mesh[sei].si))
                  {
                    // clear the bit to indicate that this segment has been processed
                    segsel.Clear(sei);

                    // Find matching segment pair on other surface
                    for(SegmentIndex sej = 0; sej < nseg; sej++)
                      {
                        PointIndex segpair_p1 = mesh[sej][1];
                        PointIndex segpair_p2 = mesh[sej][0];

                        // Find the segment pair on the neighbouring surface element
                        // Identified by: seg1[0] = seg_pair[1] and seg1[1] = seg_pair[0]
                        if(segsel.Test(sej) && ((segpair_p1 == seg_p1) && (segpair_p2 == seg_p2)))
                          {
                            // clear bit to indicate that processing of this segment is done
                            segsel.Clear(sej);

                            // Only worry about those surfaces which are not in the
                            // boundary layer list
                            if(!blp.surfid.Contains(mesh[sej].si))
                              {
                                SurfaceElementIndex pnt_commelem;
                                SetInvalid(pnt_commelem);

                                auto pnt1_elems = meshtopo.GetVertexSurfaceElements(segpair_p1);
                                auto pnt2_elems = meshtopo.GetVertexSurfaceElements(segpair_p2);

                                for(auto pnt1_sei : pnt1_elems)
                                  if(mesh[pnt1_sei].GetIndex() == mesh[sej].si)
                                    for(auto pnt2_sei : pnt2_elems)
                                      if(pnt1_sei == pnt2_sei)
                                        pnt_commelem = pnt1_sei;

                                if(IsInvalid(pnt_commelem))
                                  throw Exception("Couldn't find element on other side for " + ToString(segpair_p1) + " to " + ToString(segpair_p2));

                                const auto& commsel = mesh[pnt_commelem];
                                auto surfelem_vect = GetSurfaceNormal(mesh, commsel);
                                if(blp.outside)
                                  surfelem_vect *= -1;
                                Element2d sel(QUAD);
                                if(blp.outside)
                                  Swap(seg_p1, seg_p2);
                                sel[0] = seg_p1;
                                sel[1] = seg_p2;
                                sel[2] = mapto[seg_p2];
                                sel[3] = mapto[seg_p1];
                                auto domains = make_tuple(commsel.GetIndex(), blp.new_matnrs[layer-1], mesh.GetFaceDescriptor(commsel.GetIndex()).DomainOut());

                                if(domains_to_surf_index.find(domains) == domains_to_surf_index.end())
                                  {
                                    domains_to_surf_index[domains] = ++max_surface_index;
                                    domains_to_surf_index[make_tuple(max_surface_index, get<1>(domains), get<2>(domains))] = max_surface_index;
                                    FaceDescriptor fd(max_surface_index-1,
                                                      get<1>(domains),
                                                      get<2>(domains),
                                                      -1);
                                    fd.SetBCProperty(max_surface_index);
                                    mesh.AddFaceDescriptor(fd);
                                    mesh.SetBCName(max_surface_index-1,
                                                   mesh.GetBCName(get<0>(domains)-1));
                                  }
                                auto new_index = domains_to_surf_index[domains];
                                sel.SetIndex(new_index);
                                mesh.AddSurfaceElement(sel);

                                // Add segments
                                Segment seg_1, seg_2;
                                seg_1[0] = mapto[seg_p1];
                                seg_1[1] = seg_p1;
                                seg_2[0] = seg_p2;
                                seg_2[1] = mapto[seg_p2];
                                auto points = make_tuple(seg_p1, mapto[seg_p1]);
                                if(pi_to_edgenr.find(points) == pi_to_edgenr.end())
                                  pi_to_edgenr[points] = ++max_edge_nr;
                                seg_1.edgenr = pi_to_edgenr[points];
                                seg_1[2] = PointIndex::INVALID;
                                seg_1.si = new_index;
                                mesh.AddSegment(seg_1);

                                points = make_tuple(seg_p2, mapto[seg_p2]);
                                if(pi_to_edgenr.find(points) == pi_to_edgenr.end())
                                  pi_to_edgenr[points] = ++max_edge_nr;

                                seg_2[2] = PointIndex::INVALID;
                                seg_2.edgenr = pi_to_edgenr[points];
                                seg_2.si = new_index;
                                mesh.AddSegment(seg_2);
                                mesh[sej].si = new_index;
                              }

                            // in last layer insert new segments
                            if(layer == blp.heights.Size())
                              {
                                Segment s1 = mesh[sei];
                                Segment s2 = mesh[sej];
                                s1.edgenr = ++max_edge_nr;
                                s2.edgenr = max_edge_nr;
                                bool create_it = true;
                                if(blp.surfid.Contains(mesh[sej].si))
                                  {
                                    if(last_layer_surface_index_map.find(s1.si) != last_layer_surface_index_map.end() &&
                                       last_layer_surface_index_map.find(s2.si) != last_layer_surface_index_map.end())
                                      // edge already mapped
                                      create_it = false;
                                    s2.si = map_surface_index(s2.si);
                                  }
                                else
                                  {
                                    s2.si = domains_to_surf_index[make_tuple(s2.si,
                                                                             blp.new_matnrs[layer-1], mesh.GetFaceDescriptor(s2.si).DomainOut())];
                                  }
                                s1.si = map_surface_index(s1.si);
                                if(create_it)
                                  {
                                    mesh.AddSegment(s1);
                                    mesh.AddSegment(s2);
                                  }
                              }

                            // remap the segments to the new points
                            mesh[sei][0] = mapto[mesh[sei][0]];
                            mesh[sei][1] = mapto[mesh[sei][1]];
                            mesh[sej][1] = mapto[mesh[sej][1]];
                            mesh[sej][0] = mapto[mesh[sej][0]];

                          }
                      }
                  }
                else
                  {
                    // check if it doesn't contain the other edge as well
                    // and if it doesn't contain both mark them as done and
                    // if necessary map them
                    for(SegmentIndex sej = 0; sej<nseg; sej++)
                      {
                        if(mesh[sej][0] == mesh[sei][1] &&
                           mesh[sej][1] == mesh[sei][0])
                          {
                            if(!blp.surfid.Contains(mesh[sej].si))
                              {
                                segsel.Clear(sei);
                                segsel.Clear(sej);
                                PointIndex mapped_point = PointIndex::INVALID;
                                auto p1 = mesh[sei][0];
                                auto p2 = mesh[sei][1];
                                if(mapto[p1].IsValid())
                                  mapped_point = p1;
                                else if(mapto[p2].IsValid())
                                  mapped_point = p2;
                                else
                                  continue;
                                auto other_point = mapped_point == p1 ? p2 : p1;
                                if(growthvectors[mapped_point] * (mesh[other_point] - mesh[mapped_point]) < 0)
                                  {
                                    if(mapto[mesh[sei][0]].IsValid())
                                      mesh[sei][0] = mapto[mesh[sei][0]];
                                    if(mapto[mesh[sei][1]].IsValid())
                                      mesh[sei][1] = mapto[mesh[sei][1]];
                                    if(mapto[mesh[sej][0]].IsValid())
                                      mesh[sej][0] = mapto[mesh[sej][0]];
                                    if(mapto[mesh[sej][1]].IsValid())
                                      mesh[sej][1] = mapto[mesh[sej][1]];
                                  }
                              }
                            else
                              continue;
                          }
                      }
                  }
                  }
              }

          // add surface elements between layer and old domain
          if(layer == blp.heights.Size())
            {
              for(SurfaceElementIndex si = 0; si < nse; si++)
                {
                  if(blp.surfid.Contains(mesh[si].GetIndex()))
                    {
                      const auto& sel = mesh[si];
                      Element2d newel = sel;
                      newel.SetIndex(map_surface_index(sel.GetIndex()));
                      mesh.AddSurfaceElement(newel);
                    }
                }
            }

          // Add prismatic cells at the boundaries
          PrintMessage(3, "Generating prism boundary layer volume elements...");

          for (SurfaceElementIndex si = 0; si < nse; si++)
            {
              const auto& sel = mesh[si];
              if(blp.surfid.Contains(sel.GetIndex()))
                {
                  int classify = 0;
                  for(auto j : Range(sel.PNums()))
                    if (mapto[sel[j]].IsValid())
                      classify += (1 << j);

                  if(classify == 0)
                    continue;

                  Element el;

                  if(sel.GetType() == TRIG)
                    {
                      ELEMENT_TYPE types[] = { PRISM, TET, TET, PYRAMID,
                                               TET, PYRAMID, PYRAMID, PRISM };
                      int nums[] = { sel[0], sel[1], sel[2], mapto[sel[0]], mapto[sel[1]], mapto[sel[2]] };
                      int vertices[][6] =
                        {
                         { 0, 1, 2, 0, 1, 2 },   // should not occur
                         { 0, 2, 1, 3, 0, 0 },
                         { 0, 2, 1, 4, 0, 0 },
                         { 0, 1, 4, 3, 2, 0 },

                         { 0, 2, 1, 5, 0, 0 },
                         { 2, 0, 3, 5, 1, 0 },
                         { 1, 2, 5, 4, 0, 0 },
                         { 0, 2, 1, 3, 5, 4 }
                        };
                      if(blp.outside)
                        {
                          if(classify != 7)
                            throw Exception("Outside with non prisms not yet implemented");
                          for(auto i : Range(6))
                            vertices[7][i] = i;
                        }

                      el = Element(types[classify]);
                      for(auto i : Range(el.PNums()))
                        el.PNums()[i] = nums[vertices[classify][i]];
                    }
                  else // sel.GetType() == QUAD
                    {
                      int nums[] = { sel[0], sel[1], sel[2], sel[3],
                                     mapto[sel[0]], mapto[sel[1]],
                                     mapto[sel[2]], mapto[sel[3]] };
                      if(classify == 15)
                        {
                          int vertices[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };
                          if(!blp.outside)
                            {
                              Swap(vertices[1], vertices[3]);
                              Swap(vertices[5], vertices[7]);
                            }
                          el = Element(HEX);
                          for(auto i : Range(el.PNums()))
                            el.PNums()[i] = nums[vertices[i]];
                        }
                      else
                        {
                          throw Exception("This type of quad layer not yet implemented!");
                        }
                    }
                  el.SetIndex(blp.new_matnrs[layer-1]);
                  mesh.AddVolumeElement(el);
                }
            }

          // Finally switch the point indices of the surface elements
          // to the newly added ones
          PrintMessage(3, "Transferring boundary layer surface elements to new vertex references...");

          for(SurfaceElementIndex sei : Range(nse))
            {
              auto& sel = mesh[sei];
              bool to_move = blp.surfid.Contains(sel.GetIndex());
              if(blp.domains.Size())
                {
                  if(blp.outside)
                    to_move |= blp.domains[mesh.GetFaceDescriptor(sel.GetIndex()).DomainIn()];
                  else
                    to_move |= !blp.domains[mesh.GetFaceDescriptor(sel.GetIndex()).DomainIn()];
                }

              if(to_move)
                {
                  for(auto& pnum : sel.PNums())
                    // Check (Doublecheck) if the corresponding point has a
                    // copy available for remapping
                    if(mapto[pnum].IsValid())
                      // Map the surface elements to the new points
                      pnum = mapto[pnum];
                }
            }

          for(ElementIndex ei : Range(ne))
            {
              auto& el = mesh[ei];
              bool to_move = blp.outside ? blp.domains[el.GetIndex()] : !blp.domains[el.GetIndex()];
              if(blp.domains.Size() == 0 || to_move)
                for(auto& pnum : el.PNums())
                  // Check (Doublecheck) if the corresponding point has a
                  // copy available for remapping
                  if(mapto[pnum].IsValid())
                    // Map the volume elements to the new points
                    pnum = mapto[pnum];
            }

          // Lock all the prism points so that the rest of the mesh can be
          // optimised without invalidating the entire mesh
          // for (PointIndex pi = mesh.Points().Begin(); pi < mesh.Points().End(); pi++)
          for (PointIndex pi = 1; pi <= np; pi++)
            if(bndnodes.Test(pi)) mesh.AddLockedPoint(pi);

          // Now, actually pull back the old surface points to create
          // the actual boundary layers
          PrintMessage(3, "Moving and optimising boundary layer points...");

          for (PointIndex i = 1; i <= np; i++)
            {
              if(bndnodes.Test(i))
                {
                  MeshPoint pointtomove;
                  pointtomove = mesh.Point(i);
                  mesh.Point(i).SetPoint(pointtomove + layerht * growthvectors[i]);
                }
            }
          mesh.Compress();
        }

      for(int i=1; i <= mesh.GetNFD(); i++)
        {
          auto& fd = mesh.GetFaceDescriptor(i);
          if(blp.surfid.Contains(fd.BCProperty()))
            {
              if(blp.outside)
                fd.SetDomainOut(blp.new_matnrs[blp.new_matnrs.Size()-1]);
              else
                fd.SetDomainIn(blp.new_matnrs[blp.new_matnrs.Size()-1]);
            }
        }

      PrintMessage(3, "New NP: ", mesh.GetNP());
      PrintMessage(1, "Boundary Layer Generation....Done!");
   }
}
