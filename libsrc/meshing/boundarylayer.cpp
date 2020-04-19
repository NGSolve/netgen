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
  Vec<3> GetSurfaceNormal(Mesh & mesh, const Element2d & el)
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
      // Angle between a surface element and a growth-vector below which
      // a prism is project onto that surface as a quad
      // (in degrees)
      double angleThreshold = 5.0;

      // Monitor and print out the number of prism and quad elements
      // added to the mesh
      int numprisms = 0;
      int numquads = 0;

      PrintMessage(1, "Generating boundary layer...");
      PrintMessage(3, "Old NP: ", mesh.GetNP());
      PrintMessage(3, "Old NSE: ",mesh.GetNSE());

      for(int layer = blp.heights.Size(); layer >= 1; layer--)
        {
          PrintMessage(3, "Generating layer: ", layer);

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
          Array<Array<Vec<3>>, PointIndex> all_growthvectors(np);

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
                auto normal = GetSurfaceNormal(mesh,sel);
                if(!blp.outside)
                  normal *= -1;
                for(int j : Range(sel.PNums()))
                  {
                    // Set the bitarray to indicate that the
                    // point is part of the required set
                    bndnodes.SetBit(sel[j]);

                    // Add the surface normal to the already existent one
                    // (This gives the effective normal direction at corners
                    //  and curved areas)
                    all_growthvectors[sel[j]].Append(normal);
                  }
              }

          if (!blp.grow_edges)
            for(const auto& sel : mesh.SurfaceElements())
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
                {
                  mapto[pi] = mesh.AddPoint(mesh[pi]);
                  growthvectors[pi] = all_growthvectors[pi][0];
                  for(int i = 1; i < all_growthvectors[pi].Size(); i++)
                    {
                      auto& veci = all_growthvectors[pi][i];
                      for(auto j : Range(i))
                        {
                          auto& vecj = all_growthvectors[pi][j];
                          veci -= 1./(vecj.Length()+1e-10) * (veci * vecj) * vecj;
                        }
                      growthvectors[pi] += veci;
                    }
                  // growthvectors[pi].Normalize();
                  // growthvectors[pi] *= -1.0;
                }
              else
                {
                  mapto[pi].Invalidate();
                  growthvectors[pi] = {0,0,0};
                }
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
                if(segsel.Test(sei) && blp.surfid.Contains(mesh[sei].si))
                  {
                    // clear the bit to indicate that this segment has been processed
                    segsel.Clear(sei);

                    // Find matching segment pair on other surface
                    for (SegmentIndex sej = 0; sej < nseg; sej++)
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
                                  {
                                    const auto& pnt1_sel = mesh[pnt1_sei];
                                    for(auto pnt2_sei : pnt2_elems)
                                      {
                                        const Element2d & pnt2_sel = mesh.SurfaceElement(pnt2_sei);
                                        if((pnt1_sel.GetIndex() == mesh[sej].si)
                                           && (pnt2_sel.GetIndex() == mesh[sej].si)
                                           && (pnt1_sei == pnt2_sei))
                                          {
                                            pnt_commelem = pnt1_sei;
                                          }
                                      }
                                  }

                                const Element2d & commsel = mesh.SurfaceElement(pnt_commelem);
                                auto surfelem_vect = GetSurfaceNormal(mesh, commsel);
                                if(blp.outside)
                                  surfelem_vect *= -1;

                                double surfangle = Angle(growthvectors[segpair_p1],surfelem_vect);
                                // remap the segments to the new points
                                if(!blp.outside)
                                  {
                                    mesh[sei][0] = mapto[seg_p1];
                                    mesh[sei][1] = mapto[seg_p2];
                                    mesh[sej][1] = mapto[seg_p1];
                                    mesh[sej][0] = mapto[seg_p2];
                                  }

                                if((surfangle < (90 + angleThreshold) * 3.141592 / 180.0)
                                   && (surfangle > (90 - angleThreshold) * 3.141592 / 180.0))
                                  {
                                    // Since the surface is lower than the threshold, change the effective
                                    // prism growth vector to match with the surface vector, so that
                                    // the Quad which is created lies on the original surface
                                    //growthvectors.Elem(segpair_p1) = surfelem_vect;

                                    // Add a quad element to account for the prism volume
                                    // element which is going to be added
                                    Element2d sel(QUAD);
                                    if(blp.outside)
                                      Swap(seg_p1, seg_p2);
                                    sel.PNum(4) = mapto[seg_p1];
                                    sel.PNum(3) = mapto[seg_p2];
                                    sel.PNum(2) = seg_p2;
                                    sel.PNum(1) = seg_p1;
                                    sel.SetIndex(mesh[sej].si);
                                    mesh.AddSurfaceElement(sel);
                                    numquads++;
                                  }
                                else
                                  {
                                    for (int k = 0; k < pnt1_elems.Size(); k++)
                                      {
                                        Element2d & pnt_sel = mesh.SurfaceElement(pnt1_elems[k]);
                                        if(pnt_sel.GetIndex() == mesh[sej].si)
                                          {
                                            for(int l = 0; l < pnt_sel.GetNP(); l++)
                                              {
                                                if(pnt_sel[l] == segpair_p1)
                                                  pnt_sel[l] = mapto[seg_p1];
                                                else if (pnt_sel[l] == segpair_p2)
                                                  pnt_sel[l] = mapto[seg_p2];
                                              }
                                          }
                                      }

                                    for (int k = 0; k < pnt2_elems.Size(); k++)
                                      {
                                        Element2d & pnt_sel = mesh.SurfaceElement(pnt2_elems[k]);
                                        if(pnt_sel.GetIndex() == mesh[sej].si)
                                          {
                                            for(int l = 0; l < pnt_sel.GetNP(); l++)
                                              {
                                                if(pnt_sel[l] == segpair_p1)
                                                  pnt_sel[l] = mapto[seg_p1];
                                                else if (pnt_sel[l] == segpair_p2)
                                                  pnt_sel[l] = mapto[seg_p2];
                                              }
                                          }
                                      }
                                  }
                                // }
                              }
                            else
                              {
                                // If the code comes here, it indicates that we are at
                                // a line segment pair which is at the intersection
                                // of two surfaces, both of which have to grow boundary
                                // layers.... here too, remapping the segments to the
                                // new points is required
                                mesh[sei][0] = mapto[seg_p1];
                                mesh[sei][1] = mapto[seg_p2];
                                mesh[sej][1] = mapto[seg_p1];
                                mesh[sej][0] = mapto[seg_p2];
                              }
                          }
                      }
                  }
              }

          // Add prismatic cells at the boundaries
          PrintMessage(3, "Generating prism boundary layer volume elements...");

          for (SurfaceElementIndex si = 0; si < nse; si++)
            {
              Element2d & sel = mesh.SurfaceElement(si);
              if(blp.surfid.Contains(sel.GetIndex()))
                {
                  int classify = 0;
                  for (int j = 0; j < 3; j++)
                    if (mapto[sel[j]].IsValid())
                      classify += (1 << j);

                  // cout << "classify = " << classify << endl;

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

                  Element el(types[classify]);
                  for (int i = 0; i < 6; i++)
                    el[i] = nums[vertices[classify][i]];
                  el.SetIndex(blp.new_matnrs[layer-1]);
                  if (classify != 0)
                    mesh.AddVolumeElement(el);
                }
            }

          // Finally switch the point indices of the surface elements
          // to the newly added ones
          PrintMessage(3, "Transferring boundary layer surface elements to new vertex references...");

          for(SurfaceElementIndex sei : Range(nse))
            {
              auto& sel = mesh[sei];
              if((blp.outside && !blp.surfid.Contains(sel.GetIndex())) ||
                 (!blp.outside && blp.surfid.Contains(sel.GetIndex())))
                {
                  for(auto& pnum : sel.PNums())
                    // Check (Doublecheck) if the corresponding point has a
                    // copy available for remapping
                    if(mapto[pnum].IsValid())
                      // Map the surface elements to the new points
                      pnum = mapto[pnum];
                }
            }

          if(blp.outside)
            for(ElementIndex ei : Range(ne))
              {
                auto& el = mesh[ei];
                for(auto& pnum : el.PNums())
                  // Check (Doublecheck) if the corresponding point has a
                  // copy available for remapping
                  if(mapto[pnum].IsValid())
                    // Map the surface elements to the new points
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
          if(blp.surfid.Contains(fd.SurfNr()))
            {
              if(blp.outside)
                fd.SetDomainOut(blp.new_matnrs[0]);
              else
                fd.SetDomainIn(blp.new_matnrs[0]);
            }
        }

      PrintMessage(3, "New NP: ", mesh.GetNP());
      PrintMessage(3, "Num of Quads: ", numquads);
      PrintMessage(3, "Num of Prisms: ", numprisms);
      PrintMessage(1, "Boundary Layer Generation....Done!");
   }
}
