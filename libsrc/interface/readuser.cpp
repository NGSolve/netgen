//
//  Read user dependent output file
//


#include <mystdlib.h>
#include <myadt.hpp>
#include <linalg.hpp>
#include <csg.hpp>
#include <stlgeom.hpp>
#include <meshing.hpp>

#include "writeuser.hpp"

namespace netgen
{
  void ReadFile (Mesh & mesh,
                 const string & hfilename)
  {
    PrintMessage(3, "Read User File");

    const char * filename = hfilename.c_str();

    char reco[100];
    int np, nbe;



    // ".surf" - mesh
  
    if ( (strlen (filename) > 5) &&
         strcmp (&filename[strlen (filename)-5], ".surf") == 0 )
    
      {
        cout << "Surface file" << endl;
      
        ifstream in (filename);
      
        in >> reco;
        in >> np;
        for (int i = 1; i <= np; i++)
          {
            Point3d p;
            in >> p.X() >> p.Y() >> p.Z();
            mesh.AddPoint (p);
          }

        mesh.ClearFaceDescriptors();
        mesh.AddFaceDescriptor (FaceDescriptor(1,1,0,0));
      
        in >> nbe;
        //      int invert = globflags.GetDefineFlag ("invertsurfacemesh");
        for (int i = 1; i <= nbe; i++)
          {
            Element2d el;
            el.SetIndex(1);

            for (int j = 1; j <= 3; j++)
              {
                in >> el.PNum(j);
                // el.PNum(j)++;
                if (el.PNum(j) < PointIndex(1) || 
                    el.PNum(j) > PointIndex(np))
                  {
                    cerr << "Point Number " << el.PNum(j) << " out of range 1..."
                         << np << endl;
                    return;
                  }
              }
            /*
              if (invert)
              swap (el.PNum(2), el.PNum(3));
            */
	  
            mesh.AddSurfaceElement (el);
          }
      
      
        cout << "points: " << np << " faces: " << nbe << endl;
      }
  
  
  

  
    if ( (strlen (filename) > 4) &&
         strcmp (&filename[strlen (filename)-4], ".unv") == 0 )
      {  
        char reco[100];
        // int invert;
      
        ifstream in(filename);

        mesh.ClearFaceDescriptors();
        mesh.AddFaceDescriptor (FaceDescriptor(0,1,0,0));
        mesh.GetFaceDescriptor(1).SetBCProperty (1);
        // map from unv element nr to our element number + an index if it is vol (0), bnd(1), ...
        std::map<size_t, std::tuple<size_t, int>> element_map;
	int dim = 3;
	int bccounter = 0;

        NgArray<Segment> tmp_segments;
        while (in.good())
          {
            in >> reco;
            if (strcmp(reco, "-1") == 0)
              continue;

            else if (strcmp (reco, "2411") == 0)
              {
                cout << "nodes found" << endl;

                while (1)
                  {
                    int pi, hi;
                    Point<3> p;

                    in >> pi;
                    if (pi == -1)
                      break;
		    
                    in >> hi >> hi >> hi;
                    in >> p(0) >> p(1) >> p(2);

                    mesh.AddPoint (p);
                  }
		cout << "read " << mesh.GetNP() << " points" << endl;
                Point3d pmin, pmax;
                mesh.GetBox (pmin, pmax);
                if(fabs(pmin.Z() - pmax.Z()) < 1e-10 * Dist(pmin, pmax))
                {
                       cout << "Set Dimension to 2." << endl;
                       mesh.SetDimension(2);
                       dim = 2 ;
		}

              }

            else if (strcmp (reco, "2412") == 0)
              {
                cout << "elements found" << endl;

                while (1)
                  {
		    int label, fe_id, phys_prop, mat_prop, color, nnodes;
		    int nodes[100];
		    int hi;

		    in >> label;
		    if (label == -1) break;
		    in >> fe_id >> phys_prop >> mat_prop >> color >> nnodes;
		    
		    if (fe_id >= 11 && fe_id <= 32)
		      in >> hi >> hi >> hi;
		      

		    for (int j = 0; j < nnodes; j++)
		      in >> nodes[j];
		    
		    switch (fe_id)
		      {
                      case 11: // (Rod) SEGM
                        {
                          Segment el;
                          el[0] = nodes[0];
                          el[1] = nodes[1];
                          el[2] = -1;

                          if(dim == 3){
                            auto nr = tmp_segments.Size();
                            tmp_segments.Append(el);
                            element_map[label] = std::make_tuple(nr+1, 2);
                          }
                          else if(dim == 2){
		            el.si = -1; // add label to segment, will be changed later when BC's are assigned
                            auto nr = mesh.AddSegment(el);
                            element_map[label] = std::make_tuple(nr+1, 2);
                          }
                          break;
                        }

                      case 22: // (Tapered beam) SEGM
                        {
                          Segment el;
                          el[0] = nodes[0];
                          el[1] = nodes[2];
                          el[2] = nodes[1];
                          
                          if(dim == 3){
                            auto nr = tmp_segments.Size();
                            tmp_segments.Append(el);
                            element_map[label] = std::make_tuple(nr+1, 2);
                          }
                          else if(dim == 2){
		            el.si = -1; // add label to segment, will be changed later when BC's are assigned
                            auto nr = mesh.AddSegment(el);
                            element_map[label] = std::make_tuple(nr+1, 2);
                          }

                          break;
                        }
		      case 41: // TRIG
			{
			  Element2d el (TRIG);
			  el.SetIndex (1);
			  for (int j = 0; j < nnodes; j++)
			    el[j] = nodes[j];
			  auto nr = mesh.AddSurfaceElement (el);
                          element_map[label] = std::make_tuple(nr+1, 1);
			  break;
			}
                      case 42: // TRIG6
                        {
                          Element2d el(TRIG6);
                          el.SetIndex(1);
                          int jj = 0;
                          for(auto j : {0,2,4,3,5,1})
                              el[jj++] = nodes[j];
                          auto nr = mesh.AddSurfaceElement(el);
                          element_map[label] = std::make_tuple(nr+1, 1);
                          break;
                        }
		      case 111: // TET
			{
			  Element el (TET);
			  el.SetIndex (1);
			  for (int j = 0; j < nnodes; j++)
			    el[j] = nodes[j];
			  auto nr = mesh.AddVolumeElement (el);
			  element_map[label] = std::make_tuple(nr+1, 0);
			  break;
			}
                      case 118: // TET10
                        {
                          Element el(TET10);
                          el.SetIndex(1);
                          int jj = 0;
                          for(auto j : {0,2,4,9,1,5,6,3,7,8})
                            el[jj++] = nodes[j];
                          auto nr = mesh.AddVolumeElement(el);
                          element_map[label] = std::make_tuple(nr+1, 0);
                          break;
                        }
                      default:
                        cout << "Do not know fe_id = " << fe_id << ", skipping it." << endl;
                        break;
		      }
                  }
                cout << mesh.GetNE() << " elements found" << endl;
                cout << mesh.GetNSE() << " surface elements found" << endl;

              }
            else if(strcmp (reco, "2467") == 0)
              {
                int matnr = 1;
                cout << "Groups found" << endl;
                while(in.good())
                  {
                    int len;
                    string name;
                    in >> len;
                    if(len == -1)
                      break;
                    for(int i=0; i < 7; i++)
                      in >> len;
                    in >> name;
                    cout << len << " element are in group " << name << endl;
                    int hi, index;
                    int fdnr, ednr;

                    in >> hi >> index >> hi >> hi;
                    int codim = get<1>(element_map[index]);
                    // use first element to determine if boundary or volume
                    
                    switch (codim)
                      {
                      case 0:
                        {
                          mesh.SetMaterial(++matnr, name);
                          mesh.VolumeElement(get<0>(element_map[index])).SetIndex(matnr);
                          break;
                        }
                      case 1:
                        {
                          if(dim == 3)
                          {
                            int bcpr = mesh.GetNFD();
                            fdnr = mesh.AddFaceDescriptor(FaceDescriptor(bcpr, 0,0,0));
                            mesh.GetFaceDescriptor(fdnr).SetBCProperty(bcpr+1);
                            mesh.SetBCName(bcpr, name);
                            mesh.SurfaceElement(get<0>(element_map[index])).SetIndex(fdnr);
                            bccounter++;
                          }
                          else if(dim == 2)
                          {
                            mesh.SetMaterial(matnr, name);
                            fdnr = mesh.AddFaceDescriptor(FaceDescriptor(matnr, 0,0,0));
                            mesh.SurfaceElement(get<0>(element_map[index])).SetIndex(matnr);
                            mesh.GetFaceDescriptor(fdnr).SetBCProperty(matnr);
			    matnr++;
                          }
                          break;

                        }
                      case 2:
                        {
                         if(dim == 3)
                          {
                            int bcpr = mesh.GetNCD2Names()+1;
                            auto ed = EdgeDescriptor();
                            ed.SetSurfNr(0,bcpr);//?
                            ednr = mesh.AddEdgeDescriptor(ed);
                            mesh.SetCD2Name(bcpr, name);
                            auto nr = mesh.AddSegment(tmp_segments[get<0>(element_map[index])-1]);
                            mesh[nr].edgenr = ednr+1;
                          }
                          else if(dim == 2)
                          {
                            Segment & seg = mesh.LineSegment(get<0>(element_map[index]));
			    seg.si = bccounter + 1;
			    mesh.SetBCName(bccounter, name);
		            bccounter++;
                          }
                          break;

                        }
                      default:
                        {
                          cout << "Codim " << codim << " not implemented yet!" << endl;
                        }
                      }
                        
                    for(int i=0; i<len-1; i++)
                      {
                        in >> hi >> index >> hi >> hi;
                        switch (codim)
                          {
                          case 0:
                            mesh.VolumeElement(get<0>(element_map[index])).SetIndex(matnr);
                            break;
                          case 1:
			    if(dim == 3) mesh.SurfaceElement(get<0>(element_map[index])).SetIndex(fdnr);
			    else if (dim == 2){
                                    mesh.SurfaceElement(get<0>(element_map[index])).SetIndex(matnr-1);
				    mesh.GetFaceDescriptor(fdnr).SetBCProperty(matnr);
			    }
                            break;
                          case 2:
	   		    if(dim == 3)
                            {
                              auto nr = mesh.AddSegment(tmp_segments[get<0>(element_map[index])-1]);
                              mesh[nr].edgenr = ednr+1;
                            }
			    else if(dim == 2)
			    {
	 			    Segment & seg = mesh.LineSegment(get<0>(element_map[index]));
			            seg.si = bccounter;
			    }
                            break;
                          default:
                            break;
                          }
                      }
                  }
              }
            else
              {
                cout << "Do not know data field type " << reco << ", skipping it" << endl;
                while(in.good())
                  {
                    in >> reco;
                    if(strcmp(reco, "-1") == 0)
                      break;
                  }
              }
          }

	if(dim == 2){
		// loop through segments to assign default BC to unmarked edges
		int bccounter_tmp = bccounter;
		for(int index=1; index <= mesh.GetNSeg(); index++){
                	Segment & seg = mesh.LineSegment(index);
			if(seg.si == -1){
			  seg.si = bccounter + 1;
			  if(bccounter_tmp == bccounter) mesh.SetBCName(bccounter, "default"); // could be more efficient
			  bccounter_tmp++;
			}
		}
		if(bccounter_tmp > bccounter) bccounter++;
	}
      

        Point3d pmin, pmax;
        mesh.ComputeNVertices();
        mesh.RebuildSurfaceElementLists();
        mesh.GetBox (pmin, pmax);
        mesh.UpdateTopology();
        if(dim == 3) bccounter++;
        cout << "bounding-box = " << pmin << "-" << pmax << endl;
	cout << "Created " << bccounter << " boundaries." << endl;
	for(int i=0; i<bccounter; i++){
		cout << mesh.GetBCName(i) << endl;
	}
      }



    // fepp format2d:
  
    if ( (strlen (filename) > 7) &&
         strcmp (&filename[strlen (filename)-7], ".mesh2d") == 0 )
      {
        cout << "Reading FEPP2D Mesh" << endl;
      
        char buf[100];
        int np, ne, nseg, i, j;

        ifstream in (filename);

        in >> buf;

        in >> nseg;
        for (i = 1; i <= nseg; i++)
          {
            int bound, p1, p2;
            in >> bound >> p1 >> p2;
            // forget them
          }

        in >> ne;
        for (i = 1; i <= ne; i++)
          {
            int mat, nelp;
            in >> mat >> nelp;
            Element2d el (nelp == 3 ? TRIG : QUAD);
            el.SetIndex (mat);
            for (j = 1; j <= nelp; j++)
              in >> el.PNum(j);
            mesh.AddSurfaceElement (el);
          }

        in >> np;
        for (i = 1; i <= np; i++)
          {
            Point3d p(0,0,0);
            in >> p.X() >> p.Y();
            mesh.AddPoint (p);
          }
      }

  
    else if ( (strlen (filename) > 5) &&
              strcmp (&filename[strlen (filename)-5], ".mesh") == 0 )
      {
        cout << "Reading Neutral Format" << endl;
      
        int np, ne, nse, i, j;

        ifstream in (filename);

        in >> np;

        if (in.good())
          {
            // file starts with an integer

            for (i = 1; i <= np; i++)
              {
                Point3d p(0,0,0);
                in >> p.X() >> p.Y() >> p.Z();
                mesh.AddPoint (p);
              }
	  
            in >> ne;
            for (i = 1; i <= ne; i++)
              {
                int mat;
                in >> mat;
                Element el (4);
                el.SetIndex (mat);
                for (j = 1; j <= 4; j++)
                  in >> el.PNum(j);
                mesh.AddVolumeElement (el);
              }

            mesh.AddFaceDescriptor (FaceDescriptor (1, 1, 0, 0));
            int nfd = 1;

            in >> nse;
            for (i = 1; i <= nse; i++)
              {
                int mat; // , nelp;
                in >> mat;
                Element2d el (TRIG);
                el.SetIndex (mat);
                while(nfd<mat)
                  {
                    ++nfd;
                    mesh.AddFaceDescriptor(FaceDescriptor(nfd,nfd,0,0));
                  }
                for (j = 1; j <= 3; j++)
                  in >> el.PNum(j);
                mesh.AddSurfaceElement (el);
              }
          }
        else
          {
            char buf[100];
            in.clear();
            do
              {
                in >> buf;
                cout << "buf = " << buf << endl;
                if (strcmp (buf, "points") == 0)
                  {
                    in >> np;
                    cout << "np = " << np << endl;
                  }
              }
            while (in.good());
          }
      }


    if ( (strlen (filename) > 4) &&
         strcmp (&filename[strlen (filename)-4], ".emt") == 0 )
      {
        ifstream inemt (filename);
      
        string pktfile = filename;
        int len = strlen (filename);
        pktfile[len-3] = 'p';
        pktfile[len-2] = 'k';
        pktfile[len-1] = 't';
        cout << "pktfile = " << pktfile << endl;

        int np, nse, i;
        int bcprop;
        ifstream inpkt (pktfile.c_str());
        inpkt >> np;
        NgArray<double> values(np);
        for (i = 1; i <= np; i++)
          {
            Point3d p(0,0,0);
            inpkt >> p.X() >> p.Y() >> p.Z()
                  >> bcprop >> values.Elem(i);
            mesh.AddPoint (p);
          }      

        mesh.ClearFaceDescriptors();
        mesh.AddFaceDescriptor (FaceDescriptor(0,1,0,0));
        mesh.GetFaceDescriptor(1).SetBCProperty (1);
        mesh.AddFaceDescriptor (FaceDescriptor(0,1,0,0));
        mesh.GetFaceDescriptor(2).SetBCProperty (2);
        mesh.AddFaceDescriptor (FaceDescriptor(0,1,0,0));
        mesh.GetFaceDescriptor(3).SetBCProperty (3);
        mesh.AddFaceDescriptor (FaceDescriptor(0,1,0,0));
        mesh.GetFaceDescriptor(4).SetBCProperty (4);
        mesh.AddFaceDescriptor (FaceDescriptor(0,1,0,0));
        mesh.GetFaceDescriptor(5).SetBCProperty (5);

        int p1, p2, p3;
        double value;
        inemt >> nse;
        for (i = 1; i <= nse; i++)
          {
            inemt >> p1 >> p2 >> p3 >> bcprop >> value;

            if (bcprop < 1 || bcprop > 4)
              cerr << "bcprop out of range, bcprop = " << bcprop << endl;
            p1++;
            p2++;
            p3++;
            if (p1 < 1 || p1 > np || p2 < 1 || p2 > np || p3 < 1 || p3 > np)
              {
                cout << "p1 = " << p1 << " p2 = " << p2 << " p3 = " << p3 << endl;
              }

            if (i > 110354) Swap (p2, p3);
            if (mesh.Point(p1)(0) < 0.25)
              Swap (p2,p3);

            Element2d el(TRIG);

            if (bcprop == 1)
              {
                if (values.Get(p1) < -69999)
                  el.SetIndex(1);
                else
                  el.SetIndex(2);
              }
            else
              el.SetIndex(3);


            el.PNum(1) = p1;
            el.PNum(2) = p2;
            el.PNum(3) = p3;
            mesh.AddSurfaceElement (el);
          }


        ifstream incyl ("ngusers/guenter/cylinder.surf");
        int npcyl, nsecyl; 
        incyl >> npcyl;
        cout << "npcyl = " << npcyl << endl;
        for (i = 1; i <= npcyl; i++)
          {
            Point3d p(0,0,0);
            incyl >> p.X() >> p.Y() >> p.Z();
            mesh.AddPoint (p);
          }
        incyl >> nsecyl;
        cout << "nsecyl = " << nsecyl << endl;
        for (i = 1; i <= nsecyl; i++)
          {
            incyl >> p1 >> p2 >> p3;
            p1 += np;
            p2 += np;
            p3 += np;
            Element2d el(TRIG);
            el.SetIndex(5);
            el.PNum(1) = p1;
            el.PNum(2) = p2;
            el.PNum(3) = p3;
            mesh.AddSurfaceElement (el);
          }
      }


    // .tet mesh
    if ( (strlen (filename) > 4) &&
         strcmp (&filename[strlen (filename)-4], ".tet") == 0 )
      {
        ReadTETFormat (mesh, filename);
      }


    // .fnf mesh (FNF - PE neutral format)
    if ( (strlen (filename) > 4) &&
         strcmp (&filename[strlen (filename)-4], ".fnf") == 0 )
      {
        ReadFNFFormat (mesh, filename);
      }

    // .cgns file - CFD General Notation System
    if ( (strlen (filename) > 5) &&
         strcmp (&filename[strlen (filename)-5], ".cgns") == 0 )
      {
        ReadCGNSMesh (mesh, filename);
      }

    if ( ( (strlen (filename) > 4) && strcmp (&filename[strlen (filename)-4], ".stl") == 0 ) ||
         ( (strlen (filename) > 5) && strcmp (&filename[strlen (filename)-5], ".stlb") == 0 ) )
      {
        ifstream ist{string{filename}};
        auto geom = shared_ptr<STLGeometry>(STLGeometry::Load(ist));

        mesh.SetDimension (3);

        auto & points = geom->GetPoints();

        for (auto & p : points)
          mesh.AddPoint(MeshPoint(p));

        mesh.AddFaceDescriptor (FaceDescriptor (1, 1, 0, 1));

        for (auto ti : IntRange(geom->GetNT()))
        {
          Element2d el(TRIG);
          for (auto i : IntRange(3))
            el[i] = int((*geom)[STLTrigId(ti+IndexBASE<netgen::STLTrigId>())][i]);

          el.SetIndex(1);

          mesh.AddSurfaceElement(el);
        }
      }
  }
  
}

