
#include <mystdlib.h>
#include "meshing.hpp"


namespace netgen
{
  void GetPureBadness(Mesh & mesh, Array<double, PointIndex> & pure_badness,
                      const TBitArray<PointIndex> & isnewpoint)
  {
    //const int ne = mesh.GetNE();
    const int np = mesh.GetNP();

    pure_badness.SetSize(np+1);   // one extra slot for the maximum
    pure_badness = -1;

    Array<Point<3>, PointIndex> backup(np);

    for (PointIndex pi : mesh.Points().Range())
      {
        backup[pi] = mesh.Point(pi);

        if (isnewpoint.Test(pi) && mesh.mlbetweennodes[pi][0].IsValid())
          mesh.Point(pi) = Center (mesh.Point(mesh.mlbetweennodes[pi][0]),
                                   mesh.Point(mesh.mlbetweennodes[pi][1]));
      }
    for (auto & el : mesh.VolumeElements())
      {
        double bad = el.CalcJacobianBadness (mesh.Points());
        for(int j=0; j<el.GetNP(); j++)
          if(bad > pure_badness[el[j]])
            pure_badness[el[j]] = bad;

        // save maximum
        if(bad > pure_badness.Last())
          pure_badness.Last() = bad; 
      }
    
    for (PointIndex pi : mesh.Points().Range())
      mesh.Point(pi) = backup[pi];
  }


  double Validate(const Mesh & mesh, Array<ElementIndex> & bad_elements,
                  const Array<double, PointIndex> & pure_badness,
                  double max_worsening, const bool uselocalworsening,
                  Array<double, ElementIndex> * quality_loss)
  {
    PrintMessage(3,"!!!! Validating !!!!");
    //if(max_worsening > 0)
    //  (*testout) << "badness " << counter++ << endl;

    bad_elements.SetSize(0);

    double loc_pure_badness = -1;

    if(!uselocalworsening)
      loc_pure_badness = pure_badness.Last(); // maximum is saved at last position


    double worsening = -1;
    ElementIndex ind;

    if(quality_loss != NULL)
      quality_loss->SetSize(mesh.GetNE());

    for (ElementIndex i : mesh.VolumeElements().Range())
      {
        if(uselocalworsening)
          {
            loc_pure_badness = -1;
            for(int j=0; j<mesh[i].GetNP(); j++)
              if(pure_badness[mesh[i][j]] > loc_pure_badness)
                loc_pure_badness = pure_badness[mesh[i][j]];
          }


        double bad = mesh[i].CalcJacobianBadness (mesh.Points());
        if (bad > 1e10 || 
            (max_worsening > 0 && bad > loc_pure_badness*max_worsening))
          bad_elements.Append(i);
          

        if(max_worsening > 0)
          {
            double actw = bad/loc_pure_badness;
            if(quality_loss != NULL)
              (*quality_loss)[i] = actw;

            if(actw > worsening)
              {
                worsening = actw;
                ind = i;
              }
          }
      }
    return worsening;
  }


  void GetWorkingArea(TBitArray<ElementIndex> & working_elements, TBitArray<PointIndex> & working_points,
                      const Mesh & mesh, const Array<ElementIndex> & bad_elements,
                      const int width)
  {
    working_elements.Clear();
    working_points.Clear();

    for(int i=0; i<bad_elements.Size(); i++)
      {
        working_elements.SetBit(bad_elements[i]);
        const Element & el = mesh[bad_elements[i]];
        for (int j = 0; j < el.GetNP(); j++)
          working_points.SetBit(el[j]);
      }
    

    for(int i=0; i<width; i++)
      {
        for (ElementIndex j : mesh.VolumeElements().Range())
          {
            if(!working_elements.Test(j))
              {  
                const Element & el = mesh[j];
                bool set_active = false;
                
                for(int k=1; !set_active && k<=el.GetNP(); k++)
                  set_active = working_points.Test(el.PNum(k));
                
                if(set_active)
                  working_elements.SetBit(j);
              }
          }

        for (ElementIndex j : mesh.VolumeElements().Range())
          {
            if(working_elements.Test(j))
              {
                const Element & el = mesh[j];
                for (int k = 0; k < el.GetNP(); k++)
                  working_points.SetBit(el[k]);
              }
          }
      }
  }



  void RepairBisection(Mesh & mesh, Array<ElementIndex> & bad_elements, 
                       const TBitArray<PointIndex> & isnewpoint, const Refinement & refinement,
                       const Array<double, PointIndex> & pure_badness, 
                       double max_worsening, const bool uselocalworsening,
                       const Array< idmap_type* > & idmaps)
  {
    ostringstream ostrstr;

    const int maxtrials = 100;

    //bool doit;
    //cout << "DOIT: " << flush;
    //cin >> doit;

    int ne = mesh.GetNE();
    int np = mesh.GetNP();

    int numbadneighbours = 3;
    const int numtopimprove = 3;

    PrintMessage(1,"repairing");

    PushStatus("Repair Bisection");

    Array<Point<3>, PointIndex> should(np);
    Array<Point<3>, PointIndex> can(np);
    Array<Vec<3>* > nv(np);
    for(int i=0; i<np; i++)
      nv[i] = new Vec<3>;
    
    TBitArray<PointIndex> isboundarypoint(np),isedgepoint(np);
    isboundarypoint.Clear();
    isedgepoint.Clear();

    for (auto & seg : mesh.LineSegments())
      {
        isedgepoint.SetBit(seg[0]);
        isedgepoint.SetBit(seg[1]);
      }

    Array<int, PointIndex> surfaceindex(np);
    surfaceindex = -1;

    /*
    for (int i = 1; i <= mesh.GetNSE(); i++)
      {
        const Element2d & sel = mesh.SurfaceElement(i);
    */
    for (auto & sel : mesh.SurfaceElements())
      for (int j = 1; j <= sel.GetNP(); j++)
        if(!isedgepoint.Test(sel.PNum(j)))
          {
            isboundarypoint.SetBit(sel.PNum(j));
            surfaceindex[sel.PNum(j)] = 
              mesh.GetFaceDescriptor(sel.GetIndex()).SurfNr();
          }
    


    Validate(mesh,bad_elements,pure_badness,
             ((uselocalworsening) ?  (0.8*(max_worsening-1.) + 1.) : (0.1*(max_worsening-1.) + 1.)),
             uselocalworsening); // -> larger working area
    TBitArray<ElementIndex> working_elements(ne+1);
    TBitArray<PointIndex> working_points(np);

    GetWorkingArea(working_elements,working_points,mesh,bad_elements,numbadneighbours);
    //working_elements.Set();
    //working_points.Set();

    ostrstr.str("");
    ostrstr << "worsening: " <<
      Validate(mesh,bad_elements,pure_badness,max_worsening,uselocalworsening);
    PrintMessage(4,ostrstr.str());

    

    int auxnum=0;
    for(PointIndex pi : mesh.Points().Range())
      if(working_points.Test(pi))
        auxnum++;
    
    ostrstr.str("");
    ostrstr << "Percentage working points: " << 100.*double(auxnum)/np;
    PrintMessage(5,ostrstr.str());
    

    TBitArray<PointIndex> isworkingboundary(np);
    for (PointIndex pi : mesh.Points().Range())
      if(working_points.Test(pi) && isboundarypoint.Test(pi))
        isworkingboundary.SetBit(pi);
      else
        isworkingboundary.Clear(pi);


    for (PointIndex pi : mesh.Points().Range())
      should[pi] = mesh[pi];

    
    // for(int i=0; i<np; i++)
    for (PointIndex i = IndexBASE<PointIndex>(); i < IndexBASE<PointIndex>()+np; i++)
      {
        if(isnewpoint.Test(i) && 
           //working_points.Test(i+PointIndex::BASE) && 
           mesh.mlbetweennodes[i][0].IsValid())
          can[i] = Center(can[mesh.mlbetweennodes[i][0]],
                          can[mesh.mlbetweennodes[i][1]]);
        else
          can[i] = mesh[i];
      }


    int cnttrials = 1;
    
    double lamedge = 0.5;
    double lamface = 0.5;
    
    double facokedge = 0;
    double facokface = 0;
    double factryedge;
    double factryface = 0;

    double oldlamedge,oldlamface;

    auto geo = mesh.GetGeometry();
    if(!geo)
      {
        cerr << "No 2D Optimizer!" << endl;
        return;
      }    

    while ((facokedge < 1.-1e-8 || facokface < 1.-1e-8) && 
           cnttrials < maxtrials &&
           multithread.terminate != 1)
      {
        (*testout) << "   facokedge " << facokedge << " facokface " << facokface << " cnttrials " << cnttrials << endl
                   << " perc. " << 95. * max2( min2(facokedge,facokface),
                                               double(cnttrials)/double(maxtrials)) << endl;

        SetThreadPercent(95. * max2( min2(facokedge,facokface),
                                     double(cnttrials)/double(maxtrials)));

        ostrstr.str("");
        ostrstr << "max. worsening " << max_worsening;
        PrintMessage(5,ostrstr.str());
        oldlamedge = lamedge;
        lamedge *= 6;
        if (lamedge > 2)
          lamedge = 2;
           
        if(1==1 || facokedge < 1.-1e-8)
          {
            for(int i=0; i<nv.Size(); i++)
              *nv[i] = Vec<3>(0,0,0);
            /*
            for (int i = 1; i <= mesh.GetNSE(); i++)
              {
                const Element2d & sel = mesh.SurfaceElement(i);
            */
            for (auto & sel : mesh.SurfaceElements())
              {
                Vec<3> auxvec = Cross(mesh.Point(sel[1])-mesh.Point(sel[0]),
                                      mesh.Point(sel[2])-mesh.Point(sel[0]));
                auxvec.Normalize();
                for (int j = 0; j < sel.GetNP(); j++)
                  if(!isedgepoint.Test(sel[j]))
                    *nv[sel[j] - IndexBASE<PointIndex>()] += auxvec;
              }
            for(int i=0; i<nv.Size(); i++)
              nv[i]->Normalize();
            
            
            do  // move edges
              {
                lamedge *= 0.5;
                cnttrials++;
                if(cnttrials % 10 == 0)
                  max_worsening *= 1.1;
                
                
                factryedge = lamedge + (1.-lamedge) * facokedge;

                ostrstr.str("");
                ostrstr << "lamedge = " << lamedge << ", trying: " << factryedge;
                PrintMessage(5,ostrstr.str());
                

                for (PointIndex pi : mesh.Points().Range())
                  {
                    if (isedgepoint.Test(pi))
                      {
                        for (int j = 0; j < 3; j++)
                          mesh[pi](j) = 
                            lamedge * should[pi](j) +
                            (1.-lamedge) * can[pi](j);
                      }
                    else
                      mesh[pi] = can[pi];
                  }
                if(facokedge < 1.-1e-8)
                  {
                    ostrstr.str("");
                    ostrstr << "worsening: " <<
                      Validate(mesh,bad_elements,pure_badness,max_worsening,uselocalworsening);

                    PrintMessage(5,ostrstr.str());
                  }
                else
                  Validate(mesh,bad_elements,pure_badness,-1,uselocalworsening);


                ostrstr.str("");
                ostrstr << bad_elements.Size() << " bad elements";
                PrintMessage(5,ostrstr.str());
              }
            while (bad_elements.Size() > 0 && 
                   cnttrials < maxtrials &&
                   multithread.terminate != 1);
          }

        if(cnttrials < maxtrials &&
           multithread.terminate != 1)
          {
            facokedge = factryedge;
            
            // smooth faces
            mesh.CalcSurfacesOfNode();
            
            MeshingParameters dummymp;
            mesh.ImproveMeshJacobianOnSurface(dummymp,isworkingboundary,nv,OPT_QUALITY, &idmaps);
            
            for (PointIndex pi : mesh.Points().Range())
              can[pi] = mesh[pi];
            
            if(geo)
              for (PointIndex pi : surfaceindex.Range())
                {
                  if(surfaceindex[pi] >= 0)
                    {
                      should[pi] = can[pi];
                      geo->ProjectPoint(surfaceindex[pi],should[pi]);
                    }
                }
          }


        oldlamface = lamface;
        lamface *= 6;
        if (lamface > 2)
          lamface = 2;


        if(cnttrials < maxtrials &&
           multithread.terminate != 1)
          {

            do  // move faces
              {
                lamface *= 0.5;
                cnttrials++;
                if(cnttrials % 10 == 0)
                  max_worsening *= 1.1;
                factryface = lamface + (1.-lamface) * facokface;

                ostrstr.str("");
                ostrstr << "lamface = " << lamface << ", trying: " << factryface;
                PrintMessage(5,ostrstr.str());
                
                
                for (PointIndex pi : mesh.Points().Range())
                  {
                    if (isboundarypoint.Test(pi))
                      {
                        for (int j = 0; j < 3; j++)
                          mesh[pi](j) = 
                            lamface * should[pi](j) +
                            (1.-lamface) * can[pi](j);
                      }
                    else
                      mesh[pi] = can[pi];
                  }

                ostrstr.str("");
                ostrstr << "worsening: " <<
                  Validate(mesh,bad_elements,pure_badness,max_worsening,uselocalworsening);
                PrintMessage(5,ostrstr.str());
        

                ostrstr.str("");
                ostrstr << bad_elements.Size() << " bad elements";
                PrintMessage(5,ostrstr.str());
              }
            while (bad_elements.Size() > 0 && 
                   cnttrials < maxtrials &&
                   multithread.terminate != 1);
          }



        if(cnttrials < maxtrials &&
           multithread.terminate != 1)
          {
            facokface = factryface;
            // smooth interior
            
            mesh.CalcSurfacesOfNode();
            
            MeshingParameters dummymp;
            mesh.ImproveMeshJacobian (dummymp, OPT_QUALITY,&working_points);
            //mesh.ImproveMeshJacobian (OPT_WORSTCASE,&working_points);
          

            for (PointIndex pi : mesh.Points().Range())
              can[pi] = mesh[pi];
          }
          
        //!
        if((facokedge < 1.-1e-8 || facokface < 1.-1e-8) && 
           cnttrials < maxtrials &&
           multithread.terminate != 1)
          {
            MeshingParameters dummymp;
            MeshOptimize3d optmesh(mesh, dummymp, OPT_QUALITY);
            for(int i=0; i<numtopimprove; i++)
              {
                optmesh.SwapImproveSurface(&working_elements,&idmaps);
                optmesh.SwapImprove(&working_elements);
                
              }     

            //      mesh.mglevels = 1;
            
                
            ne = mesh.GetNE();
            working_elements.SetSize(ne);
            
            
            for (PointIndex pi : mesh.Points().Range())
              mesh[pi] = should[pi];
            
            Validate(mesh,bad_elements,pure_badness,
                     ((uselocalworsening) ?  (0.8*(max_worsening-1.) + 1.) : (0.1*(max_worsening-1.) + 1.)),
                     uselocalworsening);
            
            if(lamedge < oldlamedge || lamface < oldlamface)
              numbadneighbours++;
            GetWorkingArea(working_elements,working_points,mesh,bad_elements,numbadneighbours);
            for (PointIndex pi : mesh.Points().Range())
              if(working_points.Test(pi) && isboundarypoint.Test(pi))
                isworkingboundary.SetBit(pi);
              else
                isworkingboundary.Clear(pi);
            auxnum=0;
            for(PointIndex pi : mesh.Points().Range())
              if(working_points.Test(pi))
                auxnum++;

            
            ostrstr.str("");
            ostrstr << "Percentage working points: " << 100.*double(auxnum)/np;
            PrintMessage(5,ostrstr.str());
            
            for (PointIndex pi : mesh.Points().Range())
              mesh[pi] = can[pi];
          }
        //!

      }

    MeshingParameters dummymp;
    MeshOptimize3d optmesh(mesh, dummymp, OPT_QUALITY);
    for(int i=0; i<numtopimprove && multithread.terminate != 1; i++)
      {
        optmesh.SwapImproveSurface(NULL,&idmaps);
        optmesh.SwapImprove();
        //mesh.UpdateTopology();
      }
    mesh.UpdateTopology();
    /*
    if(cnttrials < 100)
      {
        nv = Vec<3>(0,0,0);
        for (int i = 1; i <= mesh.GetNSE(); i++)
          {
            const Element2d & sel = mesh.SurfaceElement(i);
            Vec<3> auxvec = Cross(mesh.Point(sel[1])-mesh.Point(sel[0]),
                                 mesh.Point(sel[2])-mesh.Point(sel[0]));
            auxvec.Normalize();
            for (int j = 1; j <= sel.GetNP(); j++)
              if(!isedgepoint.Test(sel.PNum(j)))
                nv[sel.PNum(j) - PointIndex::BASE] += auxvec;
          }
        for(int i=0; i<nv.Size(); i++)
          nv[i].Normalize();
        

        mesh.ImproveMeshJacobianOnSurface(isboundarypoint,nv,OPT_QUALITY);
        mesh.CalcSurfacesOfNode();
            // smooth interior
            
        
        for (int i = 1; i <= np; i++)
          if(isboundarypoint.Test(i))
            can.Elem(i) = mesh.Point(i);
            
        if(optimizer2d)
          optimizer2d->ProjectBoundaryPoints(surfaceindex,can,should);

        
        for (int i = 1; i <= np; i++)
          if(isboundarypoint.Test(i))
            for(int j=1; j<=3; j++)
              mesh.Point(i).X(j) = should.Get(i).X(j);
      }
    */


    if(cnttrials == maxtrials)
      {
        for (PointIndex pi : mesh.Points().Range())
          mesh[pi] = should[pi];

        Validate(mesh,bad_elements,pure_badness,max_worsening,uselocalworsening);
        
        for(int i=0; i<bad_elements.Size(); i++)
          {
            ostrstr.str("");
            ostrstr << "bad element:" << endl
                    << mesh[bad_elements[i]][0] << ": " << mesh.Point(mesh[bad_elements[i]][0]) << endl
                    << mesh[bad_elements[i]][1] << ": " << mesh.Point(mesh[bad_elements[i]][1]) << endl
                    << mesh[bad_elements[i]][2] << ": " << mesh.Point(mesh[bad_elements[i]][2]) << endl
                    << mesh[bad_elements[i]][3] << ": " << mesh.Point(mesh[bad_elements[i]][3]);
            PrintMessage(5,ostrstr.str());
          }
        for (PointIndex pi : mesh.Points().Range())
          mesh[pi] = can[pi];
      }

    for(int i=0; i<np; i++)
      delete nv[i];

    PopStatus();
  }
}
