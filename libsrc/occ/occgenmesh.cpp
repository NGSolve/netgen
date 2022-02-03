#ifdef OCCGEOMETRY

#include <mystdlib.h>
#include <meshing.hpp>

#include "occgeom.hpp"
#include "occmeshsurf.hpp"

#include <BRepAdaptor_Curve.hxx>
#include <BRepGProp.hxx>
#include <BRepLProp_CLProps.hxx>
#include <BRepLProp_SLProps.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRepTools.hxx>
#include <GProp_GProps.hxx>
#include <Quantity_Color.hxx>
#include <ShapeAnalysis.hxx>
#include <TopExp.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>
#include <TopoDS_Edge.hxx>

namespace netgen
{


#define TCL_OK 0
#define TCL_ERROR 1

#define DIVIDEEDGESECTIONS 10000   // better solution to come soon
#define IGNORECURVELENGTH 0
#define VSMALL 1e-10


  DLL_HEADER bool merge_solids = false;


  // can you please explain what you intend to compute here (JS) !!!
  double Line :: Dist (Line l)
  {
    Vec<3> n = p1-p0;
    Vec<3> q = l.p1-l.p0;
    double nq = n*q;

    Point<3> p = p0 + 0.5*n;
    double lambda = (p-l.p0)*n / (nq + VSMALL);

    if (lambda >= 0 && lambda <= 1)
      {
        double d = (p-l.p0-lambda*q).Length();
        //        if (d < 1e-3) d = 1e99;
        return d;
      }
    else
      return 1e99;
  }


  double ComputeH (double kappa, const MeshingParameters & mparam)
  {
    kappa *= mparam.curvaturesafety;
    /*
    double hret;

    if (mparam.maxh * kappa < 1)
      hret = mparam.maxh;
    else
      hret = 1 / (kappa + VSMALL);

    if (mparam.maxh < hret)
      hret = mparam.maxh;

    return hret;
    */
    // return min(mparam.maxh, 1/kappa);
    return (mparam.maxh*kappa < 1) ? mparam.maxh : 1/kappa;
  }



  void RestrictHTriangle (gp_Pnt2d & par0, gp_Pnt2d & par1, gp_Pnt2d & par2,
                          BRepLProp_SLProps * prop, BRepLProp_SLProps * prop2, Mesh & mesh, int depth, double h,
                          const MeshingParameters & mparam)
  {
    int ls = -1;

    gp_Pnt pnt0,pnt1,pnt2;

    prop->SetParameters (par0.X(), par0.Y());
    pnt0 = prop->Value();

    prop->SetParameters (par1.X(), par1.Y());
    pnt1 = prop->Value();

    prop->SetParameters (par2.X(), par2.Y());
    pnt2 = prop->Value();

    double aux;
    double maxside = pnt0.Distance(pnt1);
    ls = 2;
    aux = pnt1.Distance(pnt2);
    if(aux > maxside)
      {
        maxside = aux;
        ls = 0;
      }
    aux = pnt2.Distance(pnt0);
    if(aux > maxside)
      {
        maxside = aux;
        ls = 1;
      }



    gp_Pnt2d parmid;

    parmid.SetX( (par0.X()+par1.X()+par2.X()) / 3 );
    parmid.SetY( (par0.Y()+par1.Y()+par2.Y()) / 3 );

    if (depth%3 == 0)
      {
        double curvature = 0;

        prop2->SetParameters (parmid.X(), parmid.Y());
        if (!prop2->IsCurvatureDefined())
          {
            (*testout) << "curvature not defined!" << endl;
            return;
          }
        curvature = max(fabs(prop2->MinCurvature()),
                        fabs(prop2->MaxCurvature()));

        prop2->SetParameters (par0.X(), par0.Y());
        if (!prop2->IsCurvatureDefined())
          {
            (*testout) << "curvature not defined!" << endl;
            return;
          }
        curvature = max(curvature,max(fabs(prop2->MinCurvature()),
                                      fabs(prop2->MaxCurvature())));

        prop2->SetParameters (par1.X(), par1.Y());
        if (!prop2->IsCurvatureDefined())
          {
            (*testout) << "curvature not defined!" << endl;
            return;
          }
        curvature = max(curvature,max(fabs(prop2->MinCurvature()),
                                      fabs(prop2->MaxCurvature())));

        prop2->SetParameters (par2.X(), par2.Y());
        if (!prop2->IsCurvatureDefined())
          {
            (*testout) << "curvature not defined!" << endl;
            return;
          }
        curvature = max(curvature,max(fabs(prop2->MinCurvature()),
                                      fabs(prop2->MaxCurvature())));

        //(*testout) << "curvature " << curvature << endl;

        if (curvature < 1e-3)
          {
            //(*testout) << "curvature too small (" << curvature << ")!" << endl;
            return;
            // return war bis 10.2.05 auskommentiert
          }



        h = ComputeH (curvature+1e-10, mparam);

        if(h < 1e-4*maxside)
          return;


        // if (h > 30) return;
      }

    if (h < maxside && depth < 10)
      {
        //cout << "\r h " << h << flush;
        gp_Pnt2d pm;

        //cout << "h " << h << " maxside " << maxside << " depth " << depth << endl;
        //cout << "par0 " << par0.X() << " " << par0.Y()
        //<< " par1 " << par1.X() << " " << par1.Y()
        //   << " par2 " << par2.X() << " " << par2.Y()<< endl;

        if(ls == 0)
          {
            pm.SetX(0.5*(par1.X()+par2.X())); pm.SetY(0.5*(par1.Y()+par2.Y()));
            RestrictHTriangle(pm, par2, par0, prop, prop2, mesh, depth+1, h, mparam);
            RestrictHTriangle(pm, par0, par1, prop, prop2, mesh, depth+1, h, mparam);
          }
        else if(ls == 1)
          {
            pm.SetX(0.5*(par0.X()+par2.X())); pm.SetY(0.5*(par0.Y()+par2.Y()));
            RestrictHTriangle(pm, par1, par2, prop, prop2, mesh, depth+1, h, mparam);
            RestrictHTriangle(pm, par0, par1, prop, prop2, mesh, depth+1, h, mparam);
          }
        else if(ls == 2)
          {
            pm.SetX(0.5*(par0.X()+par1.X())); pm.SetY(0.5*(par0.Y()+par1.Y()));
            RestrictHTriangle(pm, par1, par2, prop, prop2, mesh, depth+1, h, mparam);
            RestrictHTriangle(pm, par2, par0, prop, prop2, mesh, depth+1, h, mparam);
          }

      }
    else
      {
        gp_Pnt pnt;
        Point3d p3d;

        prop->SetParameters (parmid.X(), parmid.Y());
        pnt = prop->Value();
        p3d = Point3d(pnt.X(), pnt.Y(), pnt.Z());
        mesh.RestrictLocalH (p3d, h);

        p3d = Point3d(pnt0.X(), pnt0.Y(), pnt0.Z());
        mesh.RestrictLocalH (p3d, h);

        p3d = Point3d(pnt1.X(), pnt1.Y(), pnt1.Z());
        mesh.RestrictLocalH (p3d, h);

        p3d = Point3d(pnt2.X(), pnt2.Y(), pnt2.Z());
        mesh.RestrictLocalH (p3d, h);

        //(*testout) << "p = " << p3d << ", h = " << h << ", maxside = " << maxside << endl;

      }
  }

  bool OCCMeshFace (const OCCGeometry & geom, Mesh & mesh, FlatArray<int, PointIndex> glob2loc,
                       const MeshingParameters & mparam, int nr, int projecttype, bool delete_on_failure)
  {
    auto k = nr+1;
    if(1==0 && !geom.fvispar[k-1].IsDrawable())
      {
        (*testout) << "ignoring face " << k << endl;
        cout << "ignoring face " << k << endl;
        return true;
      }

    // if(master_faces[k]!=k)
    //     continue;

    (*testout) << "mesh face " << k << endl;
    multithread.percent = 100 * k / (mesh.GetNFD() + VSMALL);
    geom.facemeshstatus[k-1] = -1;

    FaceDescriptor & fd = mesh.GetFaceDescriptor(k);
    auto face = TopoDS::Face(geom.fmap(k));
    auto fshape = face.TShape();

    int oldnf = mesh.GetNSE();

    Box<3> bb = geom.GetBoundingBox();

    // int projecttype = PLANESPACE;
    // int projecttype = PARAMETERSPACE;
    
    static Timer tinit("init");
    tinit.Start();
    Meshing2OCCSurfaces meshing(geom, face, bb, projecttype, mparam);
    tinit.Stop();


    static Timer tprint("print");
    tprint.Start();
    if (meshing.GetProjectionType() == PLANESPACE)
      PrintMessage (2, "Face ", k, " / ", geom.GetNFaces(), " (plane space projection)");
    else
      PrintMessage (2, "Face ", k, " / ", geom.GetNFaces(), " (parameter space projection)");
    tprint.Stop();

    //      Meshing2OCCSurfaces meshing(f2, bb);
    // meshing.SetStartTime (starttime);
    //(*testout) << "Face " << k << endl << endl;
    

    auto segments = geom.GetFace(k-1).GetBoundary(mesh);

    if (meshing.GetProjectionType() == PLANESPACE)
      {
        static Timer t("MeshSurface: Find edges and points - Physical"); RegionTimer r(t);
        int cntp = 0;
        glob2loc = 0;

        for (Segment & seg : segments)
          // if (seg.si == k)
            for (int j = 0; j < 2; j++)
              {
                PointIndex pi = seg[j];
                if (glob2loc[pi] == 0)
                  {
                    meshing.AddPoint (mesh.Point(pi), pi);
                    cntp++;
                    glob2loc[pi] = cntp;
                  }
              }

        /*
        for (int i = 1; i <= mesh.GetNSeg(); i++)
          {
            Segment & seg = mesh.LineSegment(i);
        */
        // for (Segment & seg : mesh.LineSegments())                
        for (Segment & seg : segments)
          //if (seg.si == k)
            {
              PointGeomInfo gi0, gi1;
              gi0.trignum = gi1.trignum = k;
              gi0.u = seg.epgeominfo[0].u;
              gi0.v = seg.epgeominfo[0].v;
              gi1.u = seg.epgeominfo[1].u;
              gi1.v = seg.epgeominfo[1].v;
              
              //if(orientation & 1)
              meshing.AddBoundaryElement (glob2loc[seg[0]], glob2loc[seg[1]], gi0, gi1);

            }
      }
    else
      {
        static Timer t("MeshSurface: Find edges and points - Parameter"); RegionTimer r(t);
        
        Array<PointGeomInfo> gis(2*segments.Size());
        gis.SetSize (0);

        Box<2> uv_box(Box<2>::EMPTY_BOX);
        for(auto & seg : segments)
            for(auto i : Range(2))
                uv_box.Add( {seg.epgeominfo[i].u, seg.epgeominfo[i].v } );

        BoxTree<2> uv_tree(uv_box);
        Array<int> found_points;

        for(auto & seg : segments)
        {
            PointGeomInfo gi[2];
            gi[0].trignum = gi[1].trignum = k;
            gi[0].u = seg.epgeominfo[0].u;
            gi[0].v = seg.epgeominfo[0].v;
            gi[1].u = seg.epgeominfo[1].u;
            gi[1].v = seg.epgeominfo[1].v;

            int locpnum[2] = {0, 0};

            for (int j = 0; j < 2; j++)
            {
                Point<2> uv = {gi[j].u, gi[j].v};
                uv_tree.GetIntersecting(uv, uv, found_points);

                if(found_points.Size())
                    locpnum[j] = found_points[0];
                else
                {
                    PointIndex pi = seg[j];
                    meshing.AddPoint (mesh.Point(pi), pi);

                    gis.Append (gi[j]);
                    locpnum[j] = gis.Size();
                    uv_tree.Insert(uv, locpnum[j]);
                }
            }

            meshing.AddBoundaryElement (locpnum[0], locpnum[1], gi[0], gi[1]);
        }
      }


    // Philippose - 15/01/2009
    double maxh = min2(geom.face_maxh[k-1], OCCGeometry::global_shape_properties[TopoDS::Face(geom.fmap(k)).TShape()].maxh);
    //double maxh = mparam.maxh;
    //      int noldpoints = mesh->GetNP();
    int noldsurfel = mesh.GetNSE();

    static Timer tsurfprop("surfprop");
    tsurfprop.Start();
    GProp_GProps sprops;
    BRepGProp::SurfaceProperties(TopoDS::Face(geom.fmap(k)),sprops);
    tsurfprop.Stop();
    meshing.SetMaxArea(2.*sprops.Mass());

    MESHING2_RESULT res;

    // TODO: check overlap not correctly working here
    MeshingParameters mparam_without_overlap = mparam;
    mparam_without_overlap.checkoverlap = false;
    
    try {
      static Timer t("GenerateMesh"); RegionTimer reg(t);
      res = meshing.GenerateMesh (mesh, mparam_without_overlap, maxh, k);
    }

    catch (SingularMatrixException)
      {
        // (*myerr) << "Singular Matrix" << endl;
        res = MESHING2_GIVEUP;
      }

    catch (UVBoundsException)
      {
        // (*myerr) << "UV bounds exceeded" << endl;
        res = MESHING2_GIVEUP;
      }

    static Timer t1("rest of loop"); RegionTimer reg1(t1);
      
    bool meshing_failed = res != MESHING2_OK;
    if(meshing_failed && delete_on_failure)
    {
        for (SurfaceElementIndex sei = noldsurfel; sei < mesh.GetNSE(); sei++)
            mesh.Delete(sei);

        mesh.Compress();
    }

    for (SurfaceElementIndex sei = oldnf; sei < mesh.GetNSE(); sei++)
      mesh[sei].SetIndex (k);

    auto n_illegal_trigs = mesh.FindIllegalTrigs();
    PrintMessage (3, n_illegal_trigs, " illegal triangles");
    return meshing_failed;
  }


  void OCCSetLocalMeshSize(const OCCGeometry & geom, Mesh & mesh,
                           const MeshingParameters & mparam, const OCCParameters& occparam)
  {
    static Timer t1("OCCSetLocalMeshSize");
    RegionTimer regt(t1);
    mesh.SetGlobalH (mparam.maxh);
    mesh.SetMinimalH (mparam.minh);

    NgArray<double> maxhdom;
    maxhdom.SetSize (geom.NrSolids());
    maxhdom = mparam.maxh;

    int dom = 0;
    for (TopExp_Explorer e(geom.GetShape(), TopAbs_SOLID); e.More(); e.Next(), dom++)
      maxhdom[dom] = min2(maxhdom[dom], OCCGeometry::global_shape_properties[e.Current().TShape()].maxh);

    mesh.SetMaxHDomain (maxhdom);

    Box<3> bb = geom.GetBoundingBox();
    bb.Increase (bb.Diam()/10);

    mesh.SetLocalH (bb.PMin(), bb.PMax(), 0.5);

    if (mparam.uselocalh)
      {
        const char * savetask = multithread.task;
        multithread.percent = 0;

        mesh.SetLocalH (bb.PMin(), bb.PMax(), mparam.grading);

        int nedges = geom.emap.Extent();

        double mincurvelength = IGNORECURVELENGTH;
        double maxedgelen = 0;
        double minedgelen = 1e99;

        if(occparam.resthminedgelenenable) 
          {
            mincurvelength = occparam.resthminedgelen;
            if(mincurvelength < IGNORECURVELENGTH) mincurvelength = IGNORECURVELENGTH;
          }

        multithread.task = "Setting local mesh size (elements per edge)";

        // Philippose - 23/01/2009
        // Find all the parent faces of a given edge
        // and limit the mesh size of the edge based on the
        // mesh size limit of the face
        TopTools_IndexedDataMapOfShapeListOfShape edge_face_map;
        edge_face_map.Clear();
        TopExp::MapShapesAndAncestors(geom.shape, TopAbs_EDGE, TopAbs_FACE, edge_face_map);

        // setting elements per edge
        for (int i = 1; i <= nedges && !multithread.terminate; i++)
          {
            TopoDS_Edge e = TopoDS::Edge (geom.emap(i));
            multithread.percent = 100 * (i-1)/double(nedges);
            if (BRep_Tool::Degenerated(e)) continue;

            GProp_GProps system;
            BRepGProp::LinearProperties(e, system);
            double len = system.Mass();

            if (len < mincurvelength)
              {
                (*testout) << "ignored" << endl;
                continue;
              }

            bool is_identified_edge = false;
            // TODO: change to use hash value
            const auto& gedge = geom.GetEdge(geom.edge_map.at(e.TShape()));
            auto& v0 = gedge.GetStartVertex();
            auto& v1 = gedge.GetEndVertex();
            for(auto & ident : v0.identifications)
            {
                auto other = ident.from == &v0 ? ident.to : ident.from;
                if(other->nr == v1.nr && ident.type == Identifications::CLOSESURFACES)
                {
                    is_identified_edge = true;
                    break;
                }
            }

            if(is_identified_edge)
              continue;

            double localh = len/mparam.segmentsperedge;
            double s0, s1;

            const TopTools_ListOfShape& parent_faces = edge_face_map.FindFromKey(e);

            TopTools_ListIteratorOfListOfShape parent_face_list;

            for(parent_face_list.Initialize(parent_faces); parent_face_list.More(); parent_face_list.Next())
              {
                TopoDS_Face parent_face = TopoDS::Face(parent_face_list.Value());

                int face_index = geom.fmap.FindIndex(parent_face);

                if(face_index >= 1) localh = min(localh,geom.face_maxh[face_index - 1]);
                localh = min2(localh, OCCGeometry::global_shape_properties[parent_face.TShape()].maxh);
              }

            Handle(Geom_Curve) c = BRep_Tool::Curve(e, s0, s1);

            localh = min2(localh, OCCGeometry::global_shape_properties[e.TShape()].maxh);
            maxedgelen = max (maxedgelen, len);
            minedgelen = min (minedgelen, len);
            int maxj = max((int) ceil(len/localh), 2);

            for (int j = 0; j <= maxj; j++)
              {
                gp_Pnt pnt = c->Value (s0+double(j)/maxj*(s1-s0));
                mesh.RestrictLocalH (Point3d(pnt.X(), pnt.Y(), pnt.Z()), localh);
              }
          }

        multithread.task = "Setting local mesh size (edge curvature)";

        // setting edge curvature

        int nsections = 20;

        for (int i = 1; i <= nedges && !multithread.terminate; i++)
          {
            double maxcur = 0;
            multithread.percent = 100 * (i-1)/double(nedges);
            TopoDS_Edge edge = TopoDS::Edge (geom.emap(i));
            if (BRep_Tool::Degenerated(edge)) continue;
            double s0, s1;
            Handle(Geom_Curve) c = BRep_Tool::Curve(edge, s0, s1);
            BRepAdaptor_Curve brepc(edge);
            BRepLProp_CLProps prop(brepc, 2, 1e-5);

            for (int j = 1; j <= nsections; j++)
              {
                double s = s0 + j/(double) nsections * (s1-s0);
                prop.SetParameter (s);
                double curvature = 0;
                if(prop.IsTangentDefined())
                  curvature = prop.Curvature();
                if(curvature> maxcur) maxcur = curvature;

                if (curvature >= 1e99)
                  continue;

                gp_Pnt pnt = c->Value (s);

                mesh.RestrictLocalH (Point3d(pnt.X(), pnt.Y(), pnt.Z()), ComputeH (fabs(curvature), mparam));
              }
          }

        multithread.task = "Setting local mesh size (face curvature)";

        // setting face curvature

        int nfaces = geom.fmap.Extent();

        for (int i = 1; i <= nfaces && !multithread.terminate; i++)
          {
            multithread.percent = 100 * (i-1)/double(nfaces);
            TopoDS_Face face = TopoDS::Face(geom.fmap(i));
            TopLoc_Location loc;
            Handle(Geom_Surface) surf = BRep_Tool::Surface (face);
            Handle(Poly_Triangulation) triangulation = BRep_Tool::Triangulation (face, loc);

            if (triangulation.IsNull())
              {
                BRepTools::Clean (geom.shape);
                BRepMesh_IncrementalMesh (geom.shape, 0.01, true);
                triangulation = BRep_Tool::Triangulation (face, loc);
              }

            BRepAdaptor_Surface sf(face, Standard_True);
            // one prop for evaluating and one for derivatives
            BRepLProp_SLProps prop(sf, 0, 1e-5);
            BRepLProp_SLProps prop2(sf, 2, 1e-5);

            int ntriangles = triangulation -> NbTriangles();
            for (int j = 1; j <= ntriangles; j++)
              {
                gp_Pnt p[3];
                gp_Pnt2d par[3];

                for (int k = 1; k <=3; k++)
                  {
                    // int n = triangulation->Triangles()(j)(k);
                    // p[k-1] = triangulation->Nodes()(n).Transformed(loc);
                    // par[k-1] = triangulation->UVNodes()(n);
                    // fix for OCC7.6.0-dev
                    int n = triangulation->Triangle(j)(k);
                    p[k-1] = triangulation->Node(n).Transformed(loc);
                    par[k-1] = triangulation->UVNode(n);
                  }

                //double maxside = 0;
                //maxside = max (maxside, p[0].Distance(p[1]));
                //maxside = max (maxside, p[0].Distance(p[2]));
                //maxside = max (maxside, p[1].Distance(p[2]));
                //cout << "\rFace " << i << " pos11 ntriangles " << ntriangles << " maxside " << maxside << flush;

                RestrictHTriangle (par[0], par[1], par[2], &prop, &prop2, mesh, 0, 0, mparam);
                //cout << "\rFace " << i << " pos12 ntriangles " << ntriangles << flush;
              }
          }

        // setting close edges

        if (mparam.closeedgefac.has_value())
          {
            multithread.task = "Setting local mesh size (close edges)";

            int sections = 100;

            NgArray<Line> lines(sections*nedges);

            /*
            BoxTree<3> * searchtree =
              new BoxTree<3> (bb.PMin(), bb.PMax());
            */
            BoxTree<3> searchtree(bb.PMin(), bb.PMax());
            
            int nlines = 0;
            Array<int> edgenumber;
            for (int i = 1; i <= nedges && !multithread.terminate; i++)
              {
                TopoDS_Edge edge = TopoDS::Edge (geom.emap(i));
                if (BRep_Tool::Degenerated(edge)) continue;

                double s0, s1;
                Handle(Geom_Curve) c = BRep_Tool::Curve(edge, s0, s1);
                BRepAdaptor_Curve brepc(edge);
                BRepLProp_CLProps prop(brepc, 1, 1e-5);
                prop.SetParameter (s0);

                gp_Vec d0 = prop.D1().Normalized();
                double s_start = s0;
                int count = 0;
                for (int j = 1; j <= sections; j++)
                  {
                    double s = s0 + (s1-s0)*(double)j/(double)sections;
                    prop.SetParameter (s);
                    gp_Vec d1 = prop.D1().Normalized();
                    double cosalpha = fabs(d0*d1);
                    if ((j == sections) || (cosalpha < cos(10.0/180.0*M_PI)))
                      {
                        count++;
                        gp_Pnt p0 = c->Value (s_start);
                        gp_Pnt p1 = c->Value (s);
                        lines[nlines].p0 = Point<3> (p0.X(), p0.Y(), p0.Z());
                        lines[nlines].p1 = Point<3> (p1.X(), p1.Y(), p1.Z());

                        Box3d box;
                        box.SetPoint (Point3d(lines[nlines].p0));
                        box.AddPoint (Point3d(lines[nlines].p1));

                        searchtree.Insert (box.PMin(), box.PMax(), nlines+1);
                        nlines++;
                        edgenumber.Append(i);

                        s_start = s;
                        d0 = d1;
                      }
                  }
              }

            NgArray<int> linenums;
            auto is_identified_edge = [&](int e0, int e1) {
                const auto& edge0 = geom.GetEdge(e0-1);
                const auto& edge1 = geom.GetEdge(e1-1);

                if(edge0.primary == edge1.primary)
                    return true;

                Array<const GeometryVertex *> v0 = { &edge0.GetStartVertex(), &edge0.GetEndVertex() };
                Array<const GeometryVertex *> v1 = { &edge1.GetStartVertex(), &edge1.GetEndVertex() };
                for(auto i : Range(2))
                    for(auto j : Range(2))
                        if(v0[i]->primary == v1[j]->primary)
                            return true;

                return false;
            };

            for (int i = 0; i < nlines; i++)
              {
                multithread.percent = (100*i)/double(nlines);
                Line & line = lines[i];

                Box3d box;
                box.SetPoint (Point3d(line.p0));
                box.AddPoint (Point3d(line.p1));
                double maxhline = max (mesh.GetH(box.PMin()),
                                       mesh.GetH(box.PMax()));
                box.Increase(maxhline);

                double mindist = 1e99;
                linenums.SetSize(0);
                searchtree.GetIntersecting(box.PMin(),box.PMax(),linenums);

                for (int j = 0; j < linenums.Size(); j++)
                  {
                    int num = linenums[j]-1;
                    if (i == num) continue;
                    if( is_identified_edge(edgenumber[i], edgenumber[num]) ) continue;
                    if ((line.p0-lines[num].p0).Length2() < 1e-15) continue;
                    if ((line.p0-lines[num].p1).Length2() < 1e-15) continue;
                    if ((line.p1-lines[num].p0).Length2() < 1e-15) continue;
                    if ((line.p1-lines[num].p1).Length2() < 1e-15) continue;
                    mindist = min (mindist, line.Dist(lines[num]));
                  }

                mindist /= (*mparam.closeedgefac + VSMALL);

                if (mindist < 1e-3 * bb.Diam())
                  {
                    (*testout) << "extremely small local h: " << mindist
                               << " --> setting to " << 1e-3 * bb.Diam() << endl;
                    (*testout) << "somewhere near " << line.p0 << " - " << line.p1 << endl;
                    mindist = 1e-3 * bb.Diam();
                  }

                mesh.RestrictLocalHLine(line.p0, line.p1, mindist);
              }
          }

        for (auto mspnt : mparam.meshsize_points)
          mesh.RestrictLocalH(mspnt.pnt, mspnt.h);

        multithread.task = savetask;

      }

    mesh.LoadLocalMeshSize (mparam.meshsizefilename);
  }
}

#endif
