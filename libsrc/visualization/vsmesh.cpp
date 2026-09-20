#include <mystdlib.h>

#include <myadt.hpp>
#include <meshing.hpp>
#include "../meshing/global.hpp"
// #include <csg.hpp>

#ifdef STLGEOM
#include <stlgeom.hpp>
#endif


// #include <parallel.hpp>

#include <visual.hpp>

namespace netgen
{
  // extern shared_ptr<Mesh> mesh;
  extern NetgenGeometry * ng_geometry;

  VisualSceneMesh vsmesh;

  VisualSceneMesh :: VisualSceneMesh ()
    : VisualScene()
  {
    selface = -1;
    selelement = -1;
    locpi = -2;
    selpoint = PointIndex::INVALID;
    selpoint2 = PointIndex::INVALID;
    seledge = -1;

    minh = 0.0;
    maxh = 0.0;
    user_me_handler = NULL;
    mesh = nullptr;
  }

  VisualSceneMesh :: ~VisualSceneMesh ()
  {
    ;
  }


  void VisualSceneMesh :: DrawScene ()
  {
    try
      {
        shared_ptr<Mesh> mesh = GetMesh();

    if (!mesh)
      {
        VisualScene::DrawScene();
        return;
      }

    static Timer timer("VSMesh::DrawScene");

    RegionTimer reg (timer);

    BuildScene();

    glEnable(GL_DEPTH_TEST);
    glClearColor(backcolor, backcolor, backcolor, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glEnable (GL_COLOR_MATERIAL);
    glColor3f (1.0f, 1.0f, 1.0f);
    glLineWidth (1.0f);

    SetLight();

    glPushMatrix();
    glMultMatrixd (transformationmat);

    GLdouble projmat[16];                 // brauchen wir das ?
    glGetDoublev (GL_PROJECTION_MATRIX, projmat);


    glEnable (GL_BLEND);
    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


    glInitNames ();
    glPushName (0);

    //    glEnable (GL_LINE_SMOOTH);
    //         glEnable (GL_BLEND);
    //         glEnable (GL_POLYGON_SMOOTH);
    //         glDisable (GL_DEPTH_TEST);
    //         glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    //    glHint (GL_LINE_SMOOTH_HINT, GL_DONT_CARE);

    glDisable (GL_COLOR_MATERIAL);

    GLfloat matcol0[] = { 0, 0, 0, 1 };
    GLfloat matcol1[] = { 1, 1, 1, 1 };
    GLfloat matcolf[] = { 0, 1, 0, 1 };
    GLfloat matcolb[] = { 0.5, 0, 0, 1 };
    // GLfloat matcolblue[] = { 0, 0, 1, 1 };

    glMatrixMode (GL_MODELVIEW);

    glMaterialfv(GL_FRONT, GL_EMISSION, matcol0);
    glMaterialfv(GL_BACK, GL_EMISSION, matcol0);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, matcol1);
    glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, matcolf);
    glMaterialfv(GL_BACK, GL_AMBIENT_AND_DIFFUSE, matcolb);

    glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);

    // glPolygonOffset (1,10);
    glPolygonOffset (2,2);
    glEnable (GL_POLYGON_OFFSET_FILL);

    SetClippingPlane ();

    if (vispar.drawfilledtrigs)
      {
        BuildFilledList (false);


#ifdef PARALLELGL
        if (ntasks > 1 && vispar.drawtetsdomain > 0 && vispar.drawtetsdomain < ntasks)
          glCallList (par_filledlists[vispar.drawtetsdomain]);
        else
#endif
          glCallList (filledlist);
      }

    if (vispar.drawbadels)
      glCallList (badellist);

    TBitArray<PointIndex> shownode(mesh->GetNP());
    if (vispar.clipping.enable)
      {
        shownode.Clear();
        for (PointIndex pi : mesh->Points().Range())
          {
            Point<3> p = (*mesh)[pi];

            double val =
              p[0] * clipplane[0] +
              p[1] * clipplane[1] +
              p[2] * clipplane[2] +
              clipplane[3];

            if (val > 0) shownode.SetBit (pi);
          }
      }
    else
      shownode.Set();
    if (vispar.drawprisms)
      {
        BuildPrismList (shownode);
        glCallList (prismlist);
      }

    if (vispar.drawpyramids)
      {
        BuildPyramidList (shownode);
        glCallList (pyramidlist);
      }

    if (vispar.drawhexes)
      {
        BuildHexList (shownode);
        glCallList (hexlist);
      }

    if (vispar.drawtets)
      {
        BuildTetList (shownode);
        glCallList (tetlist);
      }

    if (vispar.drawdomainsurf)
      {
        BuildDomainSurfList();
        glCallList (domainsurflist);
      }

    glDisable (GL_POLYGON_OFFSET_FILL);

    // draw lines

    glMatrixMode (GL_MODELVIEW);

    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, matcol0);
    glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, matcol0);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, matcol0);

    glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
    glLineWidth (1.0f);
    glColor3f (0.0f, 0.0f, 0.0f);
    glDisable (GL_LINE_SMOOTH);

    glDisable(GL_BLEND);


    if (vispar.drawoutline)
      {
        glPolygonOffset (1, 1);
        glEnable (GL_POLYGON_OFFSET_LINE);

        BuildLineList ();

#ifdef PARALLELGL
        if (ntasks > 1 && vispar.drawtetsdomain > 0 && vispar.drawtetsdomain < ntasks)
          glCallList (par_linelists[vispar.drawtetsdomain]);
        else
#endif
          glCallList (linelist);


        glDisable (GL_POLYGON_OFFSET_LINE);
      }

    if (vispar.drawidentified)
      {
        glPolygonOffset (1, -1);
        glEnable (GL_POLYGON_OFFSET_LINE);
        glCallList (identifiedlist);
        glDisable (GL_POLYGON_OFFSET_LINE);
      }

    if (vispar.drawpointnumbers ||
        vispar.drawedgenumbers ||
        vispar.drawfacenumbers ||
        vispar.drawsegmentnumbers ||
        vispar.drawsurfaceelementnumbers ||
        vispar.drawelementnumbers)
      glCallList (pointnumberlist);


    glPopName();

    if (vispar.drawedges)
      {
        BuildEdgeList();
        glCallList (edgelist);
      }

    DrawMarker();

    glDisable(GL_CLIP_PLANE0);

    glPopMatrix();

    if (vispar.colormeshsize)
      DrawColorBar (minh, maxh, 1);

    DrawCoordinateCross ();
    DrawNetgenLogo ();


    lock = {};
    
    glFinish();

    
      }
    catch (const bad_weak_ptr & e)
      {
        // cout << "don't have a mesh to visualize" << endl;
        VisualScene::DrawScene();      
      }

  }

  void VisualSceneMesh :: SelectCenter (int zoomall)
  {
    shared_ptr<Mesh> mesh = GetMesh();
    Point<3> pmin, pmax;
    mesh->GetBox (pmin, pmax, -1);

    // works in NGSolve, mesh view
    if (mesh->GetDimension() == 2)
      mesh->GetBox (pmin, pmax);
    else // otherwise strange zooms during mesh generation
      mesh->GetBox (pmin, pmax, SURFACEPOINT);

    if (vispar.use_center_coords && zoomall==2)
    {
      center(0) = vispar.centerx;
      center(1) = vispar.centery;
      center(2) = vispar.centerz;
    }
    else if (selpoint-IndexBASE<PointIndex>() >= 1 && zoomall==2)
      center = mesh->Point (selpoint);
    else if (marker && zoomall==2)
      center = *marker;
    else if (vispar.centerpoint.IsValid() && zoomall==2)
      center = mesh->Point (vispar.centerpoint);
    else
      center = Center (pmin, pmax);

    double oldrad = rad;
    rad = 0.5 * Dist (pmin, pmax);
    if(rad == 0) rad = 1e-6;

    if (rad > 1.2 * oldrad ||
        mesh->GetMajorTimeStamp() > vstimestamp ||
        zoomall)
    {
      CalcTransformationMatrices();
    }

    glEnable (GL_NORMALIZE);

  }

  void VisualSceneMesh :: BuildScene (int zoomall)
  {
        auto mesh = GetGlobalMesh();
        
        if (!mesh)
      {
        PrintMessage (3, "vsmesh::buildscene: don't have a mesh to visualize");
        VisualScene::BuildScene (zoomall);
        return;
      }

    if (!lock)
      {
        lock = std::unique_lock<std::mutex> (mesh->Mutex());
      }

    static Timer timer("VSMesh::BuildScene");
    RegionTimer reg (timer);



    Array<ElementFace> faces;

    int meshtimestamp = mesh->GetTimeStamp();
    if (meshtimestamp > vstimestamp || zoomall)
        SelectCenter(zoomall);

    if (pointnumberlist)
      {
        glDeleteLists (pointnumberlist, 1);
        pointnumberlist = 0;
      }

    if (badellist)
      {
        glDeleteLists (badellist, 1);
        badellist = 0;
      }
    /*
      if (prismlist)
      {
      glDeleteLists (prismlist, 1);
      prismlist = 0;
      }

      if (pyramidlist)
      {
      glDeleteLists (pyramidlist, 1);
      pyramidlist = 0;
      }

      if (hexlist)
      {
      glDeleteLists (hexlist, 1);
      hexlist = 0;
      }
    */
    if (identifiedlist)
      {
        glDeleteLists (identifiedlist, 1);
        identifiedlist = 0;
      }

    pointnumberlist = glGenLists (1);
    glNewList (pointnumberlist, GL_COMPILE);

    if (vispar.drawpointnumbers ||
        vispar.drawedgenumbers ||
        vispar.drawfacenumbers ||
        vispar.drawsegmentnumbers ||
        vispar.drawsurfaceelementnumbers ||
        vispar.drawelementnumbers)
      {
        //      glEnable (GL_COLOR_MATERIAL);
        GLfloat textcol[3] = { float(1-backcolor),
                               float(1-backcolor),
                               float(1-backcolor) };
        glColor3fv (textcol);
        glNormal3d (0, 0, 1);
        glPushAttrib (GL_LIST_BIT);
        // glListBase (fontbase);

        char buf[30];

        if (vispar.drawpointnumbers)
          for (PointIndex pi : mesh->Points().Range())
            {
              const Point<3> & p = mesh->Point(pi);
              glRasterPos3d (p(0), p(1), p(2));

              snprintf (buf, size(buf),  "%d", pi.Nr1());

              // glCallLists (strlen (buf), GL_UNSIGNED_BYTE, buf);
              MyOpenGLText (buf);
            }

        if (vispar.drawedgenumbers)
          {
            /*
              for (SegmentIndex i = 0; i < mesh->GetNSeg(); i++)
              {
              const Segment & seg = (*mesh)[i];

              const Point<3> & p1 = mesh->Point(seg[0]);
              const Point<3> & p2 = mesh->Point(seg[1]);
              const Point<3> p = Center (p1, p2);
              glRasterPos3d (p.X(), p.Y(), p.Z());

              snprintf (buf, size(buf),  "%d", seg.edgenr);
              glCallLists (strlen (buf), GL_UNSIGNED_BYTE, buf);
              }
            */

            const MeshTopology & top = mesh->GetTopology();
            for (int i = 1; i <= top.GetNEdges(); i++)
              {
                // int v1, v2;
                // top.GetEdgeVertices (i, v1, v2);
                auto [v1,v2] = top.GetEdgeVertices(i-1);
                const Point<3> & p1 = mesh->Point(v1);
                const Point<3> & p2 = mesh->Point(v2);
                const Point<3> p = Center (p1, p2);
                glRasterPos3d (p(0), p(1), p(2));

                snprintf (buf, size(buf),  "%d", i);
                // glCallLists (strlen (buf), GL_UNSIGNED_BYTE, buf);
                MyOpenGLText (buf);

              }

          }

          if (vispar.drawsegmentnumbers)
            {
              for (auto si : Range(mesh->LineSegments())) {
                const auto& seg = (*mesh)[si];
                Point<3> c = Center((*mesh)[seg[0]], (*mesh)[seg[1]]);
                glRasterPos3d (c[0], c[1], c[2]);
                snprintf (buf, size(buf),  "%d", si.Nr0());
                MyOpenGLText (buf);
              }
            }

          if (vispar.drawfacenumbers)
          {
            const MeshTopology & top = mesh->GetTopology();
            Array<int> v;
            for (int i = 1; i <= top.GetNFaces(); i++)
              {
                top.GetFaceVertices (i, v);
                auto P = [&] (int j) -> const Point<3> & { return mesh->Point(PointIndex::FromNr1(v[j])); };
                Point<3> p;
                if (v.Size() == 3)
                  {
                    p = Center (P(0), P(1), P(2));
                  }
                else
                  {
                    p = Center (Center (P(0), P(1)), Center (P(2), P(3)));
                  }

                glRasterPos3d (p(0), p(1), p(2));
                snprintf (buf, size(buf),  "%d", i);
                // glCallLists (strlen (buf), GL_UNSIGNED_BYTE, buf);
                MyOpenGLText (buf);
              }
          }

          if (vispar.drawsurfaceelementnumbers)
            {
              for (auto sei : Range(mesh->SurfaceElements()))
                {
                  const auto & sel = (*mesh)[sei];
                  Point<3> c;
                  if(sel.GetNV() == 3)
                    c = Center((*mesh)[sel[0]],
                               (*mesh)[sel[1]],
                               (*mesh)[sel[2]]);
                  else
                    c = Center((*mesh)[sel[0]],
                               (*mesh)[sel[1]],
                               (*mesh)[sel[2]],
                               (*mesh)[sel[3]]);
                  glRasterPos3d (c[0], c[1], c[2]);
                  snprintf (buf, size(buf),  "%d", sei.Nr0());
                  MyOpenGLText (buf);
                }
            }

        if (vispar.drawelementnumbers)
          {
            Array<int> v;
            // for (int i = 1; i <= mesh->GetNE(); i++)
            for (ElementIndex ei : Range(mesh->VolumeElements()))
              {
                netgen::Point<3> p;
                auto el = mesh->VolumeElement (ei);
                auto P = [&] (int j) -> const netgen::Point<3> & { return mesh->Point(el.PNum(j)); };

                switch (el.GetNV())
                  {
                  case 4:
                    p = Center (P(1), P(2), P(3), P(4));
                    break;
                  case 5:
                    {
                      auto c = Center (Center(P(1),P(3)), Center(P(2),P(4)));
                      const auto & p5 = P(5);
                      for (int j = 0; j < 3; j++)
                        p(j) = 0.3 * p5(j) + 0.7 * c(j);
                      break;
                    }
                  case 6:
                    p = Center (Center(P(1),P(2),P(3)), Center(P(4),P(5),P(6)));
                    break;
                  case 8:
                    p = Center (Center(Center(P(1),P(3)), Center(P(2),P(4))),
                                Center(Center(P(5),P(7)), Center(P(6),P(8))));
                    break;
                  }

                glRasterPos3d (p(0), p(1), p(2));
                snprintf (buf, size(buf),  "%d", ei-IndexBASE(ei));
                // glCallLists (strlen (buf), GL_UNSIGNED_BYTE, buf);
                MyOpenGLText (buf);

              }
          }


        glPopAttrib ();
        //              glDisable (GL_COLOR_MATERIAL);
      }
    glEndList ();






    badellist = glGenLists (1);
    glNewList (badellist, GL_COMPILE);

    if (vispar.drawbadels)
      {
        //  SetClippingPlane ();

        static float badelcol[] = { 1.0f, 0.0f, 1.0f, 1.0f };
        glLineWidth (1.0f);

        //for (int i = 1; i <= mesh->GetNE(); i++)
        for (ElementIndex ei : Range(mesh->VolumeElements()))
          {
            if (mesh->VolumeElement(ei).Flags().badel ||
                mesh->VolumeElement(ei).Flags().illegal ||
                (ei-IndexBASE(ei) == vispar.drawelement))
              {
                // copy to be thread-safe
                Element el (mesh->VolumeElement (ei));
                el.GetSurfaceTriangles (faces);

                glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, badelcol);


                //        if ( (el.GetNP() == 4) || (el.GetNP() == 10))
                if (el[0].IsValid())
                  {
                    glBegin (GL_TRIANGLES);

                    for (int j = 0; j < faces.Size(); j++)
                      {
                        ElementFace & face = faces[j];
                        const Point<3> & lp1 = mesh->Point (el.PNum(face[0]));
                        const Point<3> & lp2 = mesh->Point (el.PNum(face[1]));
                        const Point<3> & lp3 = mesh->Point (el.PNum(face[2]));
                        Vec<3> n = Cross (Vec<3> (lp1, lp2), Vec<3> (lp1, lp3));
                        n /= (n.Length()+1e-12);
                        glNormal3d (n(0), n(1), n(2));
                        glVertex3d (lp1(0), lp1(1), lp1(2));
                        glVertex3d (lp2(0), lp2(1), lp2(2));
                        glVertex3d (lp3(0), lp3(1), lp3(2));
                      }

                    glEnd();
                  }
              }
          }



        for (auto el2 : mesh->VolumeElements())
          {
            if (el2.Flags().badel)
              {
                // copy to be thread-safe
                Element el (el2);
                if ( (el.GetNP() == 4) || (el.GetNP() == 10))
                  {
                    glBegin (GL_LINES);
                    glVertex3d (0,0,0);
                    const Point<3> & p = mesh->Point(el[0]);
                    glVertex3d (p(0), p(1), p(2));
                    glEnd();
                  }
              }
          }


        for (ElementIndex ei : Range(mesh->VolumeElements()))
          {
            Element el (mesh->VolumeElement (ei));
            int hascp = 0;
            for (int j = 0; j < el.GetNP(); j++)
              if (el[j] == vispar.centerpoint)
                hascp = 1;

            if (hascp)
              {
                (*testout) << "draw el " << ei << " : ";
                for (int j = 0; j < el.GetNP(); j++)
                  (*testout) << el[j] << " ";
                (*testout) << endl;

                if (el.GetNP() == 4)
                  {
                    int et[6][2] =
                      { { 1, 2 },
                        { 1, 3 },
                        { 1, 4 },
                        { 2, 3 },
                        { 2, 4 },
                        { 3, 4 } } ;

                    for (int j = 0; j < 6; j++)
                      {
                        glBegin (GL_LINES);
                        const Point<3> & p1 = mesh->Point (el.PNum(et[j][0]));
                        const Point<3> & p2 = mesh->Point (el.PNum(et[j][1]));
                        glVertex3d (p1(0), p1(1), p1(2));
                        glVertex3d (p2(0), p2(1), p2(2));
                        glEnd ();
                      }
                  }


                if (el.GetNP() == 10)
                  {
                    int et[12][2] =
                      { { 1, 5 },
                        { 2, 5 },
                        { 1, 6 },
                        { 3, 6 },
                        { 1, 7 },
                        { 4, 7 },
                        { 2, 8 },
                        { 3, 8 },
                        { 2, 9 },
                        { 4, 9 },
                        { 3, 10 },
                        { 4, 10 } };

                    for (int j = 0; j < 12; j++)
                      {
                        glBegin (GL_LINES);
                        const Point<3> & p1 = mesh->Point (el.PNum(et[j][0]));
                        const Point<3> & p2 = mesh->Point (el.PNum(et[j][1]));
                        glVertex3d (p1(0), p1(1), p1(2));
                        glVertex3d (p2(0), p2(1), p2(2));
                        glEnd ();
                      }
                  }
              }
          }

        for (auto sel : mesh->SurfaceElements())
          {
            Element2d el (sel); // copy to be thread-safe
            if (!el.BadElement())
              continue;

            if (el.IsDeleted()) continue;
            
            bool drawel = true;
            for (int j = 0; j < el.GetNP(); j++)
              if (!el[j].IsValid())
                drawel = false;

            if (!drawel)
              continue;

            // cout << int (el.GetType()) << " " << flush;
            switch (el.GetType())
              {
              case TRIG:
                {
                  glBegin (GL_TRIANGLES);

                  Point<3> lp1 = mesh->Point (el[0]);
                  Point<3> lp2 = mesh->Point (el[1]);
                  Point<3> lp3 = mesh->Point (el[2]);
                  Vec<3> n = Cross (Vec<3> (lp1, lp2), Vec<3> (lp1, lp3));
                  n /= (n.Length() + 1e-12);
                  glNormal3dv (&n(0));
                  glVertex3dv (&lp1(0));
                  glVertex3dv (&lp2(0));
                  glVertex3dv (&lp3(0));
                  glEnd();
                  break;
                }
              case QUAD:
                {
                  glBegin (GL_QUADS);

                  const Point<3> & lp1 = mesh->Point (el[0]);
                  const Point<3> & lp2 = mesh->Point (el[1]);
                  const Point<3> & lp3 = mesh->Point (el[3]);
                  const Point<3> & lp4 = mesh->Point (el[2]);
                  Vec<3> n = Cross (Vec<3> (lp1, lp2),
                                   Vec<3> (lp1, Center (lp3, lp4)));
                  n /= (n.Length() + 1e-12);
                  glNormal3d (n(0), n(1), n(2));
                  glVertex3d (lp1(0), lp1(1), lp1(2));
                  glVertex3d (lp2(0), lp2(1), lp2(2));
                  glVertex3d (lp4(0), lp4(1), lp4(2));
                  glVertex3d (lp3(0), lp3(1), lp3(2));
                  glEnd();
                  break;
                }
              case TRIG6:
                {
                  int lines[6][2] = {
                    { 1, 6 }, { 2, 6 },
                    { 1, 5 }, { 3, 5 },
                    { 2, 4 }, { 3, 4 } };

                  glBegin (GL_LINES);
                  for (int j = 0; j < 6; j++)
                    {
                      glVertex3dv ( mesh->Point (el.PNum(lines[j][0])) );
                      glVertex3dv ( mesh->Point (el.PNum(lines[j][0])) );
                    }
                  glEnd();
                  break;
                }

              case QUAD6:
                {
                  int lines[6][2] = {
                    { 1, 5 }, { 2, 5 },
                    { 3, 6 }, { 4, 6 },
                    { 1, 4 }, { 2, 3 } };

                  glBegin (GL_LINES);

                  for (int j = 0; j < 6; j++)
                    {
                      const Point<3> & lp1 = mesh->Point (el.PNum(lines[j][0]));
                      const Point<3> & lp2 = mesh->Point (el.PNum(lines[j][1]));

                      glVertex3d (lp1(0), lp1(1), lp1(2));
                      glVertex3d (lp2(0), lp2(1), lp2(2));
                    }
                  glEnd ();
                  break;
                }
              default:
                PrintSysError ("Cannot draw surface element of type ",
                               int(el.GetType()));
              }
          }
        glLoadName (0);

      }
    glEndList ();


    if (1)
      {

        identifiedlist = glGenLists (1);
        glNewList (identifiedlist, GL_COMPILE);

        GLfloat identifiedcol[] = { 1, 0, 1, 1 };

        glLineWidth (3);
        glEnable (GL_COLOR_MATERIAL);
        glDisable (GL_LIGHTING);

        if (mesh -> HasIdentifications() )
          {
              {
                auto & idpts =
                  mesh->GetIdentifications().GetIdentifiedPoints();
                for (auto [hash, val] : idpts)
                  {
                    auto [hash_pts, hash_nr] = hash;
                    auto [pi1, pi2] = hash_pts;
                      // val = pts[2];   
                      Point<3> p1 = mesh->Point(pi1);
                      Point<3> p2 = mesh->Point(pi2);
                      Point<3> c = Center(p1, p2);
                      if (vispar.shrink < 1)
                        {
                          p1 = c + vispar.shrink * (p1 - c);
                          p2 = c + vispar.shrink * (p2 - c);
                        }

                      glColor3fv (identifiedcol);
                      glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,
                                    identifiedcol);

                      glBegin (GL_LINES);
                      glVertex3dv(p1);
                      glVertex3dv(p2);
                      glEnd();
                    }
              }
          }

        glDisable (GL_COLOR_MATERIAL);
        glEnable (GL_LIGHTING);
        glEndList ();
      }

    lock = {};

    vstimestamp = meshtimestamp;

  }



  void VisualSceneMesh :: BuildColorTexture ()
  {
    shared_ptr<Mesh> mesh = GetMesh();

    if(colors.texture == -1)
      glGenTextures(1, &colors.texture);

    // build color texture
    glBindTexture(GL_TEXTURE_2D, colors.texture);
    Array<float> data;
    for(auto fdi : mesh->FaceDescriptors().Range())
    {
      auto c = mesh->GetFaceDescriptor(fdi).SurfColour();
      ArrayMem<float, 4> cf{float(c[0]), float(c[1]), float(c[2]), float(c[3])};
      if(fdi.Nr1()==selface)
        cf = {1.0f, 0.0f, 0.0f, 1.0f};
      data.Append(cf);
    }
    int n = data.Size()/4;
    colors.width = max2(1,min2(n, 1024));
    colors.height = (n+colors.width-1)/colors.width;
    for([[maybe_unused]] auto i: Range(n, colors.width*colors.height))
      data.Append({0.0f, 0.0f, 0.0f, 0.0f});
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, colors.width, colors.height, 0, GL_RGBA, GL_FLOAT, data.Data());
    glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  }

  void RenderSurfaceElements (shared_ptr<Mesh> mesh,
      int subdivisions,
      std::function<bool(int)> face_init,
      std::function<bool(SurfaceElementIndex)> sel_init
    )
  {
    CurvedElements & curv = mesh->GetCurvedElements();

    int hoplotn = 1 << subdivisions;
    
    Array<SurfaceElementIndex> seia;

    for (int faceindex = 1; faceindex <= mesh->GetNFD(); faceindex++)
      {
        if(!face_init(faceindex-1))
          continue;

        mesh->GetSurfaceElementsOfFace (faceindex, seia);

        static Point<3> xa[129];
        static Vec<3> na[129];
        
        for (int hi = 0; hi < seia.Size(); hi++)
          {
            SurfaceElementIndex sei = seia[hi];
            const Element2dRef & el = (*mesh)[sei];

            bool drawel = (!el.IsDeleted() && el.IsVisible());

#ifdef STLGEOM
            if (checkvicinity)
              for (int j = 0; j < el.GetNP(); j++)
                if (!stlgeometry->Vicinity(el.GeomInfoPi(j+1).trignum))
                  drawel = 0;
#endif

            if (!drawel)
              continue;
            
            if (!sel_init(sei))
              continue;

            switch (el.GetType())
              {
              case TRIG:
                {
                  if (curv.IsHighOrder()) //  && curv.IsCurved(sei))
                    {
                      if (hoplotn > 128) hoplotn = 128;

                      for (int i = 0; i < hoplotn; i++)
                        {
                          glBegin (GL_TRIANGLE_STRIP);

                          for (int j = 0; j <= hoplotn-i; j++)
                            for (int k = 0; k < 2; k++)
                              {
                                if (j == hoplotn-i && k == 1) continue;

                                if (i > 0 && k == 0)
                                  {
                                    glNormal3dv (na[j]);
                                    glVertex3dv (xa[j]);
                                    continue;
                                  }

                                Point<2> xref (double(j) / hoplotn, double(i+k) / hoplotn);
                                Point<3> xglob;
                                Mat<3,2> dxdxi;
                                Vec<3> dx, dy, n;

                                curv.CalcSurfaceTransformation (xref, sei, xglob, dxdxi);
                                for (int i = 0; i < 3; i++)
                                  {
                                    dx(i) = dxdxi(i,0);
                                    dy(i) = dxdxi(i,1);
                                  }
                                n = Cross (dx, dy);
                                glNormal3dv (n);
                                glVertex3dv (xglob);

                                if (k == 1)
                                  {
                                    na[j] = n;
                                    xa[j] = xglob;
                                  }
                              }
                          glEnd();
                        }
                    }
                  else // not high order
                    {
                      glBegin (GL_TRIANGLES);
                      
                      const Point<3> & lp0 = (*mesh) [el[0]];
                      const Point<3> & lp1 = (*mesh) [el[1]];
                      const Point<3> & lp2 = (*mesh) [el[2]];

                      Vec<3> n = Cross (lp1-lp0, lp2-lp0).Normalize();
                      glNormal3dv (n);

                      for (int j = 0; j < 3; j++)
                          glVertex3dv ( (*mesh)[el[j]] );
                      
                      glEnd();
                    }
                  
                  break;
                }
              case QUAD:
                {
                  if (curv.IsHighOrder()) //  && curv.IsCurved(sei))
                    {
                      Point<2> xr[4];
                      Point<3> xg;
                      Vec<3> dx, dy, n;

                      glBegin (GL_QUADS);

                      for (int i = 0; i < hoplotn; i++)
                        for (int j = 0; j < hoplotn; j++)
                          {
                            xr[0](0) = (double)    i/hoplotn; xr[0](1) = (double)    j/hoplotn;
                            xr[1](0) = (double)(i+1)/hoplotn; xr[1](1) = (double)    j/hoplotn;
                            xr[2](0) = (double)(i+1)/hoplotn; xr[2](1) = (double)(j+1)/hoplotn;
                            xr[3](0) = (double)    i/hoplotn; xr[3](1) = (double)(j+1)/hoplotn;

                            for (int l=0; l<4; l++)
                              {
                                Mat<3,2> dxdxi;

                                curv.CalcSurfaceTransformation (xr[l], sei, xg, dxdxi);
                                for (int i = 0; i < 3; i++)
                                  {
                                    dx(i) = dxdxi(i,0);
                                    dy(i) = dxdxi(i,1);
                                  }

                                n = Cross (dx, dy);
                                n.Normalize();
                                glNormal3d (n(0), n(1), n(2));
                                glVertex3d (xg(0), xg(1), xg(2));
                              }

                          }

                      glEnd();
                    }

                  else // not high order

                    {
                      glBegin (GL_QUADS);

                      const Point<3> & lp1 = mesh->Point (el[0]);
                      const Point<3> & lp2 = mesh->Point (el[1]);
                      const Point<3> & lp3 = mesh->Point (el[3]);
                      const Point<3> & lp4 = mesh->Point (el[2]);

                      Vec<3> n = Cross (lp2-lp1,  Center (lp3, lp4)-lp1);
                      n.Normalize();
                      glNormal3dv (n);

                      glVertex3dv (lp1);
                      glVertex3dv (lp2);
                      glVertex3dv (lp4);
                      glVertex3dv (lp3);

                      glEnd ();
                    }
                  break;
                }

              case TRIG6:
                {
                  glBegin (GL_TRIANGLES);

                  static int trigs[4][3] = {
                    { 1, 6, 5 },
                    { 2, 4, 6 },
                    { 3, 5, 4 },
                    { 4, 5, 6 } };

                  for (int j = 0; j < 4; j++)
                    {
                      const Point<3> & lp1 = mesh->Point (el.PNum(trigs[j][0]));
                      const Point<3> & lp2 = mesh->Point (el.PNum(trigs[j][1]));
                      const Point<3> & lp3 = mesh->Point (el.PNum(trigs[j][2]));
                      // Vec<3> n = Cross (Vec<3> (lp1, lp2), Vec<3> (lp1, lp3));
                      Vec<3> n = Cross (lp2-lp1, lp3-lp1);
                      glNormal3dv (n);

                      glVertex3dv (lp1);
                      glVertex3dv (lp2);
                      glVertex3dv (lp3);
                    }
                  glEnd();
                  break;
                }

              case QUAD6:
                {
                  glBegin (GL_QUADS);
                  static int quads[2][4] = {
                    { 1, 5, 6, 4 },
                    { 5, 2, 3, 6 } };

                  for (int j = 0; j < 2; j++)
                    {
                      Point<3> lp1 = mesh->Point (el.PNum(quads[j][0]));
                      Point<3> lp2 = mesh->Point (el.PNum(quads[j][1]));
                      Point<3> lp3 = mesh->Point (el.PNum(quads[j][2]));
                      Point<3> lp4 = mesh->Point (el.PNum(quads[j][3]));
                      Vec<3> n = Cross (Vec<3> (lp1, lp2), Vec<3> (lp1, lp3));
                      n /= (n.Length() + 1e-12);
                      glNormal3dv (&n(0));
                      glVertex3dv (&lp1(0));
                      glVertex3dv (&lp2(0));
                      glVertex3dv (&lp3(0));
                      glVertex3dv (&lp4(0));
                    }
                  glEnd();
                  break;
                }

              case QUAD8:
                {
                  glBegin (GL_TRIANGLES);
                  static int boundary[] =
                    { 1, 5, 2, 8, 3, 6, 4, 7, 1 };

                  Point<3> c(0,0,0);
                  for (int j = 0; j < 4; j++)
                    {
                      const Point<3> & hp = mesh->Point (el[j]);
                      c(0) -= 0.25 * hp(0);
                      c(1) -= 0.25 * hp(1);
                      c(2) -= 0.25 * hp(2);
                    }
                  for (int j = 4; j < 8; j++)
                    {
                      const Point<3> & hp = mesh->Point (el[j]);
                      c(0) += 0.5 * hp(0);
                      c(1) += 0.5 * hp(1);
                      c(2) += 0.5 * hp(2);
                    }

                  for (int j = 0; j < 8; j++)
                    {
                      Point<3> lp1 = mesh->Point (el.PNum(boundary[j]));
                      Point<3> lp2 = mesh->Point (el.PNum(boundary[j+1]));

                      Vec<3> n = Cross (Vec<3> (c, lp1), Vec<3> (c, lp2));
                      n /= (n.Length() + 1e-12);
                      glNormal3dv (&n(0));
                      glVertex3dv (&lp1(0));
                      glVertex3dv (&lp2(0));
                      glVertex3dv (&c(0));
                    }
                  glEnd();
                  break;
                }


              default:
                PrintSysError ("Cannot draw (2) surface element of type ",
                               int(el.GetType()));
              }
          }
      }
  }

  void VisualSceneMesh :: BuildFilledList (bool build_select)
  {
    shared_ptr<Mesh> mesh = GetMesh();
    
    static Timer timer("Mesh::BuildFilledList");
    RegionTimer reg (timer);
    auto & list = build_select ? select.list : filledlist;
    auto & timestamp = build_select ? select.list_timestamp : filledtimestamp;
    if (list && timestamp > max(mesh->GetTimeStamp(), subdivision_timestamp))
      return;


#ifdef PARALLELGL
    if (id == 0 && ntasks > 1)
      {
        InitParallelGL();
        par_filledlists.SetSize (ntasks);

        MyMPI_SendCmd ("redraw");
        MyMPI_SendCmd ("filledlist");
        for ( int dest = 1; dest < ntasks; dest++ )
          MyMPI_Recv (par_filledlists[dest], dest, MPI_TAG_VIS);

        if (list)
          glDeleteLists (list, 1);

        list = glGenLists (1);
        glNewList (list, GL_COMPILE);

        for ( int dest = 1; dest < ntasks; dest++ )
          glCallList (par_filledlists[dest]);

        glEndList();

        timestamp = NextTimeStamp();
        return;
      }

#endif


    if (!lock)
      {
        lock = std::unique_lock<std::mutex> (mesh->Mutex());
      }

    timestamp = NextTimeStamp();

    if(!build_select && !vispar.colormeshsize)
      BuildColorTexture();

    if (list)
      glDeleteLists (list, 1);

    list = glGenLists (1);
    glNewList (list, GL_COMPILE);

    glBindTexture(GL_TEXTURE_2D, colors.texture);
      
#ifdef STLGEOM
    STLGeometry * stlgeometry = dynamic_cast<STLGeometry*> (ng_geometry);
    bool checkvicinity = (stlgeometry != NULL) && stldoctor.showvicinity;
#endif
    glEnable (GL_NORMALIZE);

    glLineWidth (1.0f);

    Array<double, PointIndex> locms;

    if (vispar.colormeshsize)
      {
        glEnable (GL_COLOR_MATERIAL);
        glShadeModel (GL_SMOOTH);
        locms.SetSize (mesh->GetNP());
        maxh = -1;
        minh = 1e99;
        for (PointIndex pi : mesh->Points().Range())
          {
            locms[pi] = mesh->GetH (mesh->Point(pi));
            if (locms[pi] > maxh) maxh = locms[pi];
            if (locms[pi] < minh) minh = locms[pi];
          }
        if (!locms.Size())
          { 
            minh = 1; 
            maxh = 10; 
          }
      }
    else if (build_select)
    {
      glDisable(GL_TEXTURE_1D);
      glDisable(GL_TEXTURE_2D);
      glDisable(GL_FOG);
      glDisable(GL_LIGHTING);
      glDisable (GL_COLOR_MATERIAL);
    }
    else
    {
      glDisable(GL_TEXTURE_1D);
      glEnable(GL_TEXTURE_2D);
      glEnable (GL_COLOR_MATERIAL);
      glBindTexture(GL_TEXTURE_2D, colors.texture);
    }

    // GLfloat matcol[] = { 0, 1, 0, 1 };
    // GLfloat matcolsel[] = { 1, 0, 0, 1 };

    GLint rendermode;
    glGetIntegerv (GL_RENDER_MODE, &rendermode);

    auto face_init = [&](int i)
      {
        if(!build_select && !vispar.colormeshsize)
        {
          float x = (0.5+i%colors.width)/colors.width;
          float y = (0.5+i/colors.width)/colors.height;
          glTexCoord2f(x,y);
        }
        return true;
      };

    CurvedElements & curv = mesh->GetCurvedElements();

    auto sel_init = [&](SurfaceElementIndex sei)
      {
        if (build_select)
          {
            GLushort r,g,b;
            r = sei.Nr1() % (1<<16);
            g = sei.Nr1() >> 16;
            b = 0;
            glColor3us(r,g,b);
          }
            
        if (vispar.colormeshsize)
        {
          auto el = (*mesh)[sei];
          if(el.GetType() == TRIG && !curv.IsHighOrder()) {
            if (vispar.colormeshsize)
              SetOpenGlColor  (locms[el[0]], minh, maxh, 0);
          }
        }
        return true;
      };

    RenderSurfaceElements(mesh, subdivisions, face_init, sel_init);

    glLoadName (0);
    glBindTexture(GL_TEXTURE_2D, 0);
    glEndList ();


#ifdef PARALLELGL
    glFinish();
    if (id > 0)
      MyMPI_Send (list, 0, MPI_TAG_VIS);
#endif
    lock = {};

  }


  void VisualSceneMesh :: BuildLineList()
  {
    shared_ptr<Mesh> mesh = GetMesh();
    if (linetimestamp > max(mesh->GetTimeStamp (), subdivision_timestamp))
      return;

    static Timer timer("Mesh::BuildLineList");
    RegionTimer reg (timer);

#ifdef PARALLELGL

    if (id == 0 && ntasks > 1)
      {
        InitParallelGL();

        par_linelists.SetSize (ntasks);

        MyMPI_SendCmd ("redraw");
        MyMPI_SendCmd ("linelist");

        for ( int dest = 1; dest < ntasks; dest++ )
          MyMPI_Recv (par_linelists[dest], dest, MPI_TAG_VIS);

        if (linelist)
          glDeleteLists (linelist, 1);

        linelist = glGenLists (1);
        glNewList (linelist, GL_COMPILE);

        for ( int dest = 1; dest < ntasks; dest++ )
          glCallList (par_linelists[dest]);

        glEndList();


        linetimestamp = NextTimeStamp();
        return;
      }

#endif

    if (!lock)
      {
        lock = std::unique_lock<std::mutex> (mesh->Mutex());
      }

    linetimestamp = NextTimeStamp();

#ifdef STLGEOM
    STLGeometry * stlgeometry = dynamic_cast<STLGeometry*> (ng_geometry);
    bool checkvicinity = (stlgeometry != NULL) && stldoctor.showvicinity;
#endif

    if (linelist)
      glDeleteLists (linelist, 1);

    linelist = glGenLists (1);
    glNewList (linelist, GL_COMPILE);

    // cout << "linelist = " << linelist << endl;

    glLineWidth (1.0f);


    int hoplotn = 1 << subdivisions;

    // PrintMessage (3, "nse = ", mesh->GetNSE());
    for (SurfaceElementIndex sei : mesh->SurfaceElements().Range())
      {
        const Element2dRef & el = (*mesh)[sei];

        bool drawel = (!el.IsDeleted() && el.IsVisible());

#ifdef STLGEOM
        if (checkvicinity)
          for (int j = 0; j < el.GetNP(); j++)
            if (!stlgeometry->Vicinity(el.GeomInfoPi(j+1).trignum))
              drawel = 0;
#endif

        if (!drawel)
          continue;

        switch (el.GetType())
          {
          case TRIG:
            {
              CurvedElements & curv = mesh->GetCurvedElements();
              if (curv.IsHighOrder()) //  && curv.IsCurved(sei))
                {
                  Point<3> xg;
                  glBegin (GL_LINE_LOOP);
                  for (int i = 0; i < hoplotn; i++)
                    {
                      Point<2> xr (double(i) / hoplotn, 0);
                      curv.CalcSurfaceTransformation (xr, sei, xg);
                      glVertex3dv (xg);
                    }
                  for (int i = 0; i < hoplotn; i++)
                    {
                      Point<2> xr (double(hoplotn-i) / hoplotn, double(i)/hoplotn);
                      curv.CalcSurfaceTransformation (xr, sei, xg);
                      glVertex3dv (xg);
                    }
                  for (int i = 0; i < hoplotn; i++)
                    {
                      Point<2> xr (0, double(hoplotn-i) / hoplotn);
                      curv.CalcSurfaceTransformation (xr, sei, xg);
                      glVertex3dv (xg);
                    }

                  glEnd();
                }
              else
                {
                  glBegin (GL_TRIANGLES);

                  for (int j = 0; j < 3; j++)
                    glVertex3dv ( (*mesh) [el[j]] );
                  /*
                  const Point<3> & lp0 = (*mesh) [el[0]];
                  const Point<3> & lp1 = (*mesh) [el[1]];
                  const Point<3> & lp2 = (*mesh) [el[2]];

                  glVertex3dv (lp0);
                  glVertex3dv (lp1);
                  glVertex3dv (lp2);
                  */
                  glEnd();
                }

              break;

            }

          case QUAD:
            {
              CurvedElements & curv = mesh->GetCurvedElements();
              if (curv.IsHighOrder()) //  && curv.IsCurved(sei))
                {
                  Point<2> xr;
                  Point<3> xg;

                  glBegin (GL_LINE_STRIP);

                  for (int side = 0; side < 4; side++)
                    {
                      for (int i = 0; i <= hoplotn; i++)
                        {
                          switch (side)
                            {
                            case 0:
                              xr(0) = (double) i/hoplotn;
                              xr(1) = 0.;
                              break;
                            case 1:
                              xr(0) = 1.;
                              xr(1) = (double) i/hoplotn;
                              break;
                            case 2:
                              xr(0) = (double) (hoplotn-i)/hoplotn;
                              xr(1) = 1.;
                              break;
                            case 3:
                              xr(0) = 0.;
                              xr(1) = (double) (hoplotn-i)/hoplotn;
                              break;
                            }

                          curv.CalcSurfaceTransformation (xr, sei, xg);
                          glVertex3d (xg(0), xg(1), xg(2));

                        }

                    }
                  glEnd();

                } else {

                glBegin (GL_QUADS);

                const Point<3> & lp1 = mesh->Point (el[0]);
                const Point<3> & lp2 = mesh->Point (el[1]);
                const Point<3> & lp3 = mesh->Point (el[3]);
                const Point<3> & lp4 = mesh->Point (el[2]);
                Vec<3> n = Cross (Vec<3> (lp1, lp2),
                                 Vec<3> (lp1, Center (lp3, lp4)));
                glNormal3d (n(0), n(1), n(2));
                glVertex3d (lp1(0), lp1(1), lp1(2));
                glVertex3d (lp2(0), lp2(1), lp2(2));
                glVertex3d (lp4(0), lp4(1), lp4(2));
                glVertex3d (lp3(0), lp3(1), lp3(2));
                glEnd();

              }

              break;

            }

          case TRIG6:
            {
              int lines[6][2] = {
                { 1, 6 }, { 2, 6 },
                { 1, 5 }, { 3, 5 },
                { 2, 4 }, { 3, 4 } };

              glBegin (GL_LINES);
              for (int j = 0; j < 6; j++)
                {
                  const Point<3> & lp1 = mesh->Point (el.PNum(lines[j][0]));
                  const Point<3> & lp2 = mesh->Point (el.PNum(lines[j][1]));

                  glVertex3d (lp1(0), lp1(1), lp1(2));
                  glVertex3d (lp2(0), lp2(1), lp2(2));
                }

              glEnd();
              break;
            }

          case QUAD6:
            {
              int lines[6][2] = {
                { 1, 5 }, { 2, 5 },
                { 3, 6 }, { 4, 6 },
                { 1, 4 }, { 2, 3 } };

              glBegin (GL_LINES);

              for (int j = 0; j < 6; j++)
                {
                  const Point<3> & lp1 = mesh->Point (el.PNum(lines[j][0]));
                  const Point<3> & lp2 = mesh->Point (el.PNum(lines[j][1]));

                  glVertex3d (lp1(0), lp1(1), lp1(2));
                  glVertex3d (lp2(0), lp2(1), lp2(2));
                }
              glEnd ();
              break;
            }

          case QUAD8:
            {
              int lines[8][2] = {
                { 1, 5 }, { 2, 5 }, { 3, 6 }, { 4, 6 },
                { 1, 7 }, { 4, 7 }, { 2, 8 }, { 3, 8 }
              };

              glBegin (GL_LINES);

              for (int j = 0; j < 8; j++)
                {
                  const Point<3> & lp1 = mesh->Point (el.PNum(lines[j][0]));
                  const Point<3> & lp2 = mesh->Point (el.PNum(lines[j][1]));

                  glVertex3d (lp1(0), lp1(1), lp1(2));
                  glVertex3d (lp2(0), lp2(1), lp2(2));
                }
              glEnd ();
              break;
            }



          default:
            PrintSysError ("Cannot draw (4) surface element of type ",
                           int(el.GetType()));
          }
      }

    glEndList ();


#ifdef PARALLELGL
    glFinish();
    if (id > 0)
      MyMPI_Send (linelist, 0, MPI_TAG_VIS);
#endif
  }



  void VisualSceneMesh :: BuildEdgeList()
  {
    shared_ptr<Mesh> mesh = GetMesh();

    if (!lock)
      {
        lock = std::unique_lock<std::mutex> (mesh->Mutex());
      }

    if (edgetimestamp > max(mesh->GetTimeStamp(), subdivision_timestamp) && vispar.drawtetsdomain == 0
        && vispar.shrink == 1)
      return;

    edgetimestamp = NextTimeStamp();

    if (edgelist)
      glDeleteLists (edgelist, 1);

    edgelist = glGenLists (1);
    glNewList (edgelist, GL_COMPILE);


    GLfloat matcoledge[] = { 0, 0, 1, 1 };
    GLfloat matcolsingedge[] = { 1, 0, 1, 1 };

    glEnable (GL_POLYGON_OFFSET_LINE);
    glPolygonOffset (1, -1);

    glEnable (GL_COLOR_MATERIAL);
    glDisable (GL_LIGHTING);

    for (SegmentIndex i : mesh->LineSegments().Range())
      {
        const Segment & seg = (*mesh)[i];

        /*
#ifdef PARALLEL
        if (ntasks > 1 && 
            vispar.drawtetsdomain && 
            // (vispar.drawtetsdomain != seg.GetPartition())) continue;
            (vispar.drawtetsdomain != mesh->seg_partition[i-1]) continue;
#endif
        */
        
        const Point<3> & p1 = (*mesh)[seg[0]];
        const Point<3> & p2 = (*mesh)[seg[1]];

        auto & ed = mesh->GetEdgeDescriptor(seg.GetIndex());
        if (ed.SingEdgeLeft() || ed.SingEdgeRight())
          glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,
                        matcolsingedge);
        else
          glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,
                        matcoledge);

        if (ed.SingEdgeLeft() || ed.SingEdgeRight())
          glColor3fv (matcolsingedge);
        else
          glColor3fv (matcoledge);

        if (mesh->GetEdgeDescriptor(seg).EdgeNr() == seledge)
          glLineWidth(5);
        else
          glLineWidth(2);

        if (mesh->GetCurvedElements().IsHighOrder())
          {
            int hoplotn = 1 << subdivisions;
            // mesh->GetCurvedElements().GetNVisualSubsecs();

            Point<3> x;
            glBegin (GL_LINE_STRIP);

            for (int j = 0; j <= hoplotn; j++)
              {
                mesh->GetCurvedElements().CalcSegmentTransformation ((double) j/hoplotn, i, x);
                glVertex3d (x(0), x(1), x(2));
                /*
                  cout << "x = " << x(0) << ", " << x(1) << ", " << x(2)
                  << ", norm = 1+" << sqrt(x(0)*x(0)+x(1)*x(1))-1
                  << ", phi = " << atan2(x(1), x(0))/M_PI << endl;
                */
              }

            glEnd();

          }
        else
          {
            glBegin (GL_LINES);
            Point<3> hp1 = p1;
            Point<3> hp2 = p2;
            Point<3> c = Center(p1, p2);
            if (vispar.shrink < 1)
              {
                hp1 = c + vispar.shrink * (hp1 - c);
                hp2 = c + vispar.shrink * (hp2 - c);
              }
            glVertex3dv (hp1);
            glVertex3dv (hp2); // p2.X(), p2.Y(), p2.Z());
            glEnd();
          }
      }

    glLineWidth (2);
    glDisable (GL_POLYGON_OFFSET_LINE);

    glDisable (GL_COLOR_MATERIAL);
    glEnable (GL_LIGHTING);

    glEndList();
  }




  void VisualSceneMesh :: BuildPointNumberList()
  {
    ;
  }



  // Bernstein Pol B_{n,i}(x) = n! / i! / (n-i)! (1-x)^{n-i} x^i
  static inline double Bernstein (int n, int i, double x)
  {
    double val = 1;
    for (int j = 0; j < i; j++)
      val *= x;
    for (int j = 1; j <= n-i; j++)
      val *= (1-x) * (j+i) / j;
    return val;
  }

  void ToBernstein (int order, Point<3> * pts, int stride)
  {
    static DenseMatrix mat, inv;
    static Vector vec1, vec2;

    if (mat.Height () != order+1)
      {
        mat.SetSize (order+1);
        inv.SetSize (order+1);
        vec1.SetSize (order+1);
        vec2.SetSize (order+1);
        for (int i = 0; i <= order; i++)
          {
            double x = double(i) / order;
            for (int j = 0; j <= order; j++)
              mat(i,j) = Bernstein (order, j, x);
          }

        CalcInverse (mat, inv);
      }

    for (int i = 0; i < 3; i++)
      {
        for (int j = 0; j <= order; j++)
          vec1(j) = pts[j*stride](i);

        inv.Mult (vec1, vec2);

        for (int j = 0; j <= order; j++)
          pts[j*stride](i) = vec2(j);
      }
  }














  void VisualSceneMesh :: BuildTetList(const TBitArray<PointIndex> & shownode)
  {
    shared_ptr<Mesh> mesh = GetMesh();

    if (tettimestamp > mesh->GetTimeStamp () &&
        tettimestamp > vispar.clipping.timestamp )
      return;

    if (!lock)
      {
        lock = std::unique_lock<std::mutex> (mesh->Mutex());
      }

    tettimestamp = NextTimeStamp();

    if (tetlist)
      glDeleteLists (tetlist, 1);


    tetlist = glGenLists (1);
    glNewList (tetlist, GL_COMPILE);


    Array<double, PointIndex> locms;

    // Philippose - 16/02/2010
    // Add Mesh size based coloring of 
    // meshes also for the volume elements
    if (vispar.colormeshsize)
      {
        glEnable (GL_COLOR_MATERIAL);
        locms.SetSize (mesh->GetNP());
        maxh = -1;
        minh = 1e99;
        for (PointIndex pi : mesh->Points().Range())
          {
            locms[pi] = mesh->GetH (mesh->Point(pi));
            if (locms[pi] > maxh) maxh = locms[pi];
            if (locms[pi] < minh) minh = locms[pi];
          }
        if (!locms.Size())
          { 
            minh = 1; 
            maxh = 10; 
          }
      }
    else
      glDisable (GL_COLOR_MATERIAL);



    Array<ElementFace> faces;

    static float tetcols[][4] =
      {
        { 1.0f, 1.0f, 0.0f, 1.0f },
        { 1.0f, 0.0f, 0.0f, 1.0f },
        { 0.0f, 1.0f, 0.0f, 1.0f },
        { 0.0f, 0.0f, 1.0f, 1.0f }
        /*
        { 1.0f, 1.0f, 0.0f, 0.3f },
        { 1.0f, 0.0f, 0.0f, 0.3f },
        { 0.0f, 1.0f, 0.0f, 0.3f },
        { 0.0f, 0.0f, 1.0f, 0.3f }
        */
      };

    CurvedElements & curv = mesh->GetCurvedElements();


    if (!curv.IsHighOrder())
      glShadeModel (GL_FLAT);
    else
      glShadeModel (GL_SMOOTH);

    int hoplotn = max (2, 1 << subdivisions);



    for (ElementIndex ei : mesh->VolumeElements().Range())
      {
        if (vispar.drawtetsdomain > 0)
          {
            /*
            int tetid = vispar.drawmetispartition ? 
              (*mesh)[ei].GetPartition() : (*mesh)[ei].GetIndex();
            */
            int tetid =  (*mesh)[ei].GetIndex().Nr1();
            if (vispar.drawtetsdomain != tetid) continue;
          }

        auto el = (*mesh)[ei];

        if ((el.GetType() == TET || el.GetType() == TET10) && !el.IsDeleted())
          {
            bool visible = true;
            for (auto pi: el.PNums())
              if (!shownode[pi])
                visible = false;
            if(!visible) continue;

            int ind = el.GetIndex().Nr1() % 4;

            // if (vispar.drawmetispartition && el.GetPartition()!=-1)
            // ind = el.GetPartition() % 4;

            glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, tetcols[ind]);


            if (curv.IsHighOrder()) //  && curv.IsCurved(ei))
              {
                const ELEMENT_FACE * faces = MeshTopology :: GetFaces1 (TET);
                const Point<3> * vertices = MeshTopology :: GetVertices (TET);

                /*
                  Point<3> grid[11][11];
                  Point<3> fpts[3];
                  int order = vispar.subdivisions+1;

                  for (int trig = 0; trig < 4; trig++)
                  {
                  for (int j = 0; j < 3; j++)
                  fpts[j] = vertices[faces[trig][j]-1];

                  static Point<3> c(0.25, 0.25, 0.25);
                  if (vispar.shrink < 1)
                  for (int j = 0; j < 3; j++)
                  fpts[j] += (1-vispar.shrink) * (c-fpts[j]);

                  for (int ix = 0; ix <= order; ix++)
                  for (int iy = 0; iy <= order; iy++)
                  {
                  double lami[3] =
                  { (1-double(ix)/order) * (1-double(iy)/order),
                  (  double(ix)/order) * (1-double(iy)/order),
                  double(iy)/order };

                  Point<3> xl;
                  for (int l = 0; l < 3; l++)
                  xl(l) = lami[0] * fpts[0](l) + lami[1] * fpts[1](l) +
                  lami[2] * fpts[2](l);

                  curv.CalcElementTransformation (xl, i-1, grid[ix][iy]);
                  }

                  for (int j = 0; j <= order; j++)
                  ToBernstein (order, &grid[j][0], &grid[0][1]-&grid[0][0]);
                  for (int j = 0; j <= order; j++)
                  ToBernstein (order, &grid[0][j], &grid[1][0]-&grid[0][0]);

                  glMap2d(GL_MAP2_VERTEX_3,
                  0.0, 1.0, &grid[0][1](0)-&grid[0][0](0), order+1,
                  0.0, 1.0, &grid[1][0](0)-&grid[0][0](0), order+1,
                  &grid[0][0](0));
                  glEnable(GL_MAP2_VERTEX_3);
                  glEnable(GL_AUTO_NORMAL);

                  glMapGrid2f(8, 0.0, 0.999, 8, 0.0, 1.0);
                  glEvalMesh2(GL_FILL, 0, 8, 0, 8);

                  glDisable (GL_AUTO_NORMAL);
                  glDisable (GL_MAP2_VERTEX_3);
                  }
                */



                int order = curv.GetOrder();

                Array<Point<3> > ploc ( (order+1)*(order+1) );
                Array<Point<3> > pglob ( (order+1)*(order+1) );
                Point<3> fpts[3];

                for (int trig = 0; trig < 4; trig++)
                  {
                    for (int j = 0; j < 3; j++)
                      fpts[j] = vertices[faces[trig][j]-1];

                    static Point<3> c(0.25, 0.25, 0.25);
                    if (vispar.shrink < 1)
                      for (int j = 0; j < 3; j++)
                        fpts[j] += (1-vispar.shrink) * (c-fpts[j]);

                    for (int ix = 0, ii = 0; ix <= order; ix++)
                      for (int iy = 0; iy <= order; iy++, ii++)
                        {
                          double lami[3] =
                            { (1-double(ix)/order) * (1-double(iy)/order),
                              (  double(ix)/order) * (1-double(iy)/order),
                              double(iy)/order };

                          Point<3> xl;
                          for (int l = 0; l < 3; l++)
                            xl(l) = lami[0] * fpts[0](l) + lami[1] * fpts[1](l) +
                              lami[2] * fpts[2](l);

                          ploc[ii] = xl;
                        }

                    curv.CalcMultiPointElementTransformation (&ploc, ei, &pglob, 0);

                    Point<3> grid[11][11];
                    for (int ix = 0, ii = 0; ix <= order; ix++)
                      for (int iy = 0; iy <= order; iy++, ii++)
                        grid[ix][iy] = pglob[ii];

                    for (int j = 0; j <= order; j++)
                      ToBernstein (order, &grid[j][0], &grid[0][1]-&grid[0][0]);
                    for (int j = 0; j <= order; j++)
                      ToBernstein (order, &grid[0][j], &grid[1][0]-&grid[0][0]);

                    glMap2d(GL_MAP2_VERTEX_3,
                            0.0, 1.0, &grid[0][1](0)-&grid[0][0](0), order+1,
                            0.0, 1.0, &grid[1][0](0)-&grid[0][0](0), order+1,
                            &grid[0][0](0));
                    glEnable(GL_MAP2_VERTEX_3);
                    glEnable(GL_AUTO_NORMAL);

                    glMapGrid2f(hoplotn, 0.0, 0.9999f, hoplotn, 0.0, 1.0);
                    glEvalMesh2(GL_FILL, 0, hoplotn, 0, hoplotn);

                    glDisable (GL_AUTO_NORMAL);
                    glDisable (GL_MAP2_VERTEX_3);
                  }
              }

            else // Not High Order

              {
                Point<3> pts[4];
                for (int j = 0; j < 4; j++)
                  pts[j] = (*mesh)[el[j]];

                if (vispar.shrink < 1)
                  {
                    Point<3> c = Center (pts[0], pts[1], pts[2], pts[3]);
                    for (int j = 0; j < 4; j++)
                      pts[j] = c + vispar.shrink * (pts[j]-c);
                  }


                Vec<3> n;


                // Philippose - 16/02/2010
                // Add Mesh size based coloring of 
                // meshes also for the volume elements
                if(vispar.colormeshsize)
                  {
                    glBegin (GL_TRIANGLE_STRIP);
                    n = Cross (pts[1]-pts[0], pts[2]-pts[0]);
                    glNormal3dv (n);

                    SetOpenGlColor (locms[el[0]], minh, maxh, 0);
                    glVertex3dv (pts[0]);

                    SetOpenGlColor (locms[el[1]], minh, maxh, 0);
                    glVertex3dv (pts[1]);

                    SetOpenGlColor (locms[el[2]], minh, maxh, 0);
                    glVertex3dv (pts[2]);

                    n = Cross (pts[3]-pts[1], pts[2]-pts[1]);
                    glNormal3dv (n);

                    SetOpenGlColor (locms[el[3]], minh, maxh, 0);
                    glVertex3dv (pts[3]);

                    n = Cross (pts[3]-pts[2], pts[0]-pts[2]);
                    glNormal3dv (n);

                    SetOpenGlColor (locms[el[0]], minh, maxh, 0);
                    glVertex3dv (pts[0]);

                    n = Cross (pts[1]-pts[3], pts[0]-pts[3]);
                    glNormal3dv (n);

                    SetOpenGlColor (locms[el[1]], minh, maxh, 0);
                    glVertex3dv (pts[1]);
                    glEnd();
                  }
                else // Do not color mesh based on mesh size
                  {
                    GLubyte ind[4][3] = { { 0,1,2 }, { 3,1,0 },
                                          { 1,3,2 }, { 2,3,0 } };
                    
                    glEnableClientState(GL_VERTEX_ARRAY);
                    glVertexPointer(3, GL_DOUBLE, 0, &pts[0](0));

                    for (int j = 0; j < 4; j++)
                      { 
                        glNormal3dv (Cross (pts[ind[j][1]]-pts[ind[j][0]],
                                            pts[ind[j][2]]-pts[ind[j][0]]));

                        glDrawElements(GL_TRIANGLES, 3, GL_UNSIGNED_BYTE, &ind[j][0]);
                      }
                    glDisableClientState(GL_VERTEX_ARRAY);

                    /*
                    glBegin (GL_TRIANGLE_STRIP);
                    glNormal3dv (Cross (pts[1]-pts[0], pts[2]-pts[0]));

                    glVertex3dv (pts[0]);
                    glVertex3dv (pts[1]);
                    glVertex3dv (pts[2]);

                    glNormal3dv (Cross (pts[3]-pts[1], pts[2]-pts[1]));
                    glVertex3dv (pts[3]);

                    glNormal3dv (Cross (pts[3]-pts[2], pts[0]-pts[2]));
                    glVertex3dv (pts[0]);

                    glNormal3dv (Cross (pts[1]-pts[3], pts[0]-pts[3]));
                    glVertex3dv (pts[1]);
                    glEnd();
                    */
                  }

              }
          }
      }

    glEndList ();
  }




  void VisualSceneMesh :: BuildPrismList(const TBitArray<PointIndex> & shownode)
  {
    shared_ptr<Mesh> mesh = GetMesh();
    
    if (prismtimestamp > mesh->GetTimeStamp () &&
        prismtimestamp > vispar.clipping.timestamp )
      return;

    if (!lock)
      {
        lock = std::unique_lock<std::mutex> (mesh->Mutex());
      }

    prismtimestamp = NextTimeStamp();



    if (prismlist)
      glDeleteLists (prismlist, 1);

    prismlist = glGenLists (1);
    glNewList (prismlist, GL_COMPILE);

    static float prismcol[] = { 0.0f, 1.0f, 1.0f, 1.0f };
    glLineWidth (1.0f);

    Array<ElementFace> faces;


    glDisable (GL_COLOR_MATERIAL);
    glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, prismcol);

    for (ElementIndex ei : mesh->VolumeElements().Range())
      {
        auto el = (*mesh)[ei];
        if (el.GetType() == PRISM && !el.IsDeleted())
          {
            bool visible = true;
            for (auto pi: el.PNums())
              if (!shownode[pi])
                visible = false;
            if(!visible) continue;

            int j;

            CurvedElements & curv = mesh->GetCurvedElements();
            if (curv.IsHighOrder()) //  && curv.IsCurved(ei))
              {
                const ELEMENT_FACE * faces = MeshTopology :: GetFaces1 (PRISM);
                const Point<3> * vertices = MeshTopology :: GetVertices (PRISM);

                Point<3> grid[11][11];
                Point<3> fpts[4];
                int order = subdivisions+1;

                for (int trig = 0; trig < 2; trig++)
                  {
                    for (int j = 0; j < 3; j++)
                      fpts[j] = vertices[faces[trig][j]-1];

                    static Point<3> c(1.0/3.0, 1.0/3.0, 0.5);
                    if (vispar.shrink < 1)
                      for (int j = 0; j < 3; j++)
                        fpts[j] += (1-vispar.shrink) * (c-fpts[j]);

                    for (int ix = 0; ix <= order; ix++)
                      for (int iy = 0; iy <= order; iy++)
                        {
                          double lami[3] =
                            { (1-double(ix)/order) * (1-double(iy)/order),
                              (  double(ix)/order) * (1-double(iy)/order),
                              double(iy)/order };

                          Point<3> xl;
                          for (int l = 0; l < 3; l++)
                            xl(l) = lami[0] * fpts[0](l) + lami[1] * fpts[1](l) +
                              lami[2] * fpts[2](l);

                          curv.CalcElementTransformation (xl, ei, grid[ix][iy]);
                        }

                    for (int j = 0; j <= order; j++)
                      ToBernstein (order, &grid[j][0], &grid[0][1]-&grid[0][0]);
                    for (int j = 0; j <= order; j++)
                      ToBernstein (order, &grid[0][j], &grid[1][0]-&grid[0][0]);

                    glMap2d(GL_MAP2_VERTEX_3,
                            0.0, 1.0, &grid[0][1](0)-&grid[0][0](0), order+1,
                            0.0, 1.0, &grid[1][0](0)-&grid[0][0](0), order+1,
                            &grid[0][0](0));
                    glEnable(GL_MAP2_VERTEX_3);
                    glEnable(GL_AUTO_NORMAL);

                    glMapGrid2f(8, 0.0, 0.999f, 8, 0.0, 1.0);
                    glEvalMesh2(GL_FILL, 0, 8, 0, 8);

                    glDisable (GL_AUTO_NORMAL);
                    glDisable (GL_MAP2_VERTEX_3);
                  }

                for (int quad = 2; quad < 5; quad++)
                  {
                    for (int j = 0; j < 4; j++)
                      fpts[j] = vertices[faces[quad][j]-1];

                    static Point<3> c(1.0/3.0, 1.0/3.0, 0.5);
                    if (vispar.shrink < 1)
                      for (int j = 0; j < 4; j++)
                        fpts[j] += (1-vispar.shrink) * (c-fpts[j]);

                    for (int ix = 0; ix <= order; ix++)
                      for (int iy = 0; iy <= order; iy++)
                        {
                          double lami[4] =
                            { (1-double(ix)/order) * (1-double(iy)/order),
                              (  double(ix)/order) * (1-double(iy)/order),
                              (  double(ix)/order) * (  double(iy)/order),
                              (1-double(ix)/order) * (  double(iy)/order) };

                          Point<3> xl;
                          for (int l = 0; l < 3; l++)
                            xl(l) =
                              lami[0] * fpts[0](l) + lami[1] * fpts[1](l) +
                              lami[2] * fpts[2](l) + lami[3] * fpts[3](l);

                          curv.CalcElementTransformation (xl, ei, grid[ix][iy]);
                        }

                    for (int j = 0; j <= order; j++)
                      ToBernstein (order, &grid[j][0], &grid[0][1]-&grid[0][0]);
                    for (int j = 0; j <= order; j++)
                      ToBernstein (order, &grid[0][j], &grid[1][0]-&grid[0][0]);

                    glMap2d(GL_MAP2_VERTEX_3,
                            0.0, 1.0, &grid[0][1](0)-&grid[0][0](0), order+1,
                            0.0, 1.0, &grid[1][0](0)-&grid[0][0](0), order+1,
                            &grid[0][0](0));
                    glEnable(GL_MAP2_VERTEX_3);
                    glEnable(GL_AUTO_NORMAL);

                    glMapGrid2f(8, 0.0, 1.0, 8, 0.0, 1.0);
                    glEvalMesh2(GL_FILL, 0, 8, 0, 8);

                    glDisable (GL_AUTO_NORMAL);
                    glDisable (GL_MAP2_VERTEX_3);
                  }





                /*
                  int hoplotn = 1 << subdivisions;
                  // int hoplotn = curv.GetNVisualSubsecs();

                  const Point<3> * facepoint = MeshTopology :: GetVertices (TRIG);
                  const ELEMENT_FACE * elface = MeshTopology :: GetFaces(TRIG);

                  glBegin (GL_TRIANGLES);

                  for (int trig = 0; trig<2; trig++)
                  {

                  Vec<3> x0,x1,d0,d1;
                  x0 = facepoint[1] - facepoint[2];
                  x1 = facepoint[0] - facepoint[2];
                  x0.Normalize();
                  x1.Normalize();
                  if (trig == 1) swap (x0,x1);

                  Point<3> xr[3];
                  Point<3> xg;
                  Vec<3> dx, dy, dz, n;

                  for (int i1 = 0; i1 < hoplotn; i1++)
                  for (int j1 = 0; j1 < hoplotn-i1; j1++)
                  for (int k = 0; k < 2; k++)
                  {
                  if (k == 0)
                  {
                  xr[0](0) = (double)    i1/hoplotn; xr[0](1) = (double)    j1/hoplotn;
                  xr[1](0) = (double)(i1+1)/hoplotn; xr[1](1) = (double)    j1/hoplotn;
                  xr[2](0) = (double)    i1/hoplotn; xr[2](1) = (double)(j1+1)/hoplotn;
                  } else
                  {
                  if (j1 == hoplotn-i1-1) continue;
                  xr[0](0) = (double)(i1+1)/hoplotn; xr[0](1) = (double)    j1/hoplotn;
                  xr[1](0) = (double)(i1+1)/hoplotn; xr[1](1) = (double)(j1+1)/hoplotn;
                  xr[2](0) = (double)    i1/hoplotn; xr[2](1) = (double)(j1+1)/hoplotn;
                  };

                  for (int l=0; l<3; l++)
                  {
                  Mat<3,3> dxdxi;
                  xr[l](2) = (double) trig;
                  curv.CalcElementTransformation (xr[l], i-1, xg, dxdxi);
                  for (int i = 0; i < 3; i++)
                  {
                  dx(i) = dxdxi(i,0);
                  dy(i) = dxdxi(i,1);
                  dz(i) = dxdxi(i,2);
                  }

                  Vec<3> d0 = x0(0)*dx + x0(1)*dy + x0(2)*dz;
                  Vec<3> d1 = x1(0)*dx + x1(1)*dy + x1(2)*dz;
                  n = Cross (d1, d0);
                  glNormal3d (n(0), n(1), n(2));
                  glVertex3d (xg(0), xg(1), xg(2));
                  }
                  }

                  }

                  glEnd ();

                  glBegin (GL_QUADS);

                  for (int quad = 0; quad<3; quad++)
                  {
                  const Point<3> * facepoint = MeshTopology :: GetVertices (PRISM);

                  Vec<3> x0,x1;
                  int xyz;

                  switch (quad)
                  {
                  case 0:
                  x0 = facepoint[5] - facepoint[2];
                  x1 = facepoint[0] - facepoint[2];
                  xyz = 0;
                  break;
                  case 1:
                  x0 = facepoint[4] - facepoint[0];
                  x1 = facepoint[1] - facepoint[0];
                  xyz = 0;
                  break;
                  case 2:
                  x0 = facepoint[1] - facepoint[2];
                  x1 = facepoint[5] - facepoint[2];
                  xyz = 1;
                  break;
                  }

                  x0.Normalize();
                  x1.Normalize();

                  swap (x0,x1);

                  Point<3> xr[4];
                  Point<3> xg;
                  Vec<3> dx, dy, dz, n;

                  for (int i1 = 0; i1 < hoplotn; i1++)
                  for (int j1 = 0; j1 < hoplotn; j1++)
                  {
                  xr[0](xyz) = (double)    i1/hoplotn; xr[0](2) = (double)    j1/hoplotn;
                  xr[1](xyz) = (double)(i1+1)/hoplotn; xr[1](2) = (double)    j1/hoplotn;
                  xr[2](xyz) = (double)(i1+1)/hoplotn; xr[2](2) = (double)(j1+1)/hoplotn;
                  xr[3](xyz) = (double)    i1/hoplotn; xr[3](2) = (double)(j1+1)/hoplotn;

                  for (int l=0; l<4; l++)
                  {
                  switch (quad)
                  {
                  case 0: xr[l](1) = 0; break;
                  case 1: xr[l](1) = 1-xr[l](0); break;
                  case 2: xr[l](0) = 0; break;
                  }

                  Mat<3,3> dxdxi;
                  curv.CalcElementTransformation (xr[l], i-1, xg, dxdxi);
                  for (int i = 0; i < 3; i++)
                  {
                  dx(i) = dxdxi(i,0);
                  dy(i) = dxdxi(i,1);
                  dz(i) = dxdxi(i,2);
                  }

                  Vec<3> d0 = x0(0)*dx + x0(1)*dy + x0(2)*dz;
                  Vec<3> d1 = x1(0)*dx + x1(1)*dy + x1(2)*dz;
                  n = Cross (d1, d0);
                  glNormal3d (n(0), n(1), n(2));
                  glVertex3d (xg(0), xg(1), xg(2));
                  }
                  }
                  }
                  glEnd ();
                */
              }
            else
              {
                Point<3> c(0,0,0);
                if (vispar.shrink < 1)
                  {
                    for (j = 1; j <= 6; j++)
                      {
                        Point<3> p = mesh->Point(el.PNum(j));
                        c(0) += p(0) / 6;
                        c(1) += p(1) / 6;
                        c(2) += p(2) / 6;
                      }
                  }

                el.GetSurfaceTriangles (faces);
                glBegin (GL_TRIANGLES);
                for (j = 1; j <= faces.Size(); j++)
                  {
                    ElementFace & face = faces[j-1];
                    Point<3> lp1 = mesh->Point (el.PNum(face[0]));
                    Point<3> lp2 = mesh->Point (el.PNum(face[1]));
                    Point<3> lp3 = mesh->Point (el.PNum(face[2]));
                    Vec<3> n = Cross (Vec<3> (lp1, lp3), Vec<3> (lp1, lp2));
                    n /= (n.Length()+1e-12);
                    glNormal3d (n(0), n(1), n(2));
                    if (vispar.shrink < 1)
                      {
                        lp1 = c + vispar.shrink * (lp1 - c);
                        lp2 = c + vispar.shrink * (lp2 - c);
                        lp3 = c + vispar.shrink * (lp3 - c);
                      }
                    glVertex3d (lp1(0), lp1(1), lp1(2));
                    glVertex3d (lp2(0), lp2(1), lp2(2));
                    glVertex3d (lp3(0), lp3(1), lp3(2));
                  }

                glEnd();
              }
          }
      }
    glEndList ();
  }




  void VisualSceneMesh :: BuildHexList(const TBitArray<PointIndex> & shownode)
  {
    shared_ptr<Mesh> mesh = GetMesh();
    
    if (hextimestamp > mesh->GetTimeStamp () &&
        hextimestamp > vispar.clipping.timestamp )
      return;

    if (!lock)
      {
        lock = std::unique_lock<std::mutex> (mesh->Mutex());
      }

    hextimestamp = NextTimeStamp();

    if (hexlist) glDeleteLists (hexlist, 1);

    hexlist = glGenLists (1);
    glNewList (hexlist, GL_COMPILE);


    static float hexcol[] = { 1.0f, 1.0f, 0.0f, 1.0f };
    glLineWidth (1.0f);
    glDisable (GL_COLOR_MATERIAL);
    glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, hexcol);

    Array<ElementFace> faces;
    // int hoplotn = 1 << vispar.subdivisions;

    for (ElementIndex ei : mesh->VolumeElements().Range())
      {
        auto el = (*mesh)[ei];
        if (el.GetType() == HEX && !el.IsDeleted())
          {
            bool visible = true;
            for (auto pi: el.PNums())
              if (!shownode[pi])
                visible = false;
            if(!visible) continue;
            CurvedElements & curv = mesh->GetCurvedElements();
            if (curv.IsHighOrder()) //  && curv.IsCurved(ei))
              {
                /* // classical
                   glBegin (GL_QUADS);

                   const ELEMENT_FACE * faces = MeshTopology :: GetFaces (HEX);
                   const Point<3> * vertices = MeshTopology :: GetVertices (HEX);

                   Point<3> grid[33][33];
                   Vec<3> gridn[33][33];
                   Point<3> fpts[4];
                   for (int quad = 0; quad<6; quad++)
                   {
                   for (int j = 0; j < 4; j++)
                   fpts[j] = vertices[faces[quad][j]-1];

                   static Point<3> c(0.5, 0.5, 0.5);
                   if (vispar.shrink < 1)
                   for (int j = 0; j < 4; j++)
                   fpts[j] += (1-vispar.shrink) * (c-fpts[j]);

                   Vec<3> taux = fpts[1]-fpts[0];
                   Vec<3> tauy = fpts[3]-fpts[0];

                   for (int ix = 0; ix <= hoplotn; ix++)
                   for (int iy = 0; iy <= hoplotn; iy++)
                   {
                   Point<3> xl;
                   Mat<3,3> dxdxi;
                   double lami[4] =
                   { (1-double(ix)/hoplotn) * (1-double(iy)/hoplotn),
                   (  double(ix)/hoplotn) * (1-double(iy)/hoplotn),
                   (  double(ix)/hoplotn) * (  double(iy)/hoplotn),
                   (1-double(ix)/hoplotn) * (  double(iy)/hoplotn) };
                   for (int l = 0; l < 3; l++)
                   xl(l) = lami[0] * fpts[0](l) + lami[1] * fpts[1](l) +
                   lami[2] * fpts[2](l) + lami[3] * fpts[3](l);

                   curv.CalcElementTransformation (xl, ei, grid[ix][iy], dxdxi);

                   Vec<3> gtaux = dxdxi * taux;
                   Vec<3> gtauy = dxdxi * tauy;
                   gridn[ix][iy] = Cross (gtauy, gtaux).Normalize();
                   }

                   for (int ix = 0; ix < hoplotn; ix++)
                   for (int iy = 0; iy < hoplotn; iy++)
                   {
                   glNormal3dv (gridn[ix][iy]);
                   glVertex3dv (grid[ix][iy]);

                   glNormal3dv (gridn[ix+1][iy]);
                   glVertex3dv (grid[ix+1][iy]);

                   glNormal3dv (gridn[ix+1][iy+1]);
                   glVertex3dv (grid[ix+1][iy+1]);

                   glNormal3dv (gridn[ix][iy+1]);
                   glVertex3dv (grid[ix][iy+1]);
                   }
                   }

                   glEnd ();
                */

                const ELEMENT_FACE * faces = MeshTopology :: GetFaces1 (HEX);
                const Point<3> * vertices = MeshTopology :: GetVertices (HEX);

                Point<3> grid[11][11];
                Point<3> fpts[4];
                int order = subdivisions+1;

                for (int quad = 0; quad<6; quad++)
                  {
                    for (int j = 0; j < 4; j++)
                      fpts[j] = vertices[faces[quad][j]-1];

                    static Point<3> c(0.5, 0.5, 0.5);
                    if (vispar.shrink < 1)
                      for (int j = 0; j < 4; j++)
                        fpts[j] += (1-vispar.shrink) * (c-fpts[j]);

                    for (int ix = 0; ix <= order; ix++)
                      for (int iy = 0; iy <= order; iy++)
                        {
                          double lami[4] =
                            { (1-double(ix)/order) * (1-double(iy)/order),
                              (  double(ix)/order) * (1-double(iy)/order),
                              (  double(ix)/order) * (  double(iy)/order),
                              (1-double(ix)/order) * (  double(iy)/order) };

                          Point<3> xl;
                          for (int l = 0; l < 3; l++)
                            xl(l) = lami[0] * fpts[0](l) + lami[1] * fpts[1](l) +
                              lami[2] * fpts[2](l) + lami[3] * fpts[3](l);

                          curv.CalcElementTransformation (xl, ei, grid[ix][iy]);
                        }

                    for (int j = 0; j <= order; j++)
                      ToBernstein (order, &grid[j][0], &grid[0][1]-&grid[0][0]);
                    for (int j = 0; j <= order; j++)
                      ToBernstein (order, &grid[0][j], &grid[1][0]-&grid[0][0]);

                    glMap2d(GL_MAP2_VERTEX_3,
                            0.0, 1.0, &grid[0][1](0)-&grid[0][0](0), order+1,
                            0.0, 1.0, &grid[1][0](0)-&grid[0][0](0), order+1,
                            &grid[0][0](0));
                    glEnable(GL_MAP2_VERTEX_3);
                    glEnable(GL_AUTO_NORMAL);

                    glMapGrid2f(8, 0.0, 1.0, 8, 0.0, 1.0);
                    glEvalMesh2(GL_FILL, 0, 8, 0, 8);

                    glDisable (GL_AUTO_NORMAL);
                    glDisable (GL_MAP2_VERTEX_3);
                  }
              }
            else
              {
                Point<3> c(0,0,0);
                if (vispar.shrink < 1)
                  {
                    for (int j = 0; j < 8; j++)
                      {
                        Point<3> p = mesh->Point(el[j]);
                        c(0) += p(0);
                        c(1) += p(1);
                        c(2) += p(2);
                      }
                    c(0) /= 8;
                    c(1) /= 8;
                    c(2) /= 8;
                  }

                glBegin (GL_TRIANGLES);

                el.GetSurfaceTriangles (faces);
                for (int j = 0; j < faces.Size(); j++)
                  {
                    ElementFace & face = faces[j];
                    Point<3> lp1 = mesh->Point (el.PNum(face[0]));
                    Point<3> lp2 = mesh->Point (el.PNum(face[1]));
                    Point<3> lp3 = mesh->Point (el.PNum(face[2]));
                    Vec<3> n = Cross (lp3-lp1, lp2-lp1);
                    n.Normalize();
                    glNormal3dv (n);

                    if (vispar.shrink < 1)
                      {
                        lp1 = c + vispar.shrink * (lp1 - c);
                        lp2 = c + vispar.shrink * (lp2 - c);
                        lp3 = c + vispar.shrink * (lp3 - c);
                      }

                    glVertex3dv (lp1);
                    glVertex3dv (lp2);
                    glVertex3dv (lp3);
                  }

                glEnd();
              }
          }
      }

    static float hex7col[] = { 1.0f, 0.65f, 0.0f, 1.0f };
    glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, hex7col);

    for (auto el : mesh->VolumeElements())
      {
        if (el.GetType() == HEX7 && !el.IsDeleted())
          {
            /*
            CurvedElements & curv = mesh->GetCurvedElements();
            if (curv.IsHighOrder()) 
              {
                const ELEMENT_FACE * faces = MeshTopology :: GetFaces1 (HEX);
                const Point<3> * vertices = MeshTopology :: GetVertices (HEX);

                Point<3> grid[11][11];
                Point<3> fpts[4];
                int order = subdivisions+1;

                for (int quad = 0; quad<6; quad++)
                  {
                    for (int j = 0; j < 4; j++)
                      fpts[j] = vertices[faces[quad][j]-1];

                    static Point<3> c(0.5, 0.5, 0.5);
                    if (vispar.shrink < 1)
                      for (int j = 0; j < 4; j++)
                        fpts[j] += (1-vispar.shrink) * (c-fpts[j]);

                    for (int ix = 0; ix <= order; ix++)
                      for (int iy = 0; iy <= order; iy++)
                        {
                          double lami[4] =
                            { (1-double(ix)/order) * (1-double(iy)/order),
                              (  double(ix)/order) * (1-double(iy)/order),
                              (  double(ix)/order) * (  double(iy)/order),
                              (1-double(ix)/order) * (  double(iy)/order) };

                          Point<3> xl;
                          for (int l = 0; l < 3; l++)
                            xl(l) = lami[0] * fpts[0](l) + lami[1] * fpts[1](l) +
                              lami[2] * fpts[2](l) + lami[3] * fpts[3](l);

                          curv.CalcElementTransformation (xl, ei, grid[ix][iy]);
                        }

                    for (int j = 0; j <= order; j++)
                      ToBernstein (order, &grid[j][0], &grid[0][1]-&grid[0][0]);
                    for (int j = 0; j <= order; j++)
                      ToBernstein (order, &grid[0][j], &grid[1][0]-&grid[0][0]);

                    glMap2d(GL_MAP2_VERTEX_3,
                            0.0, 1.0, &grid[0][1](0)-&grid[0][0](0), order+1,
                            0.0, 1.0, &grid[1][0](0)-&grid[0][0](0), order+1,
                            &grid[0][0](0));
                    glEnable(GL_MAP2_VERTEX_3);
                    glEnable(GL_AUTO_NORMAL);

                    glMapGrid2f(8, 0.0, 1.0, 8, 0.0, 1.0);
                    glEvalMesh2(GL_FILL, 0, 8, 0, 8);

                    glDisable (GL_AUTO_NORMAL);
                    glDisable (GL_MAP2_VERTEX_3);
                  }
              }
            else
            */
              {
                Point<3> c(0,0,0);
                if (vispar.shrink < 1)
                  {
                    for (int j = 0; j < 7; j++)
                      {
                        Point<3> p = mesh->Point(el[j]);
                        c(0) += p(0);
                        c(1) += p(1);
                        c(2) += p(2);
                      }
                    c(0) /= 7;
                    c(1) /= 7;
                    c(2) /= 7;
                  }

                glBegin (GL_TRIANGLES);

                el.GetSurfaceTriangles (faces);
                for (int j = 0; j < faces.Size(); j++)
                  {
                    ElementFace & face = faces[j];
                    Point<3> lp1 = mesh->Point (el.PNum(face[0]));
                    Point<3> lp2 = mesh->Point (el.PNum(face[1]));
                    Point<3> lp3 = mesh->Point (el.PNum(face[2]));
                    Vec<3> n = Cross (lp3-lp1, lp2-lp1);
                    n.Normalize();
                    glNormal3dv (n);

                    if (vispar.shrink < 1)
                      {
                        lp1 = c + vispar.shrink * (lp1 - c);
                        lp2 = c + vispar.shrink * (lp2 - c);
                        lp3 = c + vispar.shrink * (lp3 - c);
                      }

                    glVertex3dv (lp1);
                    glVertex3dv (lp2);
                    glVertex3dv (lp3);
                  }

                glEnd();
              }
          }
      }

    
    glEndList ();
  }









  void VisualSceneMesh :: BuildPyramidList(const TBitArray<PointIndex> & shownode)
  {
    shared_ptr<Mesh> mesh = GetMesh();
    
    if (pyramidtimestamp > mesh->GetTimeStamp () &&
        pyramidtimestamp > vispar.clipping.timestamp )
      return;

    if (!lock)
      {
        lock = std::unique_lock<std::mutex> (mesh->Mutex());
      }

    pyramidtimestamp = NextTimeStamp();


    if (pyramidlist)
      glDeleteLists (pyramidlist, 1);


    pyramidlist = glGenLists (1);
    glNewList (pyramidlist, GL_COMPILE);

    static float pyramidcol[] = { 1.0f, 0.0f, 1.0f, 1.0f };
    glDisable (GL_COLOR_MATERIAL);
    glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, pyramidcol);

    glLineWidth (1.0f);
    Array<ElementFace> faces;

    for (ElementIndex ei : mesh->VolumeElements().Range())
      {
        auto el = (*mesh)[ei];
        if ((el.GetType() == PYRAMID || el.GetType() == PYRAMID13) && !el.IsDeleted())
          {
            bool visible = true;
            for (auto pi: el.PNums())
              if (!shownode[pi])
                visible = false;
            if(!visible) continue;

            CurvedElements & curv = mesh->GetCurvedElements();
            if (curv.IsHighOrder()) //  && curv.IsCurved(ei))
              {

                const ELEMENT_FACE * faces = MeshTopology :: GetFaces1 (PYRAMID);
                const Point<3> * vertices = MeshTopology :: GetVertices (PYRAMID);

                Point<3> grid[11][11];
                Point<3> fpts[4];
                int order = subdivisions+1;

                for (int trig = 0; trig < 4; trig++)
                  {
                    for (int j = 0; j < 3; j++)
                      fpts[j] = vertices[faces[trig][j]-1];

                    static Point<3> c(0.375, 0.375, 0.25);
                    if (vispar.shrink < 1)
                      for (int j = 0; j < 3; j++)
                        fpts[j] += (1-vispar.shrink) * (c-fpts[j]);

                    for (int ix = 0; ix <= order; ix++)
                      for (int iy = 0; iy <= order; iy++)
                        {
                          double lami[3] =
                            { (1-double(ix)/order) * (1-double(iy)/order),
                              (  double(ix)/order) * (1-double(iy)/order),
                              double(iy)/order };

                          Point<3> xl;
                          for (int l = 0; l < 3; l++)
                            xl(l) = lami[0] * fpts[0](l) + lami[1] * fpts[1](l) +
                              lami[2] * fpts[2](l);

                          curv.CalcElementTransformation (xl, ei, grid[ix][iy]);
                        }

                    for (int j = 0; j <= order; j++)
                      ToBernstein (order, &grid[j][0], &grid[0][1]-&grid[0][0]);
                    for (int j = 0; j <= order; j++)
                      ToBernstein (order, &grid[0][j], &grid[1][0]-&grid[0][0]);

                    glMap2d(GL_MAP2_VERTEX_3,
                            0.0, 1.0, &grid[0][1](0)-&grid[0][0](0), order+1,
                            0.0, 1.0, &grid[1][0](0)-&grid[0][0](0), order+1,
                            &grid[0][0](0));
                    glEnable(GL_MAP2_VERTEX_3);
                    glEnable(GL_AUTO_NORMAL);

                    glMapGrid2f(8, 0.0, 0.999f, 8, 0.0, 1.0);
                    glEvalMesh2(GL_FILL, 0, 8, 0, 8);

                    glDisable (GL_AUTO_NORMAL);
                    glDisable (GL_MAP2_VERTEX_3);
                  }

                for (int quad = 4; quad < 5; quad++)
                  {
                    for (int j = 0; j < 4; j++)
                      fpts[j] = vertices[faces[quad][j]-1];

                    static Point<3> c(0.375, 0.375, 0.25);
                    if (vispar.shrink < 1)
                      for (int j = 0; j < 4; j++)
                        fpts[j] += (1-vispar.shrink) * (c-fpts[j]);

                    for (int ix = 0; ix <= order; ix++)
                      for (int iy = 0; iy <= order; iy++)
                        {
                          double lami[4] =
                            { (1-double(ix)/order) * (1-double(iy)/order),
                              (  double(ix)/order) * (1-double(iy)/order),
                              (  double(ix)/order) * (  double(iy)/order),
                              (1-double(ix)/order) * (  double(iy)/order) };

                          Point<3> xl;
                          for (int l = 0; l < 3; l++)
                            xl(l) =
                              lami[0] * fpts[0](l) + lami[1] * fpts[1](l) +
                              lami[2] * fpts[2](l) + lami[3] * fpts[3](l);

                          curv.CalcElementTransformation (xl, ei, grid[ix][iy]);
                        }

                    for (int j = 0; j <= order; j++)
                      ToBernstein (order, &grid[j][0], &grid[0][1]-&grid[0][0]);
                    for (int j = 0; j <= order; j++)
                      ToBernstein (order, &grid[0][j], &grid[1][0]-&grid[0][0]);

                    glMap2d(GL_MAP2_VERTEX_3,
                            0.0, 1.0, &grid[0][1](0)-&grid[0][0](0), order+1,
                            0.0, 1.0, &grid[1][0](0)-&grid[0][0](0), order+1,
                            &grid[0][0](0));
                    glEnable(GL_MAP2_VERTEX_3);
                    glEnable(GL_AUTO_NORMAL);

                    glMapGrid2f(8, 0.0, 1.0, 8, 0.0, 1.0);
                    glEvalMesh2(GL_FILL, 0, 8, 0, 8);

                    glDisable (GL_AUTO_NORMAL);
                    glDisable (GL_MAP2_VERTEX_3);
                  }






                /*
                  int hoplotn = 1 << vispar.subdivisions;

                  const ELEMENT_FACE * faces = MeshTopology :: GetFaces (PYRAMID);
                  const Point<3> * vertices = MeshTopology :: GetVertices (PYRAMID);

                  Point<3> grid[33][33];
                  Vec<3> gridn[33][33];


                  glBegin (GL_TRIANGLES);

                  for (int trig = 0; trig < 4; trig++)
                  {
                  Point<3> p0 = vertices[faces[trig][0]-1];
                  Point<3> p1 = vertices[faces[trig][1]-1];
                  Point<3> p2 = vertices[faces[trig][2]-1];

                  if (vispar.shrink < 1)
                  {
                  static Point<3> c(0.375, 0.375, 0.25);
                  p0 = c + vispar.shrink * (p0 - c);
                  p1 = c + vispar.shrink * (p1 - c);
                  p2 = c + vispar.shrink * (p2 - c);
                  }


                  Vec<3> taux = p0-p2;
                  Vec<3> tauy = p1-p2;
                  Vec<3> gtaux, gtauy;

                  Point<3> xl;
                  Mat<3,3> dxdxi;

                  for (int ix = 0; ix <= hoplotn; ix++)
                  for (int iy = 0; iy <= hoplotn-ix; iy++)
                  {
                  for (int l = 0; l < 3; l++)
                  xl(l) =
                  (1-double(ix+iy)/hoplotn) * p2(l) +
                  (double(ix)/hoplotn) * p0(l) +
                  (double(iy)/hoplotn) * p1(l);

                  curv.CalcElementTransformation (xl, i-1, grid[ix][iy], dxdxi);

                  gtaux = dxdxi * taux;
                  gtauy = dxdxi * tauy;
                  gridn[ix][iy] = Cross (gtauy, gtaux).Normalize();
                  }

                  for (int ix = 0; ix < hoplotn; ix++)
                  for (int iy = 0; iy < hoplotn-ix; iy++)
                  {
                  glNormal3dv (gridn[ix][iy]);
                  glVertex3dv (grid[ix][iy]);

                  glNormal3dv (gridn[ix+1][iy]);
                  glVertex3dv (grid[ix+1][iy]);

                  glNormal3dv (gridn[ix][iy+1]);
                  glVertex3dv (grid[ix][iy+1]);

                  if (iy < hoplotn-ix-1)
                  {
                  glNormal3dv (gridn[ix][iy+1]);
                  glVertex3dv (grid[ix][iy+1]);

                  glNormal3dv (gridn[ix+1][iy]);
                  glVertex3dv (grid[ix+1][iy]);

                  glNormal3dv (gridn[ix+1][iy+1]);
                  glVertex3dv (grid[ix+1][iy+1]);
                  }
                  }
                  }

                  glEnd ();




                  glBegin (GL_QUADS);

                  for (int quad = 4; quad < 5; quad++)
                  {
                  Point<3> p0 = vertices[faces[quad][0]-1];
                  Point<3> p1 = vertices[faces[quad][1]-1];
                  Point<3> p2 = vertices[faces[quad][2]-1];
                  Point<3> p3 = vertices[faces[quad][3]-1];

                  if (vispar.shrink < 1)
                  {
                  static Point<3> c(0.375, 0.375, 0.25);
                  p0 = c + vispar.shrink * (p0 - c);
                  p1 = c + vispar.shrink * (p1 - c);
                  p2 = c + vispar.shrink * (p2 - c);
                  p3 = c + vispar.shrink * (p3 - c);
                  }

                  Vec<3> taux = p1-p0;
                  Vec<3> tauy = p3-p0;
                  Vec<3> gtaux, gtauy;

                  Point<3> xl, xg;
                  Mat<3,3> dxdxi;

                  for (int ix = 0; ix <= hoplotn; ix++)
                  for (int iy = 0; iy <= hoplotn; iy++)
                  {
                  Point<3> xl;
                  for (int l = 0; l < 3; l++)
                  xl(l) =
                  (1-double(ix)/hoplotn)*(1-double(iy)/hoplotn) * p0(l) +
                  (  double(ix)/hoplotn)*(1-double(iy)/hoplotn) * p1(l) +
                  (  double(ix)/hoplotn)*(  double(iy)/hoplotn) * p2(l) +
                  (1-double(ix)/hoplotn)*(  double(iy)/hoplotn) * p3(l);

                  curv.CalcElementTransformation (xl, i-1, grid[ix][iy], dxdxi);

                  gtaux = dxdxi * taux;
                  gtauy = dxdxi * tauy;
                  gridn[ix][iy] = Cross (gtauy, gtaux).Normalize();
                  }

                  for (int ix = 0; ix < hoplotn; ix++)
                  for (int iy = 0; iy < hoplotn; iy++)
                  {
                  glNormal3dv (gridn[ix][iy]);
                  glVertex3dv (grid[ix][iy]);

                  glNormal3dv (gridn[ix+1][iy]);
                  glVertex3dv (grid[ix+1][iy]);

                  glNormal3dv (gridn[ix+1][iy+1]);
                  glVertex3dv (grid[ix+1][iy+1]);

                  glNormal3dv (gridn[ix][iy+1]);
                  glVertex3dv (grid[ix][iy+1]);
                  }
                  }

                  glEnd ();
                */


              }
            else
              {



                Point<3> c(0,0,0);
                if (vispar.shrink < 1)
                  {
                    for (int j = 0; j < 5; j++)
                      {
                        Point<3> p = mesh->Point(el[j]);
                        c(0) += p(0) / 5;
                        c(1) += p(1) / 5;
                        c(2) += p(2) / 5;
                      }
                  }


                el.GetSurfaceTriangles (faces);

                if (el[0].IsValid())
                  {
                    glBegin (GL_TRIANGLES);

                    for (int j = 0; j < faces.Size(); j++)
                      {
                        ElementFace & face = faces[j];
                        Point<3> lp1 = mesh->Point (el.PNum(face[0]));
                        Point<3> lp2 = mesh->Point (el.PNum(face[1]));
                        Point<3> lp3 = mesh->Point (el.PNum(face[2]));
                        Vec<3> n = Cross (Vec<3> (lp1, lp2), Vec<3> (lp1, lp3));
                        n /= (n.Length()+1e-12);
                        n *= -1;
                        glNormal3d (n(0), n(1), n(2));

                        if (vispar.shrink < 1)
                          {
                            lp1 = c + vispar.shrink * (lp1 - c);
                            lp2 = c + vispar.shrink * (lp2 - c);
                            lp3 = c + vispar.shrink * (lp3 - c);
                          }

                        glVertex3d (lp1(0), lp1(1), lp1(2));
                        glVertex3d (lp2(0), lp2(1), lp2(2));
                        glVertex3d (lp3(0), lp3(1), lp3(2));
                      }

                    glEnd();
                  }
              }
          }
      }
    glEndList ();
  }

  void VisualSceneMesh :: BuildBadelList()
  {
    ;
  }

  void VisualSceneMesh :: BuildIdentifiedList()
  {
    ;
  }

  void VisualSceneMesh :: BuildDomainSurfList()
  {
    shared_ptr<Mesh> mesh = GetMesh();
    
    if (domainsurflist)
      glDeleteLists (domainsurflist, 1);

    domainsurflist = glGenLists (1);
    glNewList (domainsurflist, GL_COMPILE);

    int j;
    glLineWidth (1.0f);

    glDisable (GL_COLOR_MATERIAL);

    for (auto sel : mesh->SurfaceElements())
      {
        Element2d el (sel);

        int drawel = 1;
        for (j = 1; j <= el.GetNP(); j++)
          {
            if (!el.PNum(j).IsValid())
              drawel = 0;
          }

        if (!drawel)
          continue;

        if (!mesh->HasFaceDescriptor(el))
          continue;
        int domin = mesh->GetFaceDescriptor(el.GetIndex()).DomainIn();
        int domout = mesh->GetFaceDescriptor(el.GetIndex()).DomainOut();

        int fac;
        if (domin == vispar.drawdomainsurf)
          fac = 1;
        else if (domout == vispar.drawdomainsurf)
          fac = -1;
        else
          continue;


        GLfloat matcol[] = { 1, 0, 0, 1 };
        glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, matcol);


        if (el.GetNP() == 3)
          {
            glBegin (GL_TRIANGLES);

            const Point<3> & lp1 = mesh->Point (el[0]);
            const Point<3> & lp2 = mesh->Point (el[1]);
            const Point<3> & lp3 = mesh->Point (el[2]);
            Vec<3> n = Cross (Vec<3> (lp1, lp2), Vec<3> (lp1, lp3));
            n /= ( fac * (n.Length()+1e-12));
            glNormal3d (n(0), n(1), n(2));

            if (!vispar.colormeshsize)
              {
                glVertex3d (lp1(0), lp1(1), lp1(2));
                glVertex3d (lp2(0), lp2(1), lp2(2));
                glVertex3d (lp3(0), lp3(1), lp3(2));
              }
            glEnd();
          }
        else if (el.GetNP() == 4)
          {
            glBegin (GL_QUADS);

            const Point<3> & lp1 = mesh->Point (el[0]);
            const Point<3> & lp2 = mesh->Point (el[1]);
            const Point<3> & lp3 = mesh->Point (el[3]);
            const Point<3> & lp4 = mesh->Point (el[2]);
            Vec<3> n = Cross (Vec<3> (lp1, lp2),
                             Vec<3> (lp1, Center (lp3, lp4)));
            n /= (fac * (n.Length()+1e-12));
            glNormal3d (n(0), n(1), n(2));
            glVertex3d (lp1(0), lp1(1), lp1(2));
            glVertex3d (lp2(0), lp2(1), lp2(2));
            glVertex3d (lp4(0), lp4(1), lp4(2));
            glVertex3d (lp3(0), lp3(1), lp3(2));
            glEnd();
          }
        else if (el.GetNP() == 6)
          {
            glBegin (GL_TRIANGLES);
            static int trigs[4][3] = {
              { 1, 6, 5 },
              { 2, 4, 6 },
              { 3, 5, 4 },
              { 4, 5, 6 } };

            for (j = 0; j < 4; j++)
              {
                const Point<3> & lp1 = mesh->Point (el.PNum(trigs[j][0]));
                const Point<3> & lp2 = mesh->Point (el.PNum(trigs[j][1]));
                const Point<3> & lp3 = mesh->Point (el.PNum(trigs[j][2]));
                Vec<3> n = Cross (Vec<3> (lp1, lp2), Vec<3> (lp1, lp3));
                n /= (fac * (n.Length() + 1e-12));
                glNormal3d (n(0), n(1), n(2));
                glVertex3d (lp1(0), lp1(1), lp1(2));
                glVertex3d (lp2(0), lp2(1), lp2(2));
                glVertex3d (lp3(0), lp3(1), lp3(2));
              }
            glEnd();
          }
      }
    glEndList ();
  }





  bool VisualSelect :: SelectSurfaceElement (shared_ptr<Mesh> mesh, int px, int py, Point<3> &p, bool select_on_clipping_plane)
  {
    selelement = -1;
    // marker = nullopt;
    if(px != x || py != y)
    {
      x = px;
      y = py;
    }

    glGetIntegerv (GL_VIEWPORT, viewport);
    // GLenum err;
    if(framebuffer == 0 || viewport[2] != width || viewport[3] != height)
    {
      width = viewport[2];
      height = viewport[3];
      if(framebuffer != 0)
      {
        glDeleteRenderbuffers(2, render_buffers);
        glDeleteFramebuffers(1, &framebuffer);
      }

      glGenFramebuffers(1, &framebuffer);
      glBindFramebuffer(GL_FRAMEBUFFER, framebuffer);

      // create, reserve and attach color and depth renderbuffer
      glGenRenderbuffers(2, render_buffers);
      glBindRenderbuffer(GL_RENDERBUFFER, render_buffers[0]);
      glRenderbufferStorage(GL_RENDERBUFFER, GL_RGB16, width, height);
      glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, render_buffers[0]);

      glBindRenderbuffer(GL_RENDERBUFFER, render_buffers[1]);
      glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, width, height);
      glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, render_buffers[1]);

      // check if framebuffer status is complete
      if(int fbstatus; (fbstatus = glCheckFramebufferStatus(GL_FRAMEBUFFER)) != GL_FRAMEBUFFER_COMPLETE)
        cerr << "no frame buffer " << fbstatus << endl;

    }
      glFlush();

      glBindFramebuffer(GL_FRAMEBUFFER, framebuffer);

      glEnable(GL_DEPTH_TEST);
      glClearColor(0, 0, 0, 1.0);
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

      glMatrixMode (GL_MODELVIEW);
      glPushMatrix();
      glLoadIdentity();
      glMultMatrixd (transformationmat);

      glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
      auto hy = viewport[3] - py;

      if (enable_clipping_plane)
      {
        glClipPlane(GL_CLIP_PLANE0, clipplane);
        glEnable(GL_CLIP_PLANE0);
        Vec<3> n(clipplane[0], clipplane[1], clipplane[2]);
        double len = Abs(n);
        double mu = -clipplane[3] / (len*len);
        Point<3> p (mu * n);
        n /= len;
        Vec<3> t1 = GetNormal (n);
        Vec<3> t2 = Cross (n, t1);

        double xi1mid = (center - p) * t1;
        double xi2mid = (center - p) * t2;

        if(select_on_clipping_plane)
        {
          glColor3us(0,0,0);
          glBegin (GL_QUADS);
          glVertex3dv (p + (xi1mid-rad) * t1 + (xi2mid-rad) * t2);
          glVertex3dv (p + (xi1mid+rad) * t1 + (xi2mid-rad) * t2);
          glVertex3dv (p + (xi1mid+rad) * t1 + (xi2mid+rad) * t2);
          glVertex3dv (p + (xi1mid-rad) * t1 + (xi2mid+rad) * t2);
          glEnd ();
        }
      }
      glCallList (list);
      glFinish();

      glGetDoublev (GL_PROJECTION_MATRIX, projmat);
      auto found = Unproject(px, py, p);
      if(found)
      {
        // marker = p;
        GLushort numbers[3];
        glReadPixels (px, hy, 1, 1, GL_RGB, GL_UNSIGNED_SHORT, numbers);
        selelement = numbers[0] + numbers[1]*(1<<16);
      }
      glBindFramebuffer(GL_FRAMEBUFFER, 0);
      glPopMatrix();

    return found;
  }

  bool VisualSceneMesh :: Unproject(int px, int py, Point<3> &p)
  {
    return select.Unproject(px, py, p);
  }

  ngcore::IVec<2> VisualSceneMesh :: Project(Point<3> p)
  {
    Point<3> pwin;
    gluProject(p[0], p[1], p[2], transformationmat, select.projmat, select.viewport,
        &pwin[0], &pwin[1], &pwin[2]);

    return ngcore::IVec<2>(pwin[0]+0.5, select.viewport[3]-pwin[1]+0.5);
  }


  bool VisualSceneMesh :: SelectSurfaceElement(int px, int py, Point<3> &p, bool select_on_clipping_plane) {
    BuildFilledList(true);
    memcpy(select.transformationmat, transformationmat, sizeof(transformationmat));
    memcpy(select.clipplane, clipplane, sizeof(clipplane));
    select.center = center;
    select.rad = rad;
    select.enable_clipping_plane = vispar.clipping.enable;
    bool found = select.SelectSurfaceElement(GetMesh(), px, py, p, select_on_clipping_plane);
    selelement = select.selelement;
    return found;
  }

  void VisualSceneMesh :: MouseDblClick (int px, int py)
  {
    Point<3> p;
    bool found_point = SelectSurfaceElement(px, py, p, false);

    if(selelement>0)
      {
        const Element2dRef & sel = (*GetMesh())[SurfaceElementIndex::FromNr1(selelement)];
        SetSelectedFace(sel.GetIndex().Nr1());

        auto pi_nearest = sel[0];
        double min_dist = 1e99;
        for(auto pi : sel.PNums())
          if(Dist2(GetMesh()->Point(pi), p) < min_dist)
          {
            min_dist = Dist2(GetMesh()->Point(pi), p);
            pi_nearest = pi;
          }
        auto p_win = Project(GetMesh()->Point(pi_nearest));
        if(abs(p_win[0]-px) < 5 && abs(p_win[1]-py) < 5)
        {
          marker = GetMesh()->Point(pi_nearest);
          selpoint = pi_nearest;
          cout << "select point " << pi_nearest << " at " << *marker << endl;
        }
        else
        {
            marker = p;
            cout << endl << "select element " << selelement
              << " on face " << sel.GetIndex();
            // output face name
            auto mesh = GetMesh();
            string name;
            if(mesh->GetDimension() == 3)
              name = mesh->GetFaceDescriptor(sel.GetIndex()).GetBCName();
            else
              name = mesh->GetMaterial(sel.GetIndex().Nr1());

            if(name != "")
              cout << " with name " << name;
            cout << endl;
            if(mesh->GetDimension() == 3) {
              auto & fd = mesh->GetFaceDescriptor(sel.GetIndex());
              auto domin = fd.DomainIn();
              auto domout = fd.DomainOut();
              string name_in = domin >0 ? mesh->GetMaterial(domin) : "";
              string name_out = domout >0 ? mesh->GetMaterial(domout) : "";
              cout << "\tadjacent domains " << domin << ": " << name_in << ", " << domout << ": " << name_out << endl;
            }
            cout << "\tpoint: " << p << endl;;
            cout << "\tnodes: ";
            for (int i = 0; i < sel.GetNP(); i++)
              cout << sel[i] << " ";
            cout << endl;
        }
      }

    if(found_point && user_me_handler)
    {
      if (selelement != -1)
        user_me_handler -> DblClick (selelement-1, p[0], p[1], p[2]);
    }

    lock = {};
  }



  void VisualSceneMesh :: SetSelectedFace (int asf)
  {
    if(selface != asf)
    {
      selface = asf;
      BuildColorTexture();
    }
  }


  static bool dummy_init_var = [] () {
    on_set_global_mesh = [](shared_ptr<Mesh> mesh) {
      vsmesh.SetMesh(mesh);
    };
    return true;
  }();
}




#ifdef NG_PYTHON
#include <../general/ngpython.hpp>
#include "../include/nginterface.h"

NGGUI_API void ExportMeshVis(py::module &m)
{
  using namespace netgen;
  vispar.drawcolorbar = true;
  vispar.drawnetgenlogo = true;
  vispar.drawcoordinatecross = true;
  vispar.drawfilledtrigs = true;
  vispar.drawdomainsurf = true;
  vispar.drawhexes = true;
  vispar.drawtets = true;
  vispar.drawprisms = true;
  vispar.drawoutline = true;
  py::class_<VisualSceneMesh, shared_ptr<VisualSceneMesh>>
    (m, "VisualSceneMesh")
    .def("Draw", &VisualSceneMesh::DrawScene)
    ;

  m.def("VS", FunctionPointer
          ([](shared_ptr<Mesh> mesh)
           {
             auto vs = make_shared<VisualSceneMesh>();
             // vs->SetMesh(mesh);
             SetGlobalMesh (mesh);
             return vs;
           }));

  m.def("MouseMove", FunctionPointer
          ([](VisualSceneMesh &vsmesh, int oldx, int oldy, int newx, int 
              newy, char mode)
           {
             vsmesh.MouseMove(oldx, oldy, newx, newy, mode);
           }));
  m.def("SelectFace", FunctionPointer
      ([] (int facenr) {
       vsmesh.SetSelectedFace(facenr);
       }));
  m.def("GetGlobalMesh", FunctionPointer
      ([] () {
       return vsmesh.GetMesh();
       }));
}
// BOOST_PYTHON_MODULE(libvisual)
// {
//   ExportMeshVis();
// }
#endif
