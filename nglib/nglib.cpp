/**************************************************************************/
/* File:   nglib.cpp                                                      */
/* Author: Joachim Schoeberl                                              */
/* Date:   7. May. 2000                                                   */
/**************************************************************************/

/*
  
  Interface to the netgen meshing kernel
  
*/
#include <mystdlib.h>
#include <myadt.hpp>

#include <linalg.hpp>
#include <csg.hpp>
#include <stlgeom.hpp>
#include <geometry2d.hpp>
#include <meshing.hpp>
#include <../meshing/soldata.hpp>

#include <nginterface.h>


namespace netgen {
   extern void MeshFromSpline2D (SplineGeometry2d & geometry,
                                 shared_ptr<Mesh> & mesh, 
                                 MeshingParameters & mp);
   extern MeshingParameters mparam;
   DLL_HEADER extern STLParameters stlparam;
}



/*
namespace netgen
{
  int id = 0, ntasks = 1;
}
*/


/*
// should not be needed (occ currently requires it)
namespace netgen {
#include "../libsrc/visualization/vispar.hpp"
  VisualizationParameters vispar;
  VisualizationParameters :: VisualizationParameters() { ; }
}
*/


namespace nglib {
#include "nglib.h"
}

using namespace netgen;

// constants and types:

namespace nglib
{
  inline void NOOP_Deleter(void *) { ; }

  
   // initialize, deconstruct Netgen library:
   NGLIB_API void Ng_Init ()
   {
      mycout = &cout;
      myerr = &cerr;
      // netgen::testout->SetOutStream (new ofstream ("test.out"));
      // testout = new ofstream ("test.out");
   }




   // Clean-up functions before ending usage of nglib
   NGLIB_API void Ng_Exit ()
   {
      ;
   }




   // Create a new netgen mesh object
   NGLIB_API Ng_Mesh * Ng_NewMesh ()
   {
      Mesh * mesh = new Mesh;  
      mesh->AddFaceDescriptor (FaceDescriptor (1, 1, 0, 1));
      return (Ng_Mesh*)(void*)mesh;
   }




   // Delete an existing netgen mesh object
   NGLIB_API void Ng_DeleteMesh (Ng_Mesh * mesh)
   {
      if(mesh != NULL)
      {
         // Delete the Mesh structures
         ((Mesh*)mesh)->DeleteMesh();

         // Now delete the Mesh class itself
         delete (Mesh*)mesh;

         // Set the Ng_Mesh pointer to NULL
         mesh = NULL;
      }
   }




   // Save a netgen mesh in the native VOL format 
   NGLIB_API void Ng_SaveMesh(Ng_Mesh * mesh, const char* filename)
   {
      ((Mesh*)mesh)->Save(filename);
   }




   // Load a netgen native VOL mesh from a given file
   NGLIB_API Ng_Mesh * Ng_LoadMesh(const char* filename)
   {
      Mesh * mesh = new Mesh;
      mesh->Load(filename);
      return ( (Ng_Mesh*)mesh );
   }




   // Merge another mesh file into the currently loaded one
   NGLIB_API Ng_Result Ng_MergeMesh( Ng_Mesh* mesh, const char* filename)
   {
      Ng_Result status = NG_OK;

      ifstream infile(filename);
      Mesh * m = (Mesh*)mesh;

      if(!infile.good())
      {
         status = NG_FILE_NOT_FOUND;
      }

      if(!m)
      {
         status = NG_ERROR;
      }

      if(status == NG_OK)
      {
         const int num_pts = m->GetNP();
         const int face_offset = m->GetNFD();

         m->Merge(infile, face_offset);

         if(m->GetNP() > num_pts)
         {
            status = NG_OK;
         }
         else
         {
            status = NG_ERROR;
         }
      }

      return status;
   }




   // Merge another mesh file into the currently loaded one
   NGLIB_API Ng_Result Ng_MergeMesh( Ng_Mesh* mesh1, Ng_Mesh* mesh2)
   {
      return NG_ERROR;
   }




   // Manually add a point to an existing mesh object
   NGLIB_API void Ng_AddPoint (Ng_Mesh * mesh, double * x)
   {
      Mesh * m = (Mesh*)mesh;
      m->AddPoint (Point3d (x[0], x[1], x[2]));
   }




   // Manually add a surface element of a given type to an existing mesh object
   NGLIB_API void Ng_AddSurfaceElement (Ng_Mesh * mesh, Ng_Surface_Element_Type et,
                                         int * pi)
   {
      Mesh * m = (Mesh*)mesh;
      Element2d el (3);
      el.SetIndex (1);
      el.PNum(1) = pi[0];
      el.PNum(2) = pi[1];
      el.PNum(3) = pi[2];
      m->AddSurfaceElement (el);
   }




   // Manually add a volume element of a given type to an existing mesh object
   NGLIB_API void Ng_AddVolumeElement (Ng_Mesh * mesh, Ng_Volume_Element_Type et,
                                        int * pi)
   {
      Mesh * m = (Mesh*)mesh;
      Element el (4);
      el.SetIndex (1);
      el.PNum(1) = pi[0];
      el.PNum(2) = pi[1];
      el.PNum(3) = pi[2];
      el.PNum(4) = pi[3];
      m->AddVolumeElement (el);
   }




   // Obtain the number of points in the mesh
   NGLIB_API int Ng_GetNP (Ng_Mesh * mesh)
   {
      return ((Mesh*)mesh) -> GetNP();
   }




   // Obtain the number of surface elements in the mesh
   NGLIB_API int Ng_GetNSE (Ng_Mesh * mesh)
   {
      return ((Mesh*)mesh) -> GetNSE();
   }




   // Obtain the number of volume elements in the mesh
   NGLIB_API int Ng_GetNE (Ng_Mesh * mesh)
   {
      return ((Mesh*)mesh) -> GetNE();
   }




   //  Return point coordinates of a given point index in the mesh
   NGLIB_API void Ng_GetPoint (Ng_Mesh * mesh, int num, double * x)
   {
      const Point3d & p = ((Mesh*)mesh)->Point(num);
      x[0] = p.X();
      x[1] = p.Y();
      x[2] = p.Z();
   }




   // Return the surface element at a given index "pi"
   NGLIB_API Ng_Surface_Element_Type 
      Ng_GetSurfaceElement (Ng_Mesh * mesh, int num, int * pi)
   {
      const Element2d & el = ((Mesh*)mesh)->SurfaceElement(num);
      for (int i = 1; i <= el.GetNP(); i++)
         pi[i-1] = el.PNum(i);
      Ng_Surface_Element_Type et;
      switch (el.GetNP())
      {
      case 3: et = NG_TRIG; break;
      case 4: et = NG_QUAD; break;
      case 6: 
         switch (el.GetNV())
         {
         case 3: et = NG_TRIG6; break;
         case 4: et = NG_QUAD6; break;
         default:
            et = NG_TRIG6; break;
         }
         break;
      case 8: et = NG_QUAD8; break;
      default:
         et = NG_TRIG; break; // for the compiler
      }
      return et;
   }




   // Return the volume element at a given index "pi"
   NGLIB_API Ng_Volume_Element_Type
      Ng_GetVolumeElement (Ng_Mesh * mesh, int num, int * pi)
   {
      const Element & el = ((Mesh*)mesh)->VolumeElement(num);
      for (int i = 1; i <= el.GetNP(); i++)
         pi[i-1] = el.PNum(i);
      Ng_Volume_Element_Type et;
      switch (el.GetNP())
      {
      case 4: et = NG_TET; break;
      case 5: et = NG_PYRAMID; break;
      case 6: et = NG_PRISM; break;
      case 10: et = NG_TET10; break;
      default:
         et = NG_TET; break; // for the compiler
      }
      return et;
   }




   // Set a global limit on the maximum mesh size allowed
   NGLIB_API void Ng_RestrictMeshSizeGlobal (Ng_Mesh * mesh, double h)
   {
      ((Mesh*)mesh) -> SetGlobalH (h);
   }




   // Set a local limit on the maximum mesh size allowed around the given point
   NGLIB_API void Ng_RestrictMeshSizePoint (Ng_Mesh * mesh, double * p, double h)
   {
      ((Mesh*)mesh) -> RestrictLocalH (Point3d (p[0], p[1], p[2]), h);
   }




   // Set a local limit on the maximum mesh size allowed within a given box region
   NGLIB_API void Ng_RestrictMeshSizeBox (Ng_Mesh * mesh, double * pmin, double * pmax, double h)
   {
      for (double x = pmin[0]; x < pmax[0]; x += h)
         for (double y = pmin[1]; y < pmax[1]; y += h)
            for (double z = pmin[2]; z < pmax[2]; z += h)
               ((Mesh*)mesh) -> RestrictLocalH (Point3d (x, y, z), h);
   }




   // Generates volume mesh from an existing surface mesh
   NGLIB_API Ng_Result Ng_GenerateVolumeMesh (Ng_Mesh * mesh, Ng_Meshing_Parameters * mp)
   {
      Mesh * m = (Mesh*)mesh;

      // Philippose - 30/08/2009
      // Do not locally re-define "mparam" here... "mparam" is a global 
      // object 
      //MeshingParameters mparam;
      mp->Transfer_Parameters();

      m->CalcLocalH(mparam.grading);

      MeshVolume (mparam, *m);
      RemoveIllegalElements (*m);
      OptimizeVolume (mparam, *m);

      return NG_OK;
   }




   /* ------------------ 2D Meshing Functions ------------------------- */
   NGLIB_API void Ng_AddPoint_2D (Ng_Mesh * mesh, double * x)
   {
      Mesh * m = (Mesh*)mesh;

      m->AddPoint (Point3d (x[0], x[1], 0));
   }




   NGLIB_API void Ng_AddBoundarySeg_2D (Ng_Mesh * mesh, int pi1, int pi2)
   {
      Mesh * m = (Mesh*)mesh;

      Segment seg;
      seg[0] = pi1;
      seg[1] = pi2;
      m->AddSegment (seg);
   }




   NGLIB_API int Ng_GetNP_2D (Ng_Mesh * mesh)
   {
      Mesh * m = (Mesh*)mesh;
      return m->GetNP();
   }




   NGLIB_API int Ng_GetNE_2D (Ng_Mesh * mesh)
   {
      Mesh * m = (Mesh*)mesh;
      return m->GetNSE();
   }




   NGLIB_API int Ng_GetNSeg_2D (Ng_Mesh * mesh)
   {
      Mesh * m = (Mesh*)mesh;
      return m->GetNSeg();
   }




   NGLIB_API void Ng_GetPoint_2D (Ng_Mesh * mesh, int num, double * x)
   {
      Mesh * m = (Mesh*)mesh;

      Point<3> & p = m->Point(num);
      x[0] = p(0);
      x[1] = p(1);
   }




   NGLIB_API Ng_Surface_Element_Type
      Ng_GetElement_2D (Ng_Mesh * mesh, int num, int * pi, int * matnum)
   {
      const Element2d & el = ((Mesh*)mesh)->SurfaceElement(num);
      for (int i = 1; i <= el.GetNP(); i++)
         pi[i-1] = el.PNum(i);

      Ng_Surface_Element_Type et;
      switch (el.GetNP())
      {
      case 3: et = NG_TRIG; break;
      case 4: et = NG_QUAD; break;
      case 6: 
         switch (el.GetNV())
         {
         case 3: et = NG_TRIG6; break;
         case 4: et = NG_QUAD6; break;
         default:
            et = NG_TRIG6; break;
         }
         break;
      case 8: et = NG_QUAD8; break;
      default:
         et = NG_TRIG; break; // for the compiler
      }

      if (matnum)
         *matnum = el.GetIndex();

      return et;
   }




   NGLIB_API void Ng_GetSegment_2D (Ng_Mesh * mesh, int num, int * pi, int * matnum)
   {
      const Segment & seg = ((Mesh*)mesh)->LineSegment(num);
      pi[0] = seg[0];
      pi[1] = seg[1];

      if (matnum)
         *matnum = seg.edgenr;
   }




   NGLIB_API Ng_Geometry_2D * Ng_LoadGeometry_2D (const char * filename)
   {
      SplineGeometry2d * geom = new SplineGeometry2d();
      geom -> Load (filename);
      return (Ng_Geometry_2D *)geom;
   }


   NGLIB_API Ng_Result Ng_GenerateMesh_2D (Ng_Geometry_2D * geom,
                                            Ng_Mesh ** mesh,
                                            Ng_Meshing_Parameters * mp)
   {
      // use global variable mparam
      //  MeshingParameters mparam;  
      mp->Transfer_Parameters();

      shared_ptr<Mesh> m(new Mesh, &NOOP_Deleter);
      MeshFromSpline2D (*(SplineGeometry2d*)geom, m, mparam);
      // new shared_ptr<Mesh> (m);  // hack to keep mesh m alive 
      
      cout << m->GetNSE() << " elements, " << m->GetNP() << " points" << endl;

      *mesh = (Ng_Mesh*)m.get();
      return NG_OK;
   }




   NGLIB_API void Ng_HP_Refinement (Ng_Geometry_2D * geom,
      Ng_Mesh * mesh,
      int levels)
   {
      Refinement ref(*(SplineGeometry2d*)geom);
      HPRefinement (*(Mesh*)mesh, &ref, levels);
   }




   NGLIB_API void Ng_HP_Refinement (Ng_Geometry_2D * geom,
      Ng_Mesh * mesh,
      int levels, double parameter)
   {
      Refinement ref(*(SplineGeometry2d*)geom);
      HPRefinement (*(Mesh*)mesh, &ref, levels, parameter);
   }




   NgArray<STLReadTriangle> readtrias; //only before initstlgeometry
   NgArray<Point<3> > readedges; //only before init stlgeometry

   // loads geometry from STL file
   NGLIB_API Ng_STL_Geometry * Ng_STL_LoadGeometry (const char * filename, int binary)
   {
      int i;
      STLGeometry geom;
      STLGeometry* geo;
      ifstream ist(filename);

      if (binary)
      {
         geo = geom.LoadBinary(ist);
      }
      else
      {
         geo = geom.Load(ist);
      }

      readtrias.SetSize(0);
      readedges.SetSize(0);

      Point3d p;
      Vec3d normal;
      double p1[3];
      double p2[3];
      double p3[3];
      double n[3];

      Ng_STL_Geometry * geo2 = Ng_STL_NewGeometry();

      for (i = 1; i <= geo->GetNT(); i++)
      {
         const STLTriangle& t = geo->GetTriangle(i);
         p = geo->GetPoint(t.PNum(1));
         p1[0] = p.X(); p1[1] = p.Y(); p1[2] = p.Z(); 
         p = geo->GetPoint(t.PNum(2));
         p2[0] = p.X(); p2[1] = p.Y(); p2[2] = p.Z(); 
         p = geo->GetPoint(t.PNum(3));
         p3[0] = p.X(); p3[1] = p.Y(); p3[2] = p.Z();
         normal = t.Normal();
         n[0] = normal.X(); n[1] = normal.Y(); n[2] = normal.Z();

         Ng_STL_AddTriangle(geo2, p1, p2, p3, n);
      }

      return geo2;
   }




   // generate new STL Geometry
   NGLIB_API Ng_STL_Geometry * Ng_STL_NewGeometry ()
   {
      return (Ng_STL_Geometry*)(void*)new STLGeometry;
   } 




   // after adding triangles (and edges) initialize
   NGLIB_API Ng_Result Ng_STL_InitSTLGeometry (Ng_STL_Geometry * geom)
   {
      STLGeometry* geo = (STLGeometry*)geom;
      geo->InitSTLGeometry(readtrias);
      readtrias.SetSize(0);

      if (readedges.Size() != 0)
      {
         /*
         for (int i = 1; i <= readedges.Size(); i+=2)
         {
         cout << "e(" << readedges.Get(i) << "," << readedges.Get(i+1) << ")" << endl;
         }
         */
         geo->AddEdges(readedges);
	 readedges.SetSize(0);
      }

      if (geo->GetStatus() == STLTopology::STL_GOOD || geo->GetStatus() == STLTopology::STL_WARNING) return NG_OK;
      return NG_SURFACE_INPUT_ERROR;
   }




   // automatically generates edges:
   NGLIB_API Ng_Result Ng_STL_MakeEdges (Ng_STL_Geometry * geom,
                                          Ng_Mesh* mesh,
                                          Ng_Meshing_Parameters * mp)
   {
      STLGeometry* stlgeometry = (STLGeometry*)geom;
      Mesh* me = (Mesh*)mesh;
      me->SetGeometry( shared_ptr<NetgenGeometry>(stlgeometry, &NOOP_Deleter) );

      // Philippose - 27/07/2009
      // Do not locally re-define "mparam" here... "mparam" is a global 
      // object 
      //MeshingParameters mparam;
      mp->Transfer_Parameters();

      me -> SetGlobalH (mparam.maxh);
      me -> SetLocalH (stlgeometry->GetBoundingBox().PMin() - Vec3d(10, 10, 10),
                       stlgeometry->GetBoundingBox().PMax() + Vec3d(10, 10, 10),
                       0.3);

      // cout << "meshsize = " << mp->meshsize_filename << endl;
      if (mp->meshsize_filename)
        me -> LoadLocalMeshSize (mp->meshsize_filename);

      /*
      if (mp->meshsize_filename)
        {
          ifstream infile (mp->meshsize_filename);
          if (!infile.good()) return NG_FILE_NOT_FOUND;
          me -> LoadLocalMeshSize (infile);
        }
      */

      STLMeshing (*stlgeometry, *me, mparam, stlparam);

      stlgeometry->edgesfound = 1;
      stlgeometry->surfacemeshed = 0;
      stlgeometry->surfaceoptimized = 0;
      stlgeometry->volumemeshed = 0;

      return NG_OK;
   }




   // generates mesh, empty mesh be already created.
   NGLIB_API Ng_Result Ng_STL_GenerateSurfaceMesh (Ng_STL_Geometry * geom,
                                                    Ng_Mesh* mesh,
                                                    Ng_Meshing_Parameters * mp)
   {
      STLGeometry* stlgeometry = (STLGeometry*)geom;
      Mesh* me = (Mesh*)mesh;
      me->SetGeometry( shared_ptr<NetgenGeometry>(stlgeometry, &NOOP_Deleter) );

      // Philippose - 27/07/2009
      // Do not locally re-define "mparam" here... "mparam" is a global 
      // object
      //MeshingParameters mparam;
      mp->Transfer_Parameters();


      /*
      me -> SetGlobalH (mparam.maxh);
      me -> SetLocalH (stlgeometry->GetBoundingBox().PMin() - Vec3d(10, 10, 10),
      stlgeometry->GetBoundingBox().PMax() + Vec3d(10, 10, 10),
      0.3);
      */
      /*
      STLMeshing (*stlgeometry, *me);

      stlgeometry->edgesfound = 1;
      stlgeometry->surfacemeshed = 0;
      stlgeometry->surfaceoptimized = 0;
      stlgeometry->volumemeshed = 0;
      */  
      int retval = STLSurfaceMeshing (*stlgeometry, *me, mparam, stlparam);
      if (retval == MESHING3_OK)
      {
         (*mycout) << "Success !!!!" << endl;
         stlgeometry->surfacemeshed = 1;
         stlgeometry->surfaceoptimized = 0;
         stlgeometry->volumemeshed = 0;
      } 
      else if (retval == MESHING3_OUTERSTEPSEXCEEDED)
      {
         (*mycout) << "ERROR: Give up because of too many trials. Meshing aborted!" << endl;
      }
      else if (retval == MESHING3_TERMINATE)
      {
         (*mycout) << "Meshing Stopped!" << endl;
      }
      else
      {
         (*mycout) << "ERROR: Surface meshing not successful. Meshing aborted!" << endl;
      }


      STLSurfaceOptimization (*stlgeometry, *me, mparam);

      return NG_OK;
   }




   // fills STL Geometry
   // positive orientation
   // normal vector may be null-pointer
   NGLIB_API void Ng_STL_AddTriangle (Ng_STL_Geometry * geom, 
                                       double * p1, double * p2, double * p3, 
                                       double * nv)
   {
      Point<3> apts[3];
      apts[0] = Point<3>(p1[0],p1[1],p1[2]);
      apts[1] = Point<3>(p2[0],p2[1],p2[2]);
      apts[2] = Point<3>(p3[0],p3[1],p3[2]);

      Vec<3> n;
      if (!nv)
         n = Cross (apts[0]-apts[1], apts[0]-apts[2]);
      else
         n = Vec<3>(nv[0],nv[1],nv[2]);

      readtrias.Append(STLReadTriangle(apts,n));
   }

   // add (optional) edges:
   NGLIB_API void Ng_STL_AddEdge (Ng_STL_Geometry * geom, 
      double * p1, double * p2)
   {
      readedges.Append(Point3d(p1[0],p1[1],p1[2]));
      readedges.Append(Point3d(p2[0],p2[1],p2[2]));
   }








   // ------------------ Begin - Meshing Parameters related functions ------------------
   // Constructor for the local nglib meshing parameters class
   NGLIB_API Ng_Meshing_Parameters :: Ng_Meshing_Parameters()
   {
      uselocalh = 1;

      maxh = 1000;
      minh = 0.0;

      fineness = 0.5;
      grading = 0.3;

      elementsperedge = 2.0;
      elementspercurve = 2.0;

      closeedgeenable = 0;
      closeedgefact = 2.0;

	  minedgelenenable = 0;
	  minedgelen = 1e-4;

      second_order = 0;
      quad_dominated = 0;

      meshsize_filename = 0;

      optsurfmeshenable = 1;
      optvolmeshenable = 1;

      optsteps_2d = 3;
      optsteps_3d = 3;

      invert_tets = 0;
      invert_trigs = 0;

      check_overlap = 1;
      check_overlapping_boundary = 1;
   }




   // Reset the local meshing parameters to the default values
   NGLIB_API void Ng_Meshing_Parameters :: Reset_Parameters()
   {
      uselocalh = 1;

      maxh = 1000;
      minh = 0;

      fineness = 0.5;
      grading = 0.3;

      elementsperedge = 2.0;
      elementspercurve = 2.0;

      closeedgeenable = 0;
      closeedgefact = 2.0;

  	  minedgelenenable = 0;
	  minedgelen = 1e-4;

      second_order = 0;
      quad_dominated = 0;

      meshsize_filename = 0;

      optsurfmeshenable = 1;
      optvolmeshenable = 1;

      optsteps_2d = 3;
      optsteps_3d = 3;

      invert_tets = 0;
      invert_trigs = 0;

      check_overlap = 1;
      check_overlapping_boundary = 1;
   }




   // 
   NGLIB_API void Ng_Meshing_Parameters :: Transfer_Parameters()
   {
      mparam.uselocalh = uselocalh;
      
      mparam.maxh = maxh;
      mparam.minh = minh;

      mparam.grading = grading;
      mparam.curvaturesafety = elementspercurve;
      mparam.segmentsperedge = elementsperedge;

      mparam.secondorder = second_order;
      mparam.quad = quad_dominated;

      if (meshsize_filename)
        mparam.meshsizefilename = meshsize_filename;
      else
        mparam.meshsizefilename = "";
      mparam.optsteps2d = optsteps_2d;
      mparam.optsteps3d = optsteps_3d;

      mparam.inverttets = invert_tets;
      mparam.inverttrigs = invert_trigs;

      mparam.checkoverlap = check_overlap;
      mparam.checkoverlappingboundary = check_overlapping_boundary;
   }
   // ------------------ End - Meshing Parameters related functions --------------------




   // ------------------ Begin - Second Order Mesh generation functions ----------------
   NGLIB_API void Ng_Generate_SecondOrder(Ng_Mesh * mesh)
   {
     Refinement ref(*((Mesh*) mesh)->GetGeometry());
      ref.MakeSecondOrder(*(Mesh*) mesh);
   }




   NGLIB_API void Ng_2D_Generate_SecondOrder(Ng_Geometry_2D * geom,
					  Ng_Mesh * mesh)
   {
      ( (SplineGeometry2d*)geom ) -> GetRefinement().MakeSecondOrder( * (Mesh*) mesh );
   }




   NGLIB_API void Ng_STL_Generate_SecondOrder(Ng_STL_Geometry * geom,
					   Ng_Mesh * mesh)
   {
      ((STLGeometry*)geom)->GetRefinement().MakeSecondOrder(*(Mesh*) mesh);
   }




   NGLIB_API void Ng_CSG_Generate_SecondOrder (Ng_CSG_Geometry * geom,
					   Ng_Mesh * mesh)
   {
      ((CSGeometry*)geom)->GetRefinement().MakeSecondOrder(*(Mesh*) mesh);
   }




   // ------------------ End - Second Order Mesh generation functions ------------------




   // ------------------ Begin - Uniform Mesh Refinement functions ---------------------
   NGLIB_API void Ng_Uniform_Refinement (Ng_Mesh * mesh)
   {
     Refinement ref(*((Mesh*)mesh)->GetGeometry());
     ref.Refine ( * (Mesh*) mesh );
   }




   NGLIB_API void Ng_2D_Uniform_Refinement (Ng_Geometry_2D * geom,
      Ng_Mesh * mesh)
   {
      ( (SplineGeometry2d*)geom ) -> GetRefinement().Refine ( * (Mesh*) mesh );
   }




   NGLIB_API void Ng_STL_Uniform_Refinement (Ng_STL_Geometry * geom,
      Ng_Mesh * mesh)
   {
      ( (STLGeometry*)geom ) -> GetRefinement().Refine ( * (Mesh*) mesh );
   }




   NGLIB_API void Ng_CSG_Uniform_Refinement (Ng_CSG_Geometry * geom,
      Ng_Mesh * mesh)
   {
      ( (CSGeometry*)geom ) -> GetRefinement().Refine ( * (Mesh*) mesh );
   }




   // ------------------ End - Uniform Mesh Refinement functions -----------------------
} // End of namespace nglib




// compatibility functions:
namespace netgen 
{
   char geomfilename[255];

   NGLIB_API void MyError2 (const char * ch)
   {
      cerr << ch;
   }




   //Destination for messages, errors, ...
   NGLIB_API void Ng_PrintDest2(const char * s)
   {
#ifdef PARALLEL
     int id = 0;
     NG_MPI_Comm_rank(NG_MPI_COMM_WORLD, &id);
     if (id != 0) return;
#endif
     (*mycout) << s << flush;
   }


  /*
   NGLIB_API double GetTime ()
   {
      return 0;
   }
  */

  /*
#ifndef WIN32
   void ResetTime ()
   {
      ;
   }
#endif
  */


   void MyBeep (int i)
   {
      ;
   }



  //void Render() { ; }

} // End of namespace netgen


/*

#ifndef WIN32
void Ng_Redraw () { ; }
void Ng_ClearSolutionData() { ; }
#endif
void Ng_SetSolutionData (Ng_SolutionData * soldata) 
{ 
  delete soldata->solclass;
}
void Ng_InitSolutionData (Ng_SolutionData * soldata) { ; }
*/

// Force linking libinterface to libnglib
#include <../interface/writeuser.hpp>
void MyDummyToForceLinkingLibInterface(Mesh &mesh, NetgenGeometry &geom)
{
  netgen::WriteUserFormat("", mesh, /* geom, */ "");
}

