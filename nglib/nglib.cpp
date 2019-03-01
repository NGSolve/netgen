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
#include <../visualization/soldata.hpp>
#include <../interface/writeuser.hpp>

#ifdef OCCGEOMETRY
#include <occgeom.hpp>
#endif


namespace netgen {
   extern void MeshFromSpline2D (SplineGeometry2d & geometry,
                                 shared_ptr<Mesh> & mesh, 
                                 MeshingParameters & mp);
}



#ifdef PARALLEL
#include <mpi.h>

namespace netgen
{
  // int id = 0, ntasks = 1;
  MPI_Comm mesh_comm;
}
#endif


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

   class NullStreambuf : public std::streambuf
   {
      char dummyBuffer[64];
   protected:
      virtual int overflow(int c)
      {
         setp(dummyBuffer, dummyBuffer + sizeof(dummyBuffer));
         return (c == traits_type::eof()) ? '\0' : c;
      }
   };

   // initialize, deconstruct Netgen library:
   DLL_HEADER void Ng_Init (bool cout_to_null, bool cerr_to_null, bool testout_to_null)
   {
      static ostream* null_stream = new ostream(new NullStreambuf);
      mycout  = cout_to_null ? null_stream : &cout;
      myerr   = cerr_to_null ? null_stream : &cerr;
      testout = testout_to_null ? null_stream :  new ofstream("test.out");
   }




   // Clean-up functions before ending usage of nglib
   DLL_HEADER void Ng_Exit ()
   {
      ;
   }




   // Create a new netgen mesh object
   DLL_HEADER Ng_Mesh * Ng_NewMesh ()
   {
      Mesh * mesh = new Mesh;  
      mesh->AddFaceDescriptor (FaceDescriptor (1, 1, 0, 1));
      return (Ng_Mesh*)(void*)mesh;
   }




   // Delete an existing netgen mesh object
   DLL_HEADER void Ng_DeleteMesh (Ng_Mesh * mesh)
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
   DLL_HEADER void Ng_SaveMesh(Ng_Mesh * mesh, const char* filename)
   {
      ((Mesh*)mesh)->Save(filename);
   }




   // Load a netgen native VOL mesh from a given file
   DLL_HEADER Ng_Mesh * Ng_LoadMesh(const char* filename)
   {
      Mesh * mesh = new Mesh;
      mesh->Load(filename);
      return ( (Ng_Mesh*)mesh );
   }



   DLL_HEADER void Ng_ExportMesh(Ng_Mesh * ng_mesh, Ng_Export_Formats format, const char* filename)
   {
      Mesh * mesh = (Mesh*)ng_mesh;
      switch (format)
      {
      case NG_GMSH:
         WriteUserFormat( "Gmsh Format", *mesh, filename ); break;
      case NG_GMSH2:
         WriteUserFormat( "Gmsh2 Format", *mesh, filename ); break;
      case NG_VTK:
         WriteUserFormat( "VTK Format", *mesh, filename ); break;
      case NG_FLUENT:
         WriteUserFormat( "Fluent Format", *mesh, filename ); break;
      case NG_ABAQUS:
         WriteUserFormat( "Abaqus Format", *mesh, filename ); break;
      }
   }



   // Merge another mesh file into the currently loaded one
   DLL_HEADER Ng_Result Ng_MergeMesh( Ng_Mesh* mesh, const char* filename)
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
   DLL_HEADER Ng_Result Ng_MergeMesh( Ng_Mesh* mesh1, Ng_Mesh* mesh2)
   {
      return NG_ERROR;
   }




   // Manually add a point to an existing mesh object
   DLL_HEADER void Ng_AddPoint (Ng_Mesh * mesh, double * x)
   {
      Mesh * m = (Mesh*)mesh;
      m->AddPoint (Point3d (x[0], x[1], x[2]));
   }


   // Manually lock a specific point
   DLL_HEADER void Ng_AddLockedPoint(Ng_Mesh * mesh, int pi)
   {
      Mesh * m = (Mesh*)mesh;
      m->AddLockedPoint(pi);
   }


   // Manually add a surface element of a given type to an existing mesh object
   DLL_HEADER void Ng_AddSurfaceElement (Ng_Mesh * mesh, Ng_Surface_Element_Type et,
                                         int * pi, int domain)
   {
      int n = 3;
      switch (et)
      {
      case NG_TRIG:
         n = 3; break;
      case NG_QUAD:
         n = 4; break;
      case NG_QUAD6:
         n = 6; break;
      case NG_TRIG6:
         n = 6; break;
      case NG_QUAD8:
         n = 8; break;
      default: break;
      }
      
      Mesh * m = (Mesh*)mesh;
      Element2d el (n);
      el.SetIndex (domain);
      for (int i=0; i<n; ++i)
         el.PNum(i+1) = pi[i];
      m->AddSurfaceElement (el);
   }




   // Manually add a volume element of a given type to an existing mesh object
   DLL_HEADER void Ng_AddVolumeElement (Ng_Mesh * mesh, Ng_Volume_Element_Type et,
                                        int * pi, int domain)
   {
      int n = 4;
      switch (et)
      {
      case NG_TET:
         n = 4; break;
      case NG_PYRAMID:
         n = 5; break;
      case NG_PRISM:
         n = 6; break;
      case NG_HEX:
         n = 8; break;
      case NG_TET10:
         n = 10; break;
      default: break;
      }
      
      Mesh * m = (Mesh*)mesh;
      Element el (n);
      el.SetIndex (domain);
      for (int i=0; i<n; ++i)
         el.PNum(i+1) = pi[i];
      m->AddVolumeElement (el);
   }


   // Obtain the number of points in the mesh
   DLL_HEADER int Ng_GetNP (Ng_Mesh * mesh)
   {
      return ((Mesh*)mesh) -> GetNP();
   }




   // Obtain the number of surface elements in the mesh
   DLL_HEADER int Ng_GetNSE (Ng_Mesh * mesh)
   {
      return ((Mesh*)mesh) -> GetNSE();
   }




   // Obtain the number of volume elements in the mesh
   DLL_HEADER int Ng_GetNE (Ng_Mesh * mesh)
   {
      return ((Mesh*)mesh) -> GetNE();
   }




   //  Return point coordinates of a given point index in the mesh
   DLL_HEADER void Ng_GetPoint (Ng_Mesh * mesh, int num, double * x)
   {
      const Point3d & p = ((Mesh*)mesh)->Point(num);
      x[0] = p.X();
      x[1] = p.Y();
      x[2] = p.Z();
   }




   // Return the surface element at a given index "pi"
   DLL_HEADER Ng_Surface_Element_Type 
      Ng_GetSurfaceElement (Ng_Mesh * mesh, int num, int * pi, int * domain)
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
      if (domain)
        *domain = el.GetIndex();
      return et;
   }




   // Return the volume element at a given index "pi"
   DLL_HEADER Ng_Volume_Element_Type
      Ng_GetVolumeElement (Ng_Mesh * mesh, int num, int * pi, int * domain)
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
      case 8: et = NG_HEX; break;
      case 10: et = NG_TET10; break;
      default:
         et = NG_TET; break; // for the compiler
      }
      if (domain)
        *domain = el.GetIndex();
      return et;
   }




   // Set a global limit on the maximum mesh size allowed
   DLL_HEADER void Ng_RestrictMeshSizeGlobal (Ng_Mesh * mesh, double h)
   {
      ((Mesh*)mesh) -> SetGlobalH (h);
   }




   // Set a local limit on the maximum mesh size allowed around the given point
   DLL_HEADER void Ng_RestrictMeshSizePoint (Ng_Mesh * mesh, double * p, double h)
   {
      ((Mesh*)mesh) -> RestrictLocalH (Point3d (p[0], p[1], p[2]), h);
   }




   // Set a local limit on the maximum mesh size allowed within a given box region
   DLL_HEADER void Ng_RestrictMeshSizeBox (Ng_Mesh * mesh, double * pmin, double * pmax, double h)
   {
      for (double x = pmin[0]; x < pmax[0]; x += h)
         for (double y = pmin[1]; y < pmax[1]; y += h)
            for (double z = pmin[2]; z < pmax[2]; z += h)
               ((Mesh*)mesh) -> RestrictLocalH (Point3d (x, y, z), h);
   }




   // Generates volume mesh from an existing surface mesh
   DLL_HEADER Ng_Result Ng_GenerateVolumeMesh (Ng_Mesh * mesh, Ng_Meshing_Parameters * mp)
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




   // Optimize existing mesh
   DLL_HEADER Ng_Result Ng_OptimizeVolume(Ng_Mesh *mesh, Ng_Meshing_Parameters *mp)
   {
      Mesh * m = (Mesh*)mesh;

      mp->Transfer_Parameters();

      m->CalcLocalH(mparam.grading);

      RemoveIllegalElements(*m);
      OptimizeVolume(mparam, *m);

      return NG_OK;
   }




   /* ------------------ 2D Meshing Functions ------------------------- */
   DLL_HEADER void Ng_AddPoint_2D (Ng_Mesh * mesh, double * x)
   {
      Mesh * m = (Mesh*)mesh;

      m->AddPoint (Point3d (x[0], x[1], 0));
   }




   DLL_HEADER void Ng_AddBoundarySeg_2D (Ng_Mesh * mesh, int pi1, int pi2, int domain_in, int domain_out)
   {
      Mesh * m = (Mesh*)mesh;

      Segment seg;
      seg[0] = pi1;
      seg[1] = pi2;
      seg.domin = domain_in;
      seg.domout = domain_out;
     m->AddSegment(seg);
   }




   DLL_HEADER int Ng_GetNP_2D (Ng_Mesh * mesh)
   {
      Mesh * m = (Mesh*)mesh;
      return m->GetNP();
   }




   DLL_HEADER int Ng_GetNE_2D (Ng_Mesh * mesh)
   {
      Mesh * m = (Mesh*)mesh;
      return m->GetNSE();
   }




   DLL_HEADER int Ng_GetNSeg_2D (Ng_Mesh * mesh)
   {
      Mesh * m = (Mesh*)mesh;
      return m->GetNSeg();
   }




   DLL_HEADER void Ng_GetPoint_2D (Ng_Mesh * mesh, int num, double * x)
   {
      Mesh * m = (Mesh*)mesh;

      Point<3> & p = m->Point(num);
      x[0] = p(0);
      x[1] = p(1);
   }




   DLL_HEADER Ng_Surface_Element_Type
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




   DLL_HEADER void Ng_GetSegment_2D (Ng_Mesh * mesh, int num, int * pi, int * matnum)
   {
      const Segment & seg = ((Mesh*)mesh)->LineSegment(num);
      pi[0] = seg[0];
      pi[1] = seg[1];

      if (matnum)
         *matnum = seg.edgenr;
   }




   DLL_HEADER Ng_Geometry_2D * Ng_LoadGeometry_2D (const char * filename)
   {
      SplineGeometry2d * geom = new SplineGeometry2d();
      geom -> Load (filename);
      return (Ng_Geometry_2D *)geom;
   }


   DLL_HEADER Ng_Result Ng_GenerateMesh_2D (Ng_Geometry_2D * geom,
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




   DLL_HEADER void Ng_HP_Refinement (Ng_Geometry_2D * geom,
      Ng_Mesh * mesh,
      int levels)
   {
      Refinement2d ref(*(SplineGeometry2d*)geom);
      HPRefinement (*(Mesh*)mesh, &ref, levels);
   }




   DLL_HEADER void Ng_HP_Refinement (Ng_Geometry_2D * geom,
      Ng_Mesh * mesh,
      int levels, double parameter)
   {
      Refinement2d ref(*(SplineGeometry2d*)geom);
      HPRefinement (*(Mesh*)mesh, &ref, levels, parameter);
   }




   Array<STLReadTriangle> readtrias; //only before initstlgeometry
   Array<Point<3> > readedges; //only before init stlgeometry

   // loads geometry from STL file
   DLL_HEADER Ng_STL_Geometry * Ng_STL_LoadGeometry (const char * filename, int binary)
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
   DLL_HEADER Ng_STL_Geometry * Ng_STL_NewGeometry ()
   {
      return (Ng_STL_Geometry*)(void*)new STLGeometry;
   } 


   DLL_HEADER void Ng_STL_DeleteGeometry (Ng_STL_Geometry * geom)
   {
      if (geom)
      {
         STLGeometry* geometry = (STLGeometry*)geom;
         geometry->Clear();
         delete geometry;
         geometry = NULL;
      }
   }


   // after adding triangles (and edges) initialize
   DLL_HEADER Ng_Result Ng_STL_InitSTLGeometry (Ng_STL_Geometry * geom)
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
      }

      if (geo->GetStatus() == STLTopology::STL_GOOD || geo->GetStatus() == STLTopology::STL_WARNING) return NG_OK;
      return NG_SURFACE_INPUT_ERROR;
   }




   // automatically generates edges:
   DLL_HEADER Ng_Result Ng_STL_MakeEdges (Ng_STL_Geometry * geom,
                                          Ng_Mesh* mesh,
                                          Ng_Meshing_Parameters * mp,
                                          Ng_STL_Parameters * stlp)
   {
      STLGeometry* stlgeometry = (STLGeometry*)geom;
      Mesh* me = (Mesh*)mesh;

      // Philippose - 27/07/2009
      // Do not locally re-define "mparam" here... "mparam" is a global 
      // object 
      //MeshingParameters mparam;
      mp->Transfer_Parameters();
      if (stlp) stlp->Transfer_Parameters();

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

      STLMeshing (*stlgeometry, *me);

      stlgeometry->edgesfound = 1;
      stlgeometry->surfacemeshed = 0;
      stlgeometry->surfaceoptimized = 0;
      stlgeometry->volumemeshed = 0;

      return NG_OK;
   }




   // generates mesh, empty mesh be already created.
   DLL_HEADER Ng_Result Ng_STL_GenerateSurfaceMesh (Ng_STL_Geometry * geom,
                                                    Ng_Mesh* mesh,
                                                    Ng_Meshing_Parameters * mp,
                                                    Ng_STL_Parameters * stlp)
   {
      STLGeometry* stlgeometry = (STLGeometry*)geom;
      Mesh* me = (Mesh*)mesh;

      // Philippose - 27/07/2009
      // Do not locally re-define "mparam" here... "mparam" is a global 
      // object
      //MeshingParameters mparam;
      mp->Transfer_Parameters();
      if (stlp) stlp->Transfer_Parameters();


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
      int retval = STLSurfaceMeshing (*stlgeometry, *me);
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
   DLL_HEADER void Ng_STL_AddTriangle (Ng_STL_Geometry * geom, 
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
   DLL_HEADER void Ng_STL_AddEdge (Ng_STL_Geometry * geom, 
      double * p1, double * p2)
   {
      readedges.Append(Point3d(p1[0],p1[1],p1[2]));
      readedges.Append(Point3d(p2[0],p2[1],p2[2]));
   }




#ifdef OCCGEOMETRY
   // --------------------- OCC Geometry / Meshing Utility Functions -------------------
   // Create new OCC Geometry Object
   DLL_HEADER Ng_OCC_Geometry * Ng_OCC_NewGeometry ()
   {
      return (Ng_OCC_Geometry*)(void*)new OCCGeometry;
   } 




   // Delete the OCC Geometry Object
   DLL_HEADER Ng_Result Ng_OCC_DeleteGeometry(Ng_OCC_Geometry * geom)
   {
      if (geom != NULL)
      {
         delete (OCCGeometry*)geom;
         geom = NULL;
         return NG_OK;
      }
      
      return NG_ERROR;
   }



   
   // Loads geometry from STEP File
   DLL_HEADER Ng_OCC_Geometry * Ng_OCC_Load_STEP (const char * filename)
   {
      // Call the STEP File Load function. Note.. the geometry class 
      // is created and instantiated within the load function
      OCCGeometry * occgeo = LoadOCC_STEP(filename);

      return ((Ng_OCC_Geometry *)occgeo);
   }



   
   // Loads geometry from IGES File
   DLL_HEADER Ng_OCC_Geometry * Ng_OCC_Load_IGES (const char * filename)
   {
      // Call the IGES File Load function. Note.. the geometry class 
      // is created and instantiated within the load function
      OCCGeometry * occgeo = LoadOCC_IGES(filename);

      return ((Ng_OCC_Geometry *)occgeo);
   }



   
   // Loads geometry from BREP File
   DLL_HEADER Ng_OCC_Geometry * Ng_OCC_Load_BREP (const char * filename)
   {
      // Call the BREP File Load function. Note.. the geometry class 
      // is created and instantiated within the load function
      OCCGeometry * occgeo = LoadOCC_BREP(filename);

      return ((Ng_OCC_Geometry *)occgeo);
   }




   // Locally limit the size of the mesh to be generated at various points 
   // based on the topology of the geometry
   DLL_HEADER Ng_Result Ng_OCC_SetLocalMeshSize (Ng_OCC_Geometry * geom,
                                                 Ng_Mesh * mesh,
                                                 Ng_Meshing_Parameters * mp)
   {
      OCCGeometry * occgeom = (OCCGeometry*)geom;
      Mesh * me = (Mesh*)mesh;

      me->geomtype = Mesh::GEOM_OCC;

      mp->Transfer_Parameters();
      
      occparam.resthcloseedgeenable = mp->closeedgeenable;
      occparam.resthcloseedgefac = mp->closeedgefact;

      // Delete the mesh structures in order to start with a clean 
      // slate
      me->DeleteMesh();

      OCCSetLocalMeshSize(*occgeom, *me);

      return(NG_OK);
   }



   
   // Mesh the edges and add Face descriptors to prepare for surface meshing
   DLL_HEADER Ng_Result Ng_OCC_GenerateEdgeMesh (Ng_OCC_Geometry * geom,
                                                 Ng_Mesh * mesh,
                                                 Ng_Meshing_Parameters * mp)
   {
      OCCGeometry * occgeom = (OCCGeometry*)geom;
      Mesh * me = (Mesh*)mesh;

      mp->Transfer_Parameters();

      OCCFindEdges(*occgeom, *me);

      if((me->GetNP()) && (me->GetNFD()))
      {
         return NG_OK;
      }
      else
      {
         return NG_ERROR;
      }
   }



   
   // Mesh the edges and add Face descriptors to prepare for surface meshing
   DLL_HEADER Ng_Result Ng_OCC_GenerateSurfaceMesh (Ng_OCC_Geometry * geom,
                                                    Ng_Mesh * mesh,
                                                    Ng_Meshing_Parameters * mp)
   {
      int numpoints = 0;

      OCCGeometry * occgeom = (OCCGeometry*)geom;
      Mesh * me = (Mesh*)mesh;

      // Set the internal meshing parameters structure from the nglib meshing 
      // parameters structure
      mp->Transfer_Parameters();


      // Only go into surface meshing if the face descriptors have already been added
      if(!me->GetNFD())
         return NG_ERROR;

      numpoints = me->GetNP();

      // Initially set up only for surface meshing without any optimisation
      int perfstepsend = MESHCONST_MESHSURFACE;

      // Check and if required, enable surface mesh optimisation step
      if(mp->optsurfmeshenable)
      {
         perfstepsend = MESHCONST_OPTSURFACE;
      }

      OCCMeshSurface(*occgeom, *me, perfstepsend);

      me->CalcSurfacesOfNode();
      
      if(me->GetNP() <= numpoints)
         return NG_ERROR;

      if(me->GetNSE() <= 0)
         return NG_ERROR;

      return NG_OK;
   }




   // Extract the face map from the OCC geometry
   // The face map basically gives an index to each face in the geometry, 
   // which can be used to access a specific face
   DLL_HEADER Ng_Result Ng_OCC_GetFMap(Ng_OCC_Geometry * geom, 
                                       Ng_OCC_TopTools_IndexedMapOfShape * FMap)
   {
      OCCGeometry* occgeom = (OCCGeometry*)geom;
      TopTools_IndexedMapOfShape *occfmap = (TopTools_IndexedMapOfShape *)FMap;

      // Copy the face map from the geometry to the given variable
      occfmap->Assign(occgeom->fmap);

      if(occfmap->Extent())
      {
         return NG_OK;
      }
      else
      {
         return NG_ERROR;
      }
   }

   // ------------------ End - OCC Geometry / Meshing Utility Functions ----------------
#endif




   // ------------------ Begin - Meshing Parameters related functions ------------------
   // Constructor for the local nglib meshing parameters class
   DLL_HEADER Ng_Meshing_Parameters :: Ng_Meshing_Parameters()
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
      
      optimize3d = "cmdmustm";
      optimize2d = "smsmsmSmSmSm";

      invert_tets = 0;
      invert_trigs = 0;

      check_overlap = 1;
      check_overlapping_boundary = 1;
   }




   // Reset the local meshing parameters to the default values
   DLL_HEADER void Ng_Meshing_Parameters :: Reset_Parameters()
   {
      (*this) = Ng_Meshing_Parameters();
   }




   // 
   DLL_HEADER void Ng_Meshing_Parameters :: Transfer_Parameters()
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
      
      if (strlen(optimize2d) > 0) mparam.optimize2d = optimize2d;
      if (strlen(optimize3d) > 0) mparam.optimize3d = optimize3d;

      mparam.inverttets = invert_tets;
      mparam.inverttrigs = invert_trigs;

      mparam.checkoverlap = check_overlap;
      mparam.checkoverlappingboundary = check_overlapping_boundary;
   }



   DLL_HEADER Ng_STL_Parameters :: Ng_STL_Parameters()
   {
      yangle = 30;
      contyangle = 20;
      
      chartangle = 10; // original = 15
      outerchartangle = 80; // original = 70;
      
      usesearchtree = 0;
      
      atlasminh = 1.0; // original = 1E-4
      
      resthatlasenable = 1;
      resthatlasfac = 2;
      
      resthchartdistenable = 1;
      resthchartdistfac = 0.3; // original = 1.2
      
      resthedgeangleenable = 0;
      resthedgeanglefac = 1;
      
      resthsurfmeshcurvenable = 1;
      resthsurfmeshcurvfac = 1;
      
      resthlinelengthenable = 1;
      resthlinelengthfac = 0.2; // original = 0.5
      
      resthcloseedgefac = 1;
      resthcloseedgeenable = 1;
   }



   DLL_HEADER void Ng_STL_Parameters :: Transfer_Parameters()
   {
      stlparam.yangle = yangle;
      stlparam.contyangle = contyangle;

      stlparam.chartangle = chartangle;
      stlparam.outerchartangle = outerchartangle;

      stlparam.usesearchtree = usesearchtree;

      stlparam.atlasminh = atlasminh;

      stlparam.resthatlasenable = resthatlasenable;
      stlparam.resthatlasfac = resthatlasfac;

      stlparam.resthchartdistenable = resthchartdistenable;
      stlparam.resthchartdistfac = resthchartdistfac;

      stlparam.resthedgeangleenable = resthedgeangleenable;
      stlparam.resthedgeanglefac = resthedgeanglefac;

      stlparam.resthsurfmeshcurvenable = resthsurfmeshcurvenable;
      stlparam.resthsurfmeshcurvfac = resthsurfmeshcurvfac;

      stlparam.resthlinelengthenable = resthlinelengthenable;
      stlparam.resthlinelengthfac = resthlinelengthfac;

      stlparam.resthcloseedgeenable = resthcloseedgeenable;
      stlparam.resthcloseedgefac = resthcloseedgefac;
   }
   // ------------------ End - Meshing Parameters related functions --------------------




   // ------------------ Begin - Second Order Mesh generation functions ----------------
   DLL_HEADER void Ng_Generate_SecondOrder(Ng_Mesh * mesh)
   {
      Refinement ref;
      ref.MakeSecondOrder(*(Mesh*) mesh);
   }




   DLL_HEADER void Ng_2D_Generate_SecondOrder(Ng_Geometry_2D * geom,
                 Ng_Mesh * mesh)
   {
      ( (SplineGeometry2d*)geom ) -> GetRefinement().MakeSecondOrder( * (Mesh*) mesh );
   }




   DLL_HEADER void Ng_STL_Generate_SecondOrder(Ng_STL_Geometry * geom,
                  Ng_Mesh * mesh)
   {
      ((STLGeometry*)geom)->GetRefinement().MakeSecondOrder(*(Mesh*) mesh);
   }




   DLL_HEADER void Ng_CSG_Generate_SecondOrder (Ng_CSG_Geometry * geom,
                  Ng_Mesh * mesh)
   {
      ((CSGeometry*)geom)->GetRefinement().MakeSecondOrder(*(Mesh*) mesh);
   }




#ifdef OCCGEOMETRY
   DLL_HEADER void Ng_OCC_Generate_SecondOrder (Ng_OCC_Geometry * geom,
                  Ng_Mesh * mesh)
   {
      ((OCCGeometry*)geom )->GetRefinement().MakeSecondOrder(*(Mesh*) mesh);
   }
#endif
   // ------------------ End - Second Order Mesh generation functions ------------------




   // ------------------ Begin - Uniform Mesh Refinement functions ---------------------
   DLL_HEADER void Ng_Uniform_Refinement (Ng_Mesh * ng_mesh)
   {
      Mesh * mesh = (Mesh*) ng_mesh;

      if (auto geom = mesh->GetGeometry())
         geom->GetRefinement().Refine (*mesh);
      else
         Refinement().Refine (*mesh);
   }


   DLL_HEADER void Ng_SetRefinementFlag (Ng_Mesh * ng_mesh, int ei, int flag)
   {
      Mesh * mesh = (Mesh*) ng_mesh;
      
      if (mesh->GetDimension() == 3)
      {
         mesh->VolumeElement(ei).SetRefinementFlag (flag != 0);
         mesh->VolumeElement(ei).SetStrongRefinementFlag (flag >= 10);
      }
      else
      {
         mesh->SurfaceElement(ei).SetRefinementFlag (flag != 0);
         mesh->SurfaceElement(ei).SetStrongRefinementFlag (flag >= 10);
      }
   }


   DLL_HEADER void Ng_SetSurfaceRefinementFlag (Ng_Mesh * ng_mesh, int ei, int flag)
   {
      Mesh * mesh = (Mesh*) ng_mesh;

      if (mesh->GetDimension() == 3)
      {
         mesh->SurfaceElement(ei).SetRefinementFlag (flag != 0);
         mesh->SurfaceElement(ei).SetStrongRefinementFlag (flag >= 10);
      }
   }


   DLL_HEADER void Ng_Refine (Ng_Mesh * ng_mesh)
   {
      Mesh * mesh = (Mesh*) ng_mesh;
      BisectionOptions biopt;
      biopt.usemarkedelements = 1;
      biopt.refine_p = 0; // only h-refinement
      biopt.refine_hp = 0;

      if (auto geom = mesh->GetGeometry())
         geom->GetRefinement().Bisect (*mesh, biopt);
      else
         Refinement().Bisect (*mesh, biopt);

      // \todo not sure if this is needed?
      //mesh -> UpdateTopology();
      //mesh -> GetCurvedElements().SetIsHighOrder (false);
   }



   DLL_HEADER void Ng_2D_Uniform_Refinement (Ng_Geometry_2D * geom,
      Ng_Mesh * mesh)
   {
      ( (SplineGeometry2d*)geom ) -> GetRefinement().Refine ( * (Mesh*) mesh );
   }



   DLL_HEADER void Ng_STL_Uniform_Refinement (Ng_STL_Geometry * geom,
      Ng_Mesh * mesh)
   {
      ( (STLGeometry*)geom ) -> GetRefinement().Refine ( * (Mesh*) mesh );
   }




   DLL_HEADER void Ng_CSG_Uniform_Refinement (Ng_CSG_Geometry * geom,
      Ng_Mesh * mesh)
   {
      ( (CSGeometry*)geom ) -> GetRefinement().Refine ( * (Mesh*) mesh );
   }




#ifdef OCCGEOMETRY
   DLL_HEADER void Ng_OCC_Uniform_Refinement (Ng_OCC_Geometry * geom,
      Ng_Mesh * mesh)
   {
      ( (OCCGeometry*)geom ) -> GetRefinement().Refine ( * (Mesh*) mesh );
   }
#endif
   // ------------------ End - Uniform Mesh Refinement functions -----------------------
} // End of namespace nglib




// compatibility functions:
namespace netgen 
{
   char geomfilename[255];

   DLL_HEADER void MyError2 (const char * ch)
   {
      cerr << ch;
   }




   //Destination for messages, errors, ...
   DLL_HEADER void Ng_PrintDest2(const char * s)
   {
#ifdef PARALLEL
     int id = 0;
     MPI_Comm_rank(MPI_COMM_WORLD, &id);
     if (id != 0) return;
#endif
     (*mycout) << s << flush;
   }


  /*
   DLL_HEADER double GetTime ()
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
void MyDummyToForceLinkingLibInterface(Mesh &mesh, NetgenGeometry &geom)
{
  netgen::WriteUserFormat("", mesh, /* geom, */ "");
}

