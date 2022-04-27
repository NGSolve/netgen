#include <mystdlib.h>
#include <myadt.hpp>
#include <occgeom.hpp>

#define OCCGEOMETRY 1
namespace nglib {
#include "nglib.h"
}


namespace netgen
{
   inline void NOOP_Deleter(void *) { ; }
   extern MeshingParameters mparam;
   DLL_HEADER extern OCCParameters occparam;
} // namespace netgen

using namespace netgen;

namespace nglib
{

   // --------------------- OCC Geometry / Meshing Utility Functions -------------------
   // Create new OCC Geometry Object
   NGLIB_API Ng_OCC_Geometry * Ng_OCC_NewGeometry ()
   {
      return (Ng_OCC_Geometry*)(void*)new OCCGeometry;
   }


   // Delete the OCC Geometry Object
   NGLIB_API Ng_Result Ng_OCC_DeleteGeometry(Ng_OCC_Geometry * geom)
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
   NGLIB_API Ng_OCC_Geometry * Ng_OCC_Load_STEP (const char * filename)
   {
      // Call the STEP File Load function. Note.. the geometry class
      // is created and instantiated within the load function
      OCCGeometry * occgeo = LoadOCC_STEP(filename);

      return ((Ng_OCC_Geometry *)occgeo);
   }


   // Loads geometry from IGES File
   NGLIB_API Ng_OCC_Geometry * Ng_OCC_Load_IGES (const char * filename)
   {
      // Call the IGES File Load function. Note.. the geometry class
      // is created and instantiated within the load function
      OCCGeometry * occgeo = LoadOCC_IGES(filename);

      return ((Ng_OCC_Geometry *)occgeo);
   }


   // Loads geometry from BREP File
   NGLIB_API Ng_OCC_Geometry * Ng_OCC_Load_BREP (const char * filename)
   {
      // Call the BREP File Load function. Note.. the geometry class
      // is created and instantiated within the load function
      OCCGeometry * occgeo = LoadOCC_BREP(filename);

      return ((Ng_OCC_Geometry *)occgeo);
   }


   // Locally limit the size of the mesh to be generated at various points
   // based on the topology of the geometry
   NGLIB_API Ng_Result Ng_OCC_SetLocalMeshSize (Ng_OCC_Geometry * geom,
                                                 Ng_Mesh * mesh,
                                                 Ng_Meshing_Parameters * mp)
   {
      OCCGeometry * occgeom = (OCCGeometry*)geom;
      Mesh * me = (Mesh*)mesh;
      me->SetGeometry( shared_ptr<NetgenGeometry>(occgeom, &NOOP_Deleter) );

      me->geomtype = Mesh::GEOM_OCC;

      mp->Transfer_Parameters();

      if(mp->closeedgeenable)
        mparam.closeedgefac = mp->closeedgefact;

      // Delete the mesh structures in order to start with a clean
      // slate
      me->DeleteMesh();

      OCCSetLocalMeshSize(*occgeom, *me, mparam, occparam);

      return(NG_OK);
   }




   // Mesh the edges and add Face descriptors to prepare for surface meshing
   NGLIB_API Ng_Result Ng_OCC_GenerateEdgeMesh (Ng_OCC_Geometry * geom,
                                                 Ng_Mesh * mesh,
                                                 Ng_Meshing_Parameters * mp)
   {
      OCCGeometry * occgeom = (OCCGeometry*)geom;
      Mesh * me = (Mesh*)mesh;
      me->SetGeometry( shared_ptr<NetgenGeometry>(occgeom, &NOOP_Deleter) );

      mp->Transfer_Parameters();

      occgeom->FindEdges(*me, mparam);

      if((me->GetNP()))
      {
         return NG_OK;
      }
      else
      {
         return NG_ERROR;
      }
   }




   // Mesh the edges and add Face descriptors to prepare for surface meshing
   NGLIB_API Ng_Result Ng_OCC_GenerateSurfaceMesh (Ng_OCC_Geometry * geom,
                                                    Ng_Mesh * mesh,
                                                    Ng_Meshing_Parameters * mp)
   {
      int numpoints = 0;

      OCCGeometry * occgeom = (OCCGeometry*)geom;
      Mesh * me = (Mesh*)mesh;
      me->SetGeometry( shared_ptr<NetgenGeometry>(occgeom, &NOOP_Deleter) );

      // Set the internal meshing parameters structure from the nglib meshing
      // parameters structure
      mp->Transfer_Parameters();

      numpoints = me->GetNP();

      // Initially set up only for surface meshing without any optimisation
      int perfstepsend = MESHCONST_MESHSURFACE;

      // Check and if required, enable surface mesh optimisation step
      if(mp->optsurfmeshenable)
      {
         perfstepsend = MESHCONST_OPTSURFACE;
      }

      occgeom->MeshSurface(*me, mparam);
      occgeom->OptimizeSurface(*me, mparam);

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
   NGLIB_API Ng_Result Ng_OCC_GetFMap(Ng_OCC_Geometry * geom,
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

   NGLIB_API void Ng_OCC_Generate_SecondOrder (Ng_OCC_Geometry * geom,
                  Ng_Mesh * mesh)
   {
      ((OCCGeometry*)geom )->GetRefinement().MakeSecondOrder(*(Mesh*) mesh);
   }

   NGLIB_API void Ng_OCC_Uniform_Refinement (Ng_OCC_Geometry * geom,
      Ng_Mesh * mesh)
   {
      ( (OCCGeometry*)geom ) -> GetRefinement().Refine ( * (Mesh*) mesh );
   }

} // namespace nglib
