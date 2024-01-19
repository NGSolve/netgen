#include <meshing.hpp>

#ifdef SOCKETS
#include "../sockets/sockets.hpp"
#endif

#include "nginterface.h"
#include "nginterface_v2.hpp"
// #include <visual.hpp>

#include "writeuser.hpp"

namespace netgen
{
  extern shared_ptr<Mesh> mesh;
}


namespace netgen
{
#define NGX_INLINE
#include "nginterface_v2_impl.hpp"

  shared_ptr<Mesh> Ngx_Mesh :: SelectMesh () const
  {
    shared_ptr<Mesh> hmesh = netgen::mesh;

    netgen::mesh = mesh;
    SetGlobalMesh (mesh);

    return hmesh;
  }
  
  Ngx_Mesh :: Ngx_Mesh (shared_ptr<Mesh> amesh)
  { mesh = amesh ? amesh : netgen::mesh; }
  Ngx_Mesh :: Ngx_Mesh (string filename, NgMPI_Comm acomm)
  { LoadMesh(filename, acomm); }
  
  Ngx_Mesh * LoadMesh (const string & filename, NgMPI_Comm comm)
  {
    netgen::mesh.reset();
    Ng_LoadMesh (filename.c_str(), comm);
    return new Ngx_Mesh (netgen::mesh);
  }

  void Ngx_Mesh :: LoadMesh (const string & filename, NgMPI_Comm comm)
  {
    netgen::mesh.reset();
    Ng_LoadMesh (filename.c_str(), comm);
    // mesh = move(netgen::mesh);
    mesh = netgen::mesh;
  }

  void Ngx_Mesh :: LoadMesh (istream & ist, NgMPI_Comm comm)
  {
    netgen::mesh = make_shared<Mesh>();
    netgen::mesh->SetCommunicator(comm);
    netgen::mesh -> Load (ist);
    // mesh = move(netgen::mesh);
    mesh = netgen::mesh;
    SetGlobalMesh (mesh);
  }

  const NgMPI_Comm & Ngx_Mesh :: GetCommunicator() const
  {
    // return Valid() ? mesh->GetCommunicator() : NgMPI_Comm{};
    if (!Valid()) throw Exception("Ngx_mesh::GetCommunicator: don't have a valid mesh");
    return mesh->GetCommunicator();
  }

  void Ngx_Mesh :: SaveMesh (ostream & ost) const
  {
    mesh -> Save (ost);
  }

  void Ngx_Mesh :: DoArchive (Archive & archive)
  {
#ifdef PARALLEL
    if (archive.Input()) {
      mesh = make_shared<Mesh>();
      mesh->SetCommunicator(GetCommunicator());
    }
#endif
    mesh->DoArchive(archive);
    if (archive.Input())
      {
        netgen::mesh = mesh;
        SetGlobalMesh (mesh);
      }
    /*
    if (archive.Output())
      {
        stringstream str;
        SaveMesh (str);
        string st = str.str();
        archive & st;
      }
    else
      {
        string st;
        archive & st;
        stringstream str(st);
        LoadMesh (str);
      }
    */
  }

  void Ngx_Mesh :: UpdateTopology ()
  {
    if (mesh)
      mesh -> UpdateTopology();
  }


  /*
  Ngx_Mesh :: Ngx_Mesh (Mesh * amesh)
    : mesh(amesh)
  { ; }
  */

  Ngx_Mesh :: ~Ngx_Mesh ()
  { 
    // causes crashes when global variable netgen::mesh is destructed
    // before visualization data 
    
    if (mesh == netgen::mesh)
      netgen::mesh = nullptr;
  }

  int Ngx_Mesh :: GetDimension() const
  {
    return mesh -> GetDimension();
  }

  int Ngx_Mesh :: GetNLevels() const
  {
    return max(size_t(1), mesh -> level_nv.Size());
  }

  size_t Ngx_Mesh :: GetNVLevel(int level) const
  {
    if (level >= mesh->level_nv.Size())
      return mesh->GetNV();
    else
      return mesh->level_nv[level];
  }
  
  int Ngx_Mesh :: GetNElements (int dim) const
  {
    switch (dim)
      {
      case 0: return mesh -> pointelements.Size();
      case 1: return mesh -> GetNSeg();
      case 2: return mesh -> GetNSE();
      case 3: return mesh -> GetNE();
      }
    return -1;
  }
  
  int Ngx_Mesh :: GetNNodes (int nt) const
  {
    switch (nt)
      {
      case 0: return mesh -> GetNV();
      case 1: return mesh->GetTopology().GetNEdges();
      case 2: return mesh->GetTopology().GetNFaces();
      case 3: return mesh -> GetNE();
      }
    return -1;
  }

  /*
  Ng_Point Ngx_Mesh :: GetPoint (int nr) const
  {
    return Ng_Point (&mesh->Point(nr + PointIndex::BASE)(0));
  }
  */

  /*
  template <> DLL_HEADER Ng_Element Ngx_Mesh :: GetElement<0> (int nr) const
  {
    const Element0d & el = mesh->pointelements[nr];
    
    Ng_Element ret;
    ret.type = NG_PNT;
    ret.index = el.index;

    ret.points.num = 1;
    ret.points.ptr = (int*)&el.pnum;

    ret.vertices.num = 1;
    ret.vertices.ptr = (int*)&el.pnum;

    ret.edges.num = 0;
    ret.edges.ptr = NULL;

    ret.faces.num = 0;
    ret.faces.ptr = NULL;

    return ret;
  }
  */
  /*
  template <> DLL_HEADER Ng_Element Ngx_Mesh :: GetElement<1> (int nr) const
  {
    const Segment & el = mesh->LineSegment (SegmentIndex(nr));

    Ng_Element ret;
    ret.type = NG_ELEMENT_TYPE(el.GetType());

    ret.points.num = el.GetNP();
    ret.points.ptr = (int*)&(el[0]);

    ret.vertices.num = 2;
    ret.vertices.ptr = (int*)&(el[0]);

    ret.edges.num = 1;
    ret.edges.ptr = mesh->GetTopology().GetSegmentElementEdgesPtr (nr);

    ret.faces.num = 0;
    ret.faces.ptr = NULL;

    return ret;
  }

  template <> DLL_HEADER Ng_Element Ngx_Mesh :: GetElement<2> (int nr) const
  {
    const Element2d & el = mesh->SurfaceElement (SurfaceElementIndex (nr));
  
    Ng_Element ret;
    ret.type = NG_ELEMENT_TYPE(el.GetType());
    ret.points.num = el.GetNP();
    ret.points.ptr  = (int*)&el[0];

    ret.vertices.num = el.GetNV();
    ret.vertices.ptr = (int*)&(el[0]);

    ret.edges.num = MeshTopology::GetNEdges (el.GetType());
    ret.edges.ptr = mesh->GetTopology().GetSurfaceElementEdgesPtr (nr);

    ret.faces.num = MeshTopology::GetNFaces (el.GetType());
    ret.faces.ptr = mesh->GetTopology().GetSurfaceElementFacesPtr (nr);

    return ret;
  }

  template <> DLL_HEADER Ng_Element Ngx_Mesh :: GetElement<3> (int nr) const
  {
    const Element & el = mesh->VolumeElement (ElementIndex (nr));
  
    Ng_Element ret;
    ret.type = NG_ELEMENT_TYPE(el.GetType());
    ret.points.num = el.GetNP();
    ret.points.ptr = (int*)&el[0];

    ret.vertices.num = el.GetNV();
    ret.vertices.ptr = (int*)&(el[0]);

    ret.edges.num = MeshTopology::GetNEdges (el.GetType());
    ret.edges.ptr = mesh->GetTopology().GetElementEdgesPtr (nr);

    ret.faces.num = MeshTopology::GetNFaces (el.GetType());
    ret.faces.ptr = mesh->GetTopology().GetElementFacesPtr (nr);

    return ret;
  }
  */

  
  /*
  template <>
  DLL_HEADER int Ngx_Mesh :: GetElementIndex<0> (int nr) const
  {
    return 0;
  }

  template <>
  DLL_HEADER int Ngx_Mesh :: GetElementIndex<1> (int nr) const
  {
    return (*mesh)[SegmentIndex(nr)].si;
  }
  
  template <>
  DLL_HEADER int Ngx_Mesh :: GetElementIndex<2> (int nr) const
  {
    int ind = (*mesh)[SurfaceElementIndex(nr)].GetIndex(); 
    return mesh->GetFaceDescriptor(ind).BCProperty();
  }
  
  template <>
  DLL_HEADER int Ngx_Mesh :: GetElementIndex<3> (int nr) const
  {
    return (*mesh)[ElementIndex(nr)].GetIndex();
  }
  */

















  /*
  DLL_HEADER Ng_Point Ng_GetPoint (int nr)
  {
    Ng_Point ret;
    ret.pt = &mesh->Point(nr + PointIndex::BASE)(0);
    return ret;
  }


  template <>
  DLL_HEADER int Ng_GetElementIndex<1> (int nr)
  {
    return (*mesh)[SegmentIndex(nr)].si;
  }
  
  template <>
  DLL_HEADER int Ng_GetElementIndex<2> (int nr)
  {
    int ind = (*mesh)[SurfaceElementIndex(nr)].GetIndex(); 
    return mesh->GetFaceDescriptor(ind).BCProperty();
  }
  
  template <>
  DLL_HEADER int Ng_GetElementIndex<3> (int nr)
  {
    return (*mesh)[ElementIndex(nr)].GetIndex();
  }

  template <> int DLL_HEADER Ng_GetNElements<0> ()
  {
    return 0;
  }

  template <> int DLL_HEADER Ng_GetNElements<1> ()
  {
    return mesh->GetNSeg();
  }

  template <> DLL_HEADER int Ng_GetNElements<2> ()
  {
    return mesh->GetNSE();
  }

  template <> DLL_HEADER int Ng_GetNElements<3> ()
  {
    return mesh->GetNE();
  }

  template <> DLL_HEADER Ng_Element Ng_GetElement<0> (int nr)
  {
    cout << "Netgen does not support 0-D elements" << endl;
	Ng_Element ret;
	return ret;
  }

  template <> DLL_HEADER Ng_Element Ng_GetElement<1> (int nr)
  {
    const Segment & el = mesh->LineSegment (SegmentIndex(nr));

    Ng_Element ret;
    ret.type = NG_ELEMENT_TYPE(el.GetType());

    ret.points.num = el.GetNP();
    ret.points.ptr = (int*)&(el[0]);

    ret.vertices.num = 2;
    ret.vertices.ptr = (int*)&(el[0]);

    ret.edges.num = 1;
    ret.edges.ptr = mesh->GetTopology().GetSegmentElementEdgesPtr (nr);

    ret.faces.num = 0;
    ret.faces.ptr = NULL;

    return ret;
  }

  template <> DLL_HEADER Ng_Element Ng_GetElement<2> (int nr)
  {
    const Element2d & el = mesh->SurfaceElement (SurfaceElementIndex (nr));
  
    Ng_Element ret;
    ret.type = NG_ELEMENT_TYPE(el.GetType());
    ret.points.num = el.GetNP();
    ret.points.ptr  = (int*)&el[0];

    ret.vertices.num = el.GetNV();
    ret.vertices.ptr = (int*)&(el[0]);

    ret.edges.num = MeshTopology::GetNEdges (el.GetType());
    ret.edges.ptr = mesh->GetTopology().GetSurfaceElementEdgesPtr (nr);

    ret.faces.num = MeshTopology::GetNFaces (el.GetType());
    ret.faces.ptr = mesh->GetTopology().GetSurfaceElementFacesPtr (nr);

    return ret;
  }

  template <> DLL_HEADER Ng_Element Ng_GetElement<3> (int nr) 
  {
    const Element & el = mesh->VolumeElement (ElementIndex (nr));
  
    Ng_Element ret;
    ret.type = NG_ELEMENT_TYPE(el.GetType());
    ret.points.num = el.GetNP();
    ret.points.ptr = (int*)&el[0];

    ret.vertices.num = el.GetNV();
    ret.vertices.ptr = (int*)&(el[0]);

    ret.edges.num = MeshTopology::GetNEdges (el.GetType());
    ret.edges.ptr = mesh->GetTopology().GetElementEdgesPtr (nr);

    ret.faces.num = MeshTopology::GetNFaces (el.GetType());
    ret.faces.ptr = mesh->GetTopology().GetElementFacesPtr (nr);

    return ret;
  }
  */







  
  
  template <> DLL_HEADER void Ngx_Mesh :: 
  ElementTransformation<3,3> (int elnr, 
                              const double * xi,
                              double * x, 
                              double * dxdxi) const
  {
    Point<3> xl(xi[0], xi[1], xi[2]);
    Point<3> xg;
    Mat<3,3> dx;
    mesh->GetCurvedElements().CalcElementTransformation (xl, elnr, xg, dx);

    if (x)
      for (int i = 0; i < 3; i++) x[i] = xg(i);
    
    if (dxdxi)
      for (int i=0; i<3; i++)
        {
          dxdxi[3*i] = dx(i,0);
          dxdxi[3*i+1] = dx(i,1);
          dxdxi[3*i+2] = dx(i,2);
        }
  }

  template <> DLL_HEADER void Ngx_Mesh ::
  ElementTransformation<2,3> (int elnr,
                              const double * xi, 
                              double * x, 
                              double * dxdxi) const
  {
    Point<2> xl(xi[0], xi[1]);
    Point<3> xg;
    Mat<3,2> dx;

    mesh->GetCurvedElements().CalcSurfaceTransformation (xl, elnr, xg, dx);
    
    if (x)
      for (int i = 0; i < 3; i++) x[i] = xg(i);
    
    if (dxdxi)
      for (int i=0; i<3; i++)
        {
          dxdxi[2*i] = dx(i,0);
          dxdxi[2*i+1] = dx(i,1);
        }
  }

  template <> DLL_HEADER void Ngx_Mesh ::
  ElementTransformation<1,3> (int elnr,
			      const double * xi,
			      double * x,
			      double * dxdxi) const
  {
    Point<3> xg;
    Vec<3> dx;
    mesh->GetCurvedElements().CalcSegmentTransformation(xi[0],elnr,xg,dx);
    if(x)
      for(int i=0;i<3;i++) x[i] = xg(i);

    if(dxdxi)
      for(int i=0;i<3;i++) dxdxi[i] = dx(i);
  }

  template <> DLL_HEADER void Ngx_Mesh ::
  ElementTransformation<0,3> (int elnr,
			      const double * xi,
			      double * x,
			      double * dxdxi) const
  {
    PointIndex pi = mesh->pointelements[elnr].pnum;
    Point<3> xg = mesh->Point(pi);
    if (x)
      for(int i=0;i<3;i++) x[i] = xg(i);
  }

  
  template <> DLL_HEADER void Ngx_Mesh ::
  ElementTransformation<2,2> (int elnr,
                                 const double * xi, 
                                 double * x, 
                                 double * dxdxi) const
  {
    Point<2> xl(xi[0], xi[1]);
    Point<3> xg;
    Mat<3,2> dx;

    mesh->GetCurvedElements().CalcSurfaceTransformation (xl, elnr, xg, dx);
    
    if (x)
      for (int i = 0; i < 2; i++) x[i] = xg(i);
    
    if (dxdxi)
      for (int i=0; i<2; i++)
        {
          dxdxi[2*i] = dx(i,0);
          dxdxi[2*i+1] = dx(i,1);
        }
  }




  template <> DLL_HEADER void Ngx_Mesh :: 
  ElementTransformation<1,2> (int elnr,
                              const double * xi,
                              double * x,
                              double * dxdxi) const
  {
    Point<3> xg;
    Vec<3> dx;

    mesh->GetCurvedElements().CalcSegmentTransformation (xi[0], elnr, xg, dx);
    
    if (x)
      for (int i = 0; i < 2; i++) x[i] = xg(i);
    
    if (dxdxi)
      for (int i=0; i < 2; i++)
        dxdxi[i] = dx(i);
  }



  template <> DLL_HEADER void Ngx_Mesh :: 
  ElementTransformation<1,1> (int elnr,
                              const double * xi,
                              double * x,
                              double * dxdxi) const
  {
    Point<3> xg;
    Vec<3> dx;

    mesh->GetCurvedElements().CalcSegmentTransformation (xi[0], elnr, xg, dx);
    
    if (x) x[0] = xg(0);
    if (dxdxi) dxdxi[0] = dx(0);
  }

  template <> DLL_HEADER void Ngx_Mesh ::
  ElementTransformation<0,2> (int elnr,
			      const double *xi,
			      double * x,
			      double * dxdxi) const
  {
    PointIndex pnum = mesh->pointelements[elnr].pnum;
    if (x)
      for (int i = 0; i< 2; i++) x[i] = (*mesh)[pnum](i);
  }


  template <> DLL_HEADER void Ngx_Mesh :: 
  ElementTransformation<0,1> (int elnr, 
                              const double * xi,
                              double * x,
                              double * dxdxi) const
  {
    PointIndex pnum = mesh->pointelements[elnr].pnum;
    if (x) x[0] = (*mesh)[pnum](0);
    // if (dxdxi) dxdxi[0] = 0;
    // Jacobi-matrix is 1 x 0 !!!
  }





  
  
  template <> DLL_HEADER void Ngx_Mesh :: 
  MultiElementTransformation<3,3> (int elnr, int npts,
                                   const double * xi, size_t sxi,
                                   double * x, size_t sx,
                                   double * dxdxi, size_t sdxdxi) const
  {
    mesh->GetCurvedElements().CalcMultiPointElementTransformation (elnr, npts, xi, sxi, x, sx, dxdxi, sdxdxi);
  }
  
  template <> DLL_HEADER void Ngx_Mesh ::
  MultiElementTransformation<2,2> (int elnr, int npts,
                                   const double * xi, size_t sxi,
                                   double * x, size_t sx,
                                   double * dxdxi, size_t sdxdxi) const
  {
    mesh->GetCurvedElements().CalcMultiPointSurfaceTransformation<2> (elnr, npts, xi, sxi, x, sx, dxdxi, sdxdxi);
  }

  template <> DLL_HEADER void Ngx_Mesh :: 
  MultiElementTransformation<2,3> (int elnr, int npts,
                                   const double * xi, size_t sxi,
                                   double * x, size_t sx,
                                   double * dxdxi, size_t sdxdxi) const
  {
    mesh->GetCurvedElements().CalcMultiPointSurfaceTransformation<3> (elnr, npts, xi, sxi, x, sx, dxdxi, sdxdxi);
  }

  template <> DLL_HEADER void Ngx_Mesh ::
  MultiElementTransformation<1,3> (int elnr, int npts,
				   const double * xi, size_t sxi,
                                   double * x, size_t sx,
                                   double * dxdxi, size_t sdxdxi) const
  {
    mesh->GetCurvedElements().CalcMultiPointSegmentTransformation<3> (elnr, npts, xi, sxi, x, sx, dxdxi, sdxdxi);
  }

  template <> DLL_HEADER void Ngx_Mesh ::
  MultiElementTransformation<0,3> (int elnr, int npts,
				   const double * xi, size_t sxi,
                                   double * x, size_t sx,
                                   double * dxdxi, size_t sdxdxi) const
  {
    for (int i = 0; i < npts; i++)
      ElementTransformation<0,3> (elnr, xi+i*sxi, x+i*sx, dxdxi+i*sdxdxi);
  }
    
  template <> DLL_HEADER void Ngx_Mesh :: 
  MultiElementTransformation<1,2> (int elnr, int npts,
                                   const double * xi, size_t sxi,
                                   double * x, size_t sx,
                                   double * dxdxi, size_t sdxdxi) const
  {
    mesh->GetCurvedElements().CalcMultiPointSegmentTransformation<2> (elnr, npts, xi, sxi, x, sx, dxdxi, sdxdxi);
  }

  template <> DLL_HEADER void Ngx_Mesh :: 
  MultiElementTransformation<1,1> (int elnr, int npts,
                                   const double * xi, size_t sxi,
                                   double * x, size_t sx,
                                   double * dxdxi, size_t sdxdxi) const
  {
    for (int i = 0; i < npts; i++)
      ElementTransformation<1,1> (elnr, xi + i*sxi, x+i*sx, dxdxi+i*sdxdxi);
  }

  template <> DLL_HEADER void Ngx_Mesh ::
  MultiElementTransformation<0,2> (int elnr, int npts,
                                   const double * xi, size_t sxi,
                                   double * x, size_t sx,
                                   double * dxdxi, size_t sdxdxi) const
  {
    for (int i = 0; i < npts; i++)
      ElementTransformation<0,2> (elnr, xi + i*sxi, x+i*sx, dxdxi+i*sdxdxi);
  }


  template <> DLL_HEADER void Ngx_Mesh :: 
  MultiElementTransformation<0,1> (int elnr, int npts,
                                   const double * xi, size_t sxi,
                                   double * x, size_t sx,
                                   double * dxdxi, size_t sdxdxi) const
  {
    for (int i = 0; i < npts; i++)
      ElementTransformation<0,1> (elnr, xi + i*sxi, x+i*sx, dxdxi+i*sdxdxi);
  }


  int Ngx_Mesh :: GetHPElementLevel (int ei, int dir) const
  {
    ei++;
    int level = -1;
    
    if (mesh->hpelements)
      {
	int hpelnr = -1;
	if (mesh->GetDimension() == 2)
	  hpelnr = mesh->SurfaceElement(ei).GetHpElnr();
	else
	  hpelnr = mesh->VolumeElement(ei).GetHpElnr();

        if (hpelnr < 0)
          throw NgException("Ngx_Mesh::GetHPElementLevel: Wrong hp-element number!");
        
        if (dir == 1)
          level = (*mesh->hpelements)[hpelnr].levelx;
        else if (dir == 2)
          level = (*mesh->hpelements)[hpelnr].levely;
        else if (dir == 3)
          level = (*mesh->hpelements)[hpelnr].levelz;
        else
          throw NgException("Ngx_Mesh::GetHPElementLevel: dir has to be 1, 2 or 3!");
      }
    //else
    //  throw NgException("Ngx_Mesh::GetHPElementLevel only for HPRefinement implemented!");

    return level;	  
  }
  
  int Ngx_Mesh :: GetParentElement (int ei) const
  {
      ei++;
      if (mesh->GetDimension() == 3)
      {
          if (ei <= mesh->mlparentelement.Size())
              return mesh->mlparentelement.Get(ei)-1;
      }
      else
      {
          if (ei <= mesh->mlparentsurfaceelement.Size())
              return mesh->mlparentsurfaceelement.Get(ei)-1;
      }
      return -1;
  }


  int Ngx_Mesh :: GetParentSElement (int ei) const
  {
      ei++;
      if (mesh->GetDimension() == 3)
      {
          if (ei <= mesh->mlparentsurfaceelement.Size())
              return mesh->mlparentsurfaceelement.Get(ei)-1;
      }
      else
      {
          return -1;
      }
      return -1;
  }

  int Ngx_Mesh :: GetNIdentifications () const
  {
    return mesh->GetIdentifications().GetMaxNr();
  }

  int Ngx_Mesh :: GetIdentificationType(int idnr) const
  {
    return mesh->GetIdentifications().GetType(idnr+1);
  }

  Ng_BufferMS<int,4> Ngx_Mesh::GetFaceEdges (int fnr) const
  {
    const MeshTopology & topology = mesh->GetTopology();
    NgArrayMem<int,4> ia;
    topology.GetFaceEdges (fnr+1, ia);
    Ng_BufferMS<int,4> res(ia.Size());
    for (size_t i = 0; i < ia.Size(); i++)
      res[i] = ia[i]-1;
    return res;
  }




  template<> DLL_HEADER void Ngx_Mesh :: 
  MultiElementTransformation<1,1> (int elnr, int npts,
                                   const SIMD<double> * xi, size_t sxi,
                                   SIMD<double> * x, size_t sx,
                                   SIMD<double> * dxdxi, size_t sdxdxi) const
  {
    cout << "multi-eltrafo simd called, 1,1,simd" << endl;
  }

  template<> DLL_HEADER void Ngx_Mesh :: 
  MultiElementTransformation<2,2> (int elnr, int npts,
                                   const SIMD<double> * xi, size_t sxi,
                                   SIMD<double> * x, size_t sx,
                                   SIMD<double> * dxdxi, size_t sdxdxi) const
  {
    mesh->GetCurvedElements().CalcMultiPointSurfaceTransformation<2>
      (elnr, npts, xi, sxi, x, sx, dxdxi, sdxdxi);
    /*
    for (int i = 0; i < npts; i++)
      {
        double hxi[4][2];
        double hx[4][2];
        double hdxdxi[4][4];
        for (int j = 0; j < 4; j++)
          for (int k = 0; k < 2; k++)
            hxi[j][k] = ((double*)&(xi[k]))[j];
        MultiElementTransformation<2,2> (elnr, 4, &hxi[0][0], 2, &hx[0][0], 2, &hdxdxi[0][0], 4);
        for (int j = 0; j < 4; j++)
          for (int k = 0; k < 2; k++)
            ((double*)&(x[k]))[j] = hx[j][k];
        for (int j = 0; j < 4; j++)
          for (int k = 0; k < 4; k++)
            ((double*)&(dxdxi[k]))[j] = hdxdxi[j][k];
        
        xi += sxi;
        x += sx;
        dxdxi += sdxdxi;
      }
    */
  }

  template<> DLL_HEADER void Ngx_Mesh :: 
  MultiElementTransformation<3,3> (int elnr, int npts,
                                   const SIMD<double> * xi, size_t sxi,
                                   SIMD<double> * x, size_t sx,
                                   SIMD<double> * dxdxi, size_t sdxdxi) const
  {
    mesh->GetCurvedElements().CalcMultiPointElementTransformation
      (elnr, npts,
       xi, sxi,
       x, sx,
       dxdxi, sdxdxi);
    /*
    for (int i = 0; i < npts; i++)
      {
        double hxi[4][3];
        double hx[4][3];
        double hdxdxi[4][9];
        for (int j = 0; j < 4; j++)
          for (int k = 0; k < 3; k++)
            hxi[j][k] = ((double*)&(xi[k]))[j];
        MultiElementTransformation<3,3> (elnr, 4, &hxi[0][0], 3, &hx[0][0], 3, &hdxdxi[0][0], 9);
        for (int j = 0; j < 4; j++)
          for (int k = 0; k < 3; k++)
            ((double*)&(x[k]))[j] = hx[j][k];
        for (int j = 0; j < 4; j++)
          for (int k = 0; k < 9; k++)
            ((double*)&(dxdxi[k]))[j] = hdxdxi[j][k];
        
        xi += sxi;
        x += sx;
        dxdxi += sdxdxi;
      }
    */
  }

  template<> DLL_HEADER void Ngx_Mesh ::
  MultiElementTransformation<0,2> (int elnr, int npts,
				   const SIMD<double> *xi, size_t sxi,
				   SIMD<double> * x, size_t sx,
				   SIMD<double> * dxdxi, size_t sdxdxi) const
  {
    //cout << "MultiElementtransformation<0,2> simd not implemented" << endl;

    PointIndex pi = mesh->pointelements[elnr].pnum;
    Point<3> xg = mesh->Point(pi);
    if (x)
      for (int j = 0; j < npts; j++)
	for (int i = 0; i < 2; i++)
	  x[j*sx+i] = xg(i);
  }
  
  template<> DLL_HEADER void Ngx_Mesh :: 
  MultiElementTransformation<0,1> (int elnr, int npts,
                                   const SIMD<double> * xi, size_t sxi,
                                   SIMD<double> * x, size_t sx,
                                   SIMD<double> * dxdxi, size_t sdxdxi) const
  {
    //cout << "multi-eltrafo simd called, 0,1,simd" << endl;
    PointIndex pi = mesh->pointelements[elnr].pnum;
    Point<3> xg = mesh->Point(pi);
    if (x)
      for (int j = 0; j < npts; j++)
	for (int i = 0; i < 1; i++)
	  x[j*sx+i] = xg(i);
  }

  template<> DLL_HEADER void Ngx_Mesh :: 
  MultiElementTransformation<1,3> (int elnr, int npts,
                                   const SIMD<double> * xi, size_t sxi,
                                   SIMD<double> * x, size_t sx,
                                   SIMD<double> * dxdxi, size_t sdxdxi) const
  {
    mesh->GetCurvedElements().CalcMultiPointSegmentTransformation<3>
      (elnr, npts, xi, sxi, x, sx, dxdxi, sdxdxi);
    /*
    double hxi[4][1];
    double hx[4][3];
    double hdxdxi[4][3];
    for (int j = 0; j<4;j++)
      hxi[j][0] = ((double*)&(xi[0]))[j];
    MultiElementTransformation<1,3> (elnr, 4, &hxi[0][0], 1, &hx[0][0], 3, &hdxdxi[0][0],3);
    for(int j=0; j<4; j++)
      for(int k=0; k<3; k++)
	((double*)&(x[k]))[j] = hx[j][k];
    for(int j=0; j< 4; j++)
      for (int k = 0; k<3; k++)
	((double*) & (dxdxi[k]))[j] = hdxdxi[j][k];

    xi += sxi;
    x += sx;
    dxdxi += sdxdxi;
    */
  }
  
  template<> DLL_HEADER void Ngx_Mesh :: 
  MultiElementTransformation<1,2> (int elnr, int npts,
                                   const SIMD<double> * xi, size_t sxi,
                                   SIMD<double> * x, size_t sx,
                                   SIMD<double> * dxdxi, size_t sdxdxi) const
  {
    mesh->GetCurvedElements().CalcMultiPointSegmentTransformation<2>
      (elnr, npts, xi, sxi, x, sx, dxdxi, sdxdxi);
    /*
    for (int i = 0; i < npts; i++)
      {
        double hxi[4][1];
        double hx[4][2];
        double hdxdxi[4][2];
        for (int j = 0; j < 4; j++)
          for (int k = 0; k < 1; k++)
            hxi[j][k] = ((double*)&(xi[k]))[j];
        MultiElementTransformation<1,2> (elnr, 4, &hxi[0][0], 1, &hx[0][0], 2, &hdxdxi[0][0], 2);
        for (int j = 0; j < 4; j++)
          for (int k = 0; k < 2; k++)
            ((double*)&(x[k]))[j] = hx[j][k];
        for (int j = 0; j < 4; j++)
          for (int k = 0; k < 2; k++)
            ((double*)&(dxdxi[k]))[j] = hdxdxi[j][k];
        
        xi += sxi;
        x += sx;
        dxdxi += sdxdxi;
      }
    */
  }

  template<> DLL_HEADER void Ngx_Mesh :: 
  MultiElementTransformation<2,3> (int elnr, int npts,
                                   const SIMD<double> * xi, size_t sxi,
                                   SIMD<double> * x, size_t sx,
                                   SIMD<double> * dxdxi, size_t sdxdxi) const
  {
    mesh->GetCurvedElements().CalcMultiPointSurfaceTransformation<3>
      (elnr, npts, xi, sxi, x, sx, dxdxi, sdxdxi);
    /*
    for (int i = 0; i < npts; i++)
      {
        double hxi[4][2];
        double hx[4][3];
        double hdxdxi[4][6];
        for (int j = 0; j < 4; j++)
          for (int k = 0; k < 2; k++)
            hxi[j][k] = ((double*)&(xi[k]))[j];
        MultiElementTransformation<2,3> (elnr, 4, &hxi[0][0], 2, &hx[0][0], 3, &hdxdxi[0][0], 6);
        for (int j = 0; j < 4; j++)
          for (int k = 0; k < 3; k++)
            ((double*)&(x[k]))[j] = hx[j][k];
        for (int j = 0; j < 4; j++)
          for (int k = 0; k < 6; k++)
            ((double*)&(dxdxi[k]))[j] = hdxdxi[j][k];
        
        xi += sxi;
        x += sx;
        dxdxi += sdxdxi;
      }
    */
  }

  template<> DLL_HEADER void Ngx_Mesh :: 
  MultiElementTransformation<0,3> (int elnr, int npts,
                                   const SIMD<double> * xi, size_t sxi,
                                   SIMD<double> * x, size_t sx,
                                   SIMD<double> * dxdxi, size_t sdxdxi) const
  {
    for (int i = 0; i < npts; i++)
      {
        double hxi[4][1];
        double hx[4][3];
        for (int j = 0; j < 4; j++)
          for (int k = 0; k < 1; k++)
            hxi[j][k] = ((double*)&(xi[k]))[j];
        MultiElementTransformation<0,3> (elnr, 4, &hxi[0][0], 2, &hx[0][0], 3, (double*)nullptr, 0);
        for (int j = 0; j < 4; j++)
          for (int k = 0; k < 3; k++)
            ((double*)&(x[k]))[j] = hx[j][k];
        xi += sxi;
        x += sx;
        dxdxi += sdxdxi;
      }
  }



  


  template <>
  DLL_HEADER int Ngx_Mesh :: FindElementOfPoint <1> 
  (double * hp, double * lami,
   bool build_searchtree, 
   int * const indices, int numind) const

  {
    switch (mesh->GetDimension())
      {
      case 1:
        {
          Point<3> p(hp[0], 0,0);
          for (SegmentIndex si = 0; si < mesh->GetNSeg(); si++)
            {
              auto & seg = (*mesh)[si];
              Point<3> p1 = (*mesh)[seg[0]];
              Point<3> p2 = (*mesh)[seg[1]];
              double lam = (p(0)-p1(0)) / (p2(0)-p1(0));
              if (lam >= -1e-10 && lam <= 1+1e-10)
                {
                  lami[0] = 1-lam;
                  return si;
                }
            }
        }
        break;
      case 2:
        {
          Point<3> p(hp[0], hp[1],0);
          try
            {
              auto ind = mesh->GetSurfaceElementOfPoint(p, lami, nullptr,
                                                        build_searchtree);
              return ind - 1;
            }
          catch(NgException e) // quads not implemented curved yet
            {
          for (SegmentIndex si = 0; si < mesh->GetNSeg(); si++)
            {
              auto & seg = (*mesh)[si];
              Point<3> p1 = (*mesh)[seg[0]];
              Point<3> p2 = (*mesh)[seg[1]];
              double lam;
              double r;
              if (fabs(p2[0]-p1[0]) >= fabs(p2[1]-p1[1]))
                {
                  lam = (p[0]-p1[0])/(p2[0]-p1[0]);
                  r = p[1] - p1[1] - lam*(p2[1]-p1[1]);
                }
              else
                {
                  lam = (p[1]-p1[1])/(p2[1]-p1[1]);
                  r = p[0] - p1[0] - lam*(p2[0]-p1[0]);
                }
              if ( lam >= -1e-10 && lam <= 1+1e-10 && fabs(r) <= 1e-10 )
                {
                  lami[0] = 1-lam;
                  return si;
                }
            }
            }
        }
        break;
      case 3:
      default:
        throw Exception("FindElementOfPoint<1> only implemented for mesh-dimension 1 and 2!");
        break;
      }
 
    return -1;
  }

  template <>
  DLL_HEADER int Ngx_Mesh :: FindElementOfPoint <2> 
  (double * p, double * lami,
   bool build_searchtree, 
   int * const indices, int numind) const

  {
    NgArray<int> dummy(numind);
    for (int i = 0; i < numind; i++) dummy[i] = indices[i]+1;
    
    double lam3[3];
    int ind;

    if (mesh->GetDimension() == 2)
      {
        Point<3> p2d(p[0], p[1], 0);
        ind = mesh->GetElementOfPoint(p2d, lam3, &dummy, build_searchtree);
      }
    else
      {
        Point3d p3d(p[0], p[1], p[2]);
        ind = mesh->GetSurfaceElementOfPoint(p3d, lam3, &dummy, build_searchtree);
      }
    
    if (ind > 0)
      {
        if(mesh->SurfaceElement(ind).GetType()==QUAD || mesh->SurfaceElement(ind).GetType()==TRIG6)
          {
            lami[0] = lam3[0];
            lami[1] = lam3[1];
          }
        else 
          {
            lami[0] = 1-lam3[0]-lam3[1];
            lami[1] = lam3[0];
          }
      }
    return ind-1;
  }


  template <>
  DLL_HEADER int Ngx_Mesh :: FindElementOfPoint <3>
  (double * p, double * lami,
   bool build_searchtree, 
   int * const indices, int numind) const

  {
    NgArray<int> dummy(numind);
    for (int i = 0; i < numind; i++) dummy[i] = indices[i]+1;
    
    Point<3> p3d(p[0], p[1], p[2]);
    int ind = 
      mesh->GetElementOfPoint(p3d, lami, &dummy, build_searchtree);
    return ind-1;
  }

  void Ngx_Mesh :: Curve (int order)
  {
    NgLock meshlock (mesh->MajorMutex(), true);
    mesh->BuildCurvedElements(order);
  }

  int Ngx_Mesh :: GetCurveOrder ()
  {
    return mesh->GetCurvedElements().GetOrder();
  }

  void Ngx_Mesh :: EnableTable (string name, bool set)
  {
    mesh->GetTopology().EnableTable (name, set);
    mesh->SetNextTimeStamp();  // update topology will do work
  }


  
  template <>
  DLL_HEADER void Ngx_Mesh :: SetRefinementFlag<2> (size_t elnr, bool flag)
  {
    mesh->SurfaceElement(elnr+1).SetRefinementFlag(flag);
  }

  template <>
  DLL_HEADER void Ngx_Mesh :: SetRefinementFlag<3> (size_t elnr, bool flag)
  {
    mesh->VolumeElement(elnr+1).SetRefinementFlag(flag);    
  }
  
  void Ngx_Mesh :: Refine (NG_REFINEMENT_TYPE reftype, bool onlyonce,
                           void (*task_manager)(function<void(int,int)>),
                           NgTracer tracer)
  {
    NgLock meshlock (mesh->MajorMutex(), 1);
    
    BisectionOptions biopt;
    biopt.usemarkedelements = 1;
    biopt.refine_p = 0;
    biopt.refine_hp = 0;
    biopt.onlyonce = onlyonce;
    if (reftype == NG_REFINE_P)
      biopt.refine_p = 1;
    if (reftype == NG_REFINE_HP)
      biopt.refine_hp = 1;
    biopt.task_manager = task_manager;
    biopt.tracer = tracer;

    mesh->GetGeometry()->GetRefinement().Bisect (*mesh, biopt);
    (*tracer)("call updatetop", false);
    mesh -> UpdateTopology(task_manager, tracer);
    (*tracer)("call updatetop", true);
    if(mesh->GetCurvedElements().IsHighOrder())
      mesh->GetCurvedElements()
        .BuildCurvedElements(&mesh->GetGeometry()->GetRefinement(),
                                mesh->GetCurvedElements().GetOrder());
  }





  // just copied with redesign

  size_t Ngx_Mesh::GetNP() const
  {
    return mesh->GetNP();
  }

  
  int Ngx_Mesh::GetSurfaceElementSurfaceNumber (size_t ei) const
  {
    if (mesh->GetDimension() == 3)
      return mesh->GetFaceDescriptor(mesh->SurfaceElement(ei).GetIndex()).SurfNr();
    else
      return mesh->LineSegment(ei).si;
  }
  int Ngx_Mesh::GetSurfaceElementFDNumber (size_t ei) const
  {
    if (mesh->GetDimension() == 3)
      return mesh->SurfaceElement(ei).GetIndex();
    else
      return -1;
  }

  
  void Ngx_Mesh::HPRefinement (int levels, double parameter, bool setorders,
                               bool ref_level)
  {
    NgLock meshlock (mesh->MajorMutex(), true);
    Refinement & ref = const_cast<Refinement&> (mesh->GetGeometry()->GetRefinement());
    ::netgen::HPRefinement (*mesh, &ref, SPLIT_HP, levels, parameter, setorders, ref_level);
  }

  void Ngx_Mesh::SplitAlfeld ()
  {
    NgLock meshlock (mesh->MajorMutex(), true);
    Refinement & ref = const_cast<Refinement&> (mesh->GetGeometry()->GetRefinement());
    ::netgen::HPRefinement (*mesh, &ref, SPLIT_ALFELD, 1, 1.0/3.0, true, true);
  }

  
int Ngx_Mesh::GetElementOrder (int enr) const
{
  if (mesh->GetDimension() == 3)
    return mesh->VolumeElement(enr).GetOrder();
  else
    return mesh->SurfaceElement(enr).GetOrder();
}

void Ngx_Mesh::GetElementOrders (int enr, int * ox, int * oy, int * oz) const
{
  if (mesh->GetDimension() == 3)
    mesh->VolumeElement(enr).GetOrder(*ox, *oy, *oz);
  else
    mesh->SurfaceElement(enr).GetOrder(*ox, *oy, *oz);
}

void Ngx_Mesh::SetElementOrder (int enr, int order)
{
  if (mesh->GetDimension() == 3)
    return mesh->VolumeElement(enr).SetOrder(order);
  else
    return mesh->SurfaceElement(enr).SetOrder(order);
}

void Ngx_Mesh::SetElementOrders (int enr, int ox, int oy, int oz)
{
  if (mesh->GetDimension() == 3)
    mesh->VolumeElement(enr).SetOrder(ox, oy, oz);
  else
    mesh->SurfaceElement(enr).SetOrder(ox, oy);
}


int Ngx_Mesh::GetSurfaceElementOrder (int enr) const
{
  return mesh->SurfaceElement(enr).GetOrder();
}

int Ngx_Mesh::GetClusterRepVertex (int pi) const
{
  return mesh->GetClusters().GetVertexRepresentant(pi);
}

int Ngx_Mesh::GetClusterRepEdge (int pi) const
{
  return mesh->GetClusters().GetEdgeRepresentant(pi);
}

int Ngx_Mesh::GetClusterRepFace (int pi) const
{
  return mesh->GetClusters().GetFaceRepresentant(pi);
}

int Ngx_Mesh::GetClusterRepElement (int pi) const
{
  return mesh->GetClusters().GetElementRepresentant(pi);
}




//HERBERT: falsche Anzahl von Argumenten
//void Ngx_Mesh::GetSurfaceElementOrders (int enr, int * ox, int * oy, int * oz)
void Ngx_Mesh::GetSurfaceElementOrders (int enr, int * ox, int * oy) const
{
  int d; 
  mesh->SurfaceElement(enr).GetOrder(*ox, *oy, d);
}

void Ngx_Mesh::SetSurfaceElementOrder (int enr, int order)
{
  return mesh->SurfaceElement(enr).SetOrder(order);
}

void Ngx_Mesh::SetSurfaceElementOrders (int enr, int ox, int oy)
{
  mesh->SurfaceElement(enr).SetOrder(ox, oy);
}

  

size_t Ngx_Mesh :: GetGlobalVertexNum (int locnum) const
{
#ifdef PARALLEL  
  return mesh->GetParallelTopology().GetGlobalPNum (locnum+1)-1;
#else
  return locnum;
#endif
}

  
FlatArray<int>  Ngx_Mesh :: GetDistantProcs (int nodetype, int locnum) const
  {
#ifdef PARALLEL
    if (mesh->GetCommunicator().Size() == 1)
      return FlatArray<int>(0,nullptr);
    
    switch (nodetype)
      {
      case 0:
        // return mesh->GetParallelTopology().GetDistantPNums(locnum);
        return mesh->GetParallelTopology().GetDistantProcs(locnum+PointIndex::BASE);
      case 1:
        // return mesh->GetParallelTopology().GetDistantEdgeNums(locnum);
        return mesh->GetParallelTopology().GetDistantEdgeProcs(locnum);
      case 2:
        // return mesh->GetParallelTopology().GetDistantFaceNums(locnum);
        return mesh->GetParallelTopology().GetDistantFaceProcs(locnum);
      default:
	return FlatArray<int>(0, nullptr);
      }
#else
    return FlatArray<int>(0,nullptr);
#endif
  }
}


int link_it_nginterface_v2;

