#ifndef NETGEN_CURVEDELEMS_HPP
#define NETGEN_CURVEDELEMS_HPP

/**************************************************************************/
/* File:   curvedelems.hpp                                                */
/* Author: Robert Gaisbauer (first version)                               */
/*         redesign by Joachim Schoeberl                                  */
/* Date:   27. Sep. 02, Feb 2006                                          */
/**************************************************************************/

#include <mydefs.hpp>
#include <general/ngarray.hpp>
#include <gprim/geomobjects.hpp>

#include "meshtype.hpp"
#include "meshclass.hpp"

namespace netgen
{
class Refinement;
class Mesh;

class CurvedElements
{
  const Mesh & mesh;

  NgArray<int> edgeorder;
  NgArray<int> faceorder;

  NgArray<int> edgecoeffsindex;
  NgArray<int> facecoeffsindex;

  NgArray< Vec<3> > edgecoeffs;
  NgArray< Vec<3> > facecoeffs;

  NgArray< double > edgeweight;  // for rational 2nd order splines

  int order;
  bool rational;

  bool ishighorder;

public:
  DLL_HEADER CurvedElements (const Mesh & amesh);
  DLL_HEADER ~CurvedElements();

  // bool IsHighOrder() const { return order > 1; }
  bool IsHighOrder() const { return ishighorder; }

  // void SetHighOrder (int aorder) { order=aorder; }
  void SetIsHighOrder (bool ho) { ishighorder = ho; }
  
  DLL_HEADER void BuildCurvedElements(const Refinement * ref, int aorder, bool arational = false);

  int GetOrder () { return order; }

  void DoArchive(Archive& ar)
  {
    ar & edgeorder & faceorder & edgecoeffsindex & facecoeffsindex & edgecoeffs & facecoeffs
      & edgeweight & order & rational & ishighorder;
  }

  DLL_HEADER bool IsSegmentCurved (SegmentIndex segnr) const;
  DLL_HEADER bool IsSurfaceElementCurved (SurfaceElementIndex sei) const;
  DLL_HEADER bool IsElementCurved (ElementIndex ei) const;
  DLL_HEADER bool IsElementHighOrder (ElementIndex ei) const;


  void CalcSegmentTransformation (double xi, SegmentIndex segnr,
				  Point<3> & x)
  { CalcSegmentTransformation<double> (xi, segnr, &x, NULL); };

  void CalcSegmentTransformation (double xi, SegmentIndex segnr,
				  Vec<3> & dxdxi)
  { CalcSegmentTransformation<double> (xi, segnr, NULL, &dxdxi); };

  void CalcSegmentTransformation (double xi, SegmentIndex segnr,
				  Point<3> & x, Vec<3> & dxdxi)
  { CalcSegmentTransformation<double> (xi, segnr, &x, &dxdxi, NULL); };

  void CalcSegmentTransformation (double xi, SegmentIndex segnr,
				  Point<3> & x, Vec<3> & dxdxi, bool & curved)
  { CalcSegmentTransformation (xi, segnr, &x, &dxdxi, &curved); };



  void CalcSurfaceTransformation (const Point<2> & xi, SurfaceElementIndex elnr,
				  Point<3> & x)
  { CalcSurfaceTransformation (xi, elnr, &x, NULL); };

  void CalcSurfaceTransformation (const Point<2> & xi, SurfaceElementIndex elnr,
				  Mat<3,2> & dxdxi)
  { CalcSurfaceTransformation (xi, elnr, NULL, &dxdxi); };

  void CalcSurfaceTransformation (const Point<2> & xi, SurfaceElementIndex elnr,
				  Point<3> & x, Mat<3,2> & dxdxi)
  { CalcSurfaceTransformation (xi, elnr, &x, &dxdxi, NULL); };

  void CalcSurfaceTransformation (const Point<2> & xi, SurfaceElementIndex elnr,
				  Point<3> & x, Mat<3,2> & dxdxi, bool & curved)
  { CalcSurfaceTransformation (xi, elnr, &x, &dxdxi, &curved); };





  void CalcElementTransformation (const Point<3> & xi, ElementIndex elnr,
				  Point<3> & x)
  { CalcElementTransformation (xi, elnr, &x, NULL); };

  void CalcElementTransformation (const Point<3> & xi, ElementIndex elnr,
				  Mat<3,3> & dxdxi)
  { CalcElementTransformation (xi, elnr, NULL, &dxdxi); };

  void CalcElementTransformation (const Point<3> & xi, ElementIndex elnr,
				  Point<3> & x, Mat<3,3> & dxdxi)
  { CalcElementTransformation (xi, elnr, &x, &dxdxi /* , NULL */ ); };

  void CalcElementTransformation (const Point<3> & xi, ElementIndex elnr,
				  Point<3> & x, Mat<3,3> & dxdxi,
                                  void * buffer, bool valid)
  { CalcElementTransformation (xi, elnr, &x, &dxdxi, /* NULL, */ buffer, valid ); };

  // void CalcElementTransformation (const Point<3> & xi, ElementIndex elnr,
  // 				  Point<3> & x, Mat<3,3> & dxdxi) // , bool & curved)
  //   { CalcElementTransformation (xi, elnr, &x, &dxdxi /* , &curved * ); }


  /*
  void CalcMultiPointSegmentTransformation (NgArray<double> * xi, SegmentIndex segnr,
					    NgArray<Point<3> > * x,
					    NgArray<Vec<3> > * dxdxi);
  */
  
  template <int DIM_SPACE, typename T>
  void CalcMultiPointSegmentTransformation (SegmentIndex elnr, int n,
                                            const T * xi, size_t sxi,
                                            T * x, size_t sx,
                                            T * dxdxi, size_t sdxdxi);

  DLL_HEADER void CalcMultiPointSurfaceTransformation (NgArray< Point<2> > * xi, SurfaceElementIndex elnr,
					    NgArray< Point<3> > * x,
					    NgArray< Mat<3,2> > * dxdxi);

  template <int DIM_SPACE, typename T>
  void CalcMultiPointSurfaceTransformation (SurfaceElementIndex elnr, int n,
                                            const T * xi, size_t sxi,
                                            T * x, size_t sx,
                                            T * dxdxi, size_t sdxdxi);

  DLL_HEADER void CalcMultiPointElementTransformation (NgArray< Point<3> > * xi, ElementIndex elnr,
					    NgArray< Point<3> > * x,
					    NgArray< Mat<3,3> > * dxdxi);

  template <typename T>
  void CalcMultiPointElementTransformation (ElementIndex elnr, int n,
                                            const T * xi, size_t sxi,
                                            T * x, size_t sx,
                                            T * dxdxi, size_t sdxdxi);




private:

  template <typename T>
  DLL_HEADER void CalcSegmentTransformation (const T & xi, SegmentIndex segnr,
				  Point<3,T> * x = NULL, Vec<3,T> * dxdxi = NULL, bool * curved = NULL);

  DLL_HEADER void CalcSurfaceTransformation (Point<2> xi, SurfaceElementIndex elnr,
				  Point<3> * x = NULL, Mat<3,2> * dxdxi = NULL, bool * curved = NULL);

  DLL_HEADER void CalcElementTransformation (Point<3> xi, ElementIndex elnr,
				  Point<3> * x = NULL, Mat<3,3> * dxdxi = NULL, // bool * curved = NULL,
                                  void * buffer = NULL, bool valid = 0);






  class SegmentInfo
  {
  public:
    SegmentIndex elnr;
    int order;
    int nv;
    int ndof;
    int edgenr;
  };

  template <typename T>
  void CalcElementShapes (SegmentInfo &  elnr, T xi, TFlatVector<T> shapes) const;
  void GetCoefficients (SegmentInfo & elnr, NgArray<Vec<3> > & coefs) const;
  template <typename T>
  void CalcElementDShapes (SegmentInfo & elnr, T xi, TFlatVector<T> dshapes) const;


  class ElementInfo
  {
  public:
    ElementIndex elnr;
    int order;
    int nv;
    int ndof;
    int nedges;
    int nfaces;
    int edgenrs[12];
    int facenrs[6];
    Mat<3> hdxdxi;
    Vec<3> hcoefs[10]; // enough for second order tets

    void SetEdges (FlatArray<T_EDGE> edges)
    {
      nedges = edges.Size();
      for (int i = 0; i < edges.Size(); i++)
        edgenrs[i] = edges[i];
    }
    
    auto GetEdges() const
    { return FlatArray(nedges, edgenrs); }

    void SetFaces (FlatArray<T_FACE> faces)
    {
      nfaces = faces.Size();
      for (int i = 0; i < faces.Size(); i++)
        facenrs[i] = faces[i];
    }

    auto GetFaces() const
    { return FlatArray(nfaces, facenrs); }
  };

  template <typename T>
  void CalcElementShapes (ElementInfo & info, Point<3,T> xi, TFlatVector<T> shapes) const;
  void GetCoefficients (ElementInfo & info, Vec<3> * coefs) const;
  template <typename T>  
  void CalcElementDShapes (ElementInfo & info, const Point<3,T> xi, MatrixFixWidth<3,T> & dshapes) const;

  template <typename T>
  bool EvaluateMapping (ElementInfo & info, const Point<3,T> xi, Point<3,T> & x, Mat<3,3,T> & jac) const;  
  
  class SurfaceElementInfo
  {
  public:
    SurfaceElementIndex elnr;
    int order;
    int nv;
    int ndof;
    NgArrayMem<int,4> edgenrs;
    int facenr;

    void SetEdges (FlatArray<T_EDGE> edges)
    {
      edgenrs.SetSize(edges.Size());
      for (int i = 0; i < edges.Size(); i++)
        edgenrs[i] = edges[i];
    }
    
  };

  template <typename T>
  void CalcElementShapes (SurfaceElementInfo & elinfo, const Point<2,T> xi, TFlatVector<T> shapes) const;
  template <int DIM_SPACE>
  void GetCoefficients (SurfaceElementInfo & elinfo, NgArray<Vec<DIM_SPACE> > & coefs) const;
  template <typename T>
  void CalcElementDShapes (SurfaceElementInfo & elinfo, const Point<2,T> xi, MatrixFixWidth<2,T> & dshapes) const;

  template <int DIM_SPACE, typename T>
  bool EvaluateMapping (SurfaceElementInfo & info, const Point<2,T> xi, Point<DIM_SPACE,T> & x, Mat<DIM_SPACE,2,T> & jac) const;  
};

} //namespace netgen
#endif // NETGEN_CURVEDELEMS_HPP
