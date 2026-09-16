#ifndef NETGEN_GEOMSEARCH_HPP
#define NETGEN_GEOMSEARCH_HPP

/**************************************************************************/
/* File:   geomsearch.hh                                                  */
/* Author: Johannes Gerstmayr                                             */
/* Date:   19. Nov. 97                                                    */
/**************************************************************************/

#include "meshtype.hpp"

namespace netgen
{

class FrontPoint3;
class FrontFace;


  /// class for quick access of 3D-elements; class cannot delete elements, but only append
class GeomSearch3d
{

public:
  ///
  GeomSearch3d();
  ///
  virtual ~GeomSearch3d();

  ///
  void Init (Array <FrontPoint3,Front3PointIndex> *pointsi, Array <FrontFace> *facesi);

  ///get elements max extension
  void ElemMaxExt(Point<3>& minp, Point<3>& maxp, const FrontElement2d& elem);
  
  ///get minimum coordinates of two points ->p2
  void MinCoords(const Point<3>& p1, Point<3>& p2);

  ///get minimum coordinates of two points ->p2
  void MaxCoords(const Point<3>& p1, Point<3>& p2);

  ///create a hashtable from an existing array of triangles
  ///sizei = number of pieces in one direction
  void Create();

  ///add new element to Hashtable
  void AddElem(const FrontElement2d& elem, int elemnum);

  ///GetLocal faces in sphere with radius xh and middlepoint p
  void GetLocals(Array<FrontElement2d> & locfaces,  Array<int> & findex,
                 int fstind, const Point<3>& p0, double xh);

private:
  
  Array <FrontFace> *faces; // Pointers to Arrays in Adfront
  Array <FrontPoint3,Front3PointIndex> *points;

  Array <Array <int>*> hashtable;

  Point<3> minext = Point<3>(0,0,0); //extension of Hashdomain
  Point<3> maxext = Point<3>(0,0,0);
  Point<3> maxextreal = Point<3>(0,0,0);
  Vec<3> elemsize = Vec<3>(0,0,0);  //size of one Hash-Element

  IVec<3> size; // size of Hashtable in each direction
  int reset;
  int hashcount;
};
} // namespace netgen
#endif // NETGEN_GEOMSEARCH_HPP
