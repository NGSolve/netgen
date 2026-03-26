#ifndef FILE_IMPROVE3
#define FILE_IMPROVE3

namespace netgen
{


///
class MeshOptimize3d
{
  Mesh & mesh;
  const MeshingParameters & mp;
  OPTIMIZEGOAL goal = OPT_QUALITY;
  double min_badness = 0;

  bool HasBadElement(FlatArray<ElementIndex> els);
  bool HasIllegalElement(FlatArray<ElementIndex> els);
  bool NeedsOptimization(FlatArray<ElementIndex> els);

public:

  MeshOptimize3d (Mesh & m, const MeshingParameters & amp, OPTIMIZEGOAL agoal = OPT_QUALITY) :
      mesh(m), mp(amp), goal(agoal) { ; }

  void SetGoal(OPTIMIZEGOAL agoal) { goal = agoal; }
  void SetMinBadness(double badness) { min_badness = badness; }

  tuple<double, double, int> UpdateBadness();

  double CombineImproveEdge (
            Table<ElementIndex, PointIndex> & elements_of_point,
            PointIndex pi0, PointIndex pi1,
            FlatArray<bool, PointIndex> is_point_removed, bool check_only=false);

  void CombineImprove ();

  void SplitImprove ();
  double SplitImproveEdge (Table<ElementIndex,PointIndex> & elementsonnode, NgArray<PointIndices<3>> &locfaces, double badmax, PointIndex pi1, PointIndex pi2, PointIndex ptmp, bool check_only=false);

  void SplitImprove2 ();
  double SplitImprove2Element (ElementIndex ei, const Table<ElementIndex, PointIndex> & elements_of_point, bool check_only);
  

  double SwapImproveEdge (const TBitArray<ElementIndex> * working_elements, Table<ElementIndex,PointIndex> & elementsonnode, INDEX_3_HASHTABLE<int> & faces, PointIndex pi1, PointIndex pi2, bool check_only=false);
  void SwapImprove (const TBitArray<ElementIndex> * working_elements = NULL);
  void SwapImproveSurface (const TBitArray<ElementIndex> * working_elements = NULL,
			   const NgArray< idmap_type* > * idmaps = NULL);
  void SwapImprove2 (bool conform_segments = false);
  double SwapImprove2 (ElementIndex eli1, int face, Table<ElementIndex, PointIndex> & elementsonnode, DynamicTable<SurfaceElementIndex, PointIndex> & belementsonnode, bool conform_segments, bool check_only=false );

  void ImproveMesh() { mesh.ImproveMesh(mp, goal); }

  double 
  CalcBad (const Mesh::T_POINTS & points, const Element & elem, double h)
  {
    if (elem.GetType() == TET)
      return CalcTetBadness (points[elem[0]], points[elem[1]],  
			     points[elem[2]], points[elem[3]], h, mp);  
    return 0;
  }


  double GetLegalPenalty()
  {
    return goal == OPT_LEGAL ? 1e15 : 1e6;
  }

};


inline double 
CalcBad (const Mesh::T_POINTS & points, const Element & elem, double h, const MeshingParameters & mp)
{
  if (elem.GetType() == TET)
    return CalcTetBadness (points[elem[0]], points[elem[1]],  
			   points[elem[2]], points[elem[3]], h, mp);  
  return 0;
}



extern int WrongOrientation (const Mesh::T_POINTS & points, const Element & el);


/* Functional depending of inner point inside triangular surface */


class MinFunctionSum : public MinFunction
{
protected:
  NgArray<MinFunction*> functions;
 
public:
  
  virtual double Func (const Vector & x) const;
  virtual void Grad (const Vector & x, Vector & g) const;
  virtual double FuncGrad (const Vector & x, Vector & g) const;
  virtual double FuncDeriv (const Vector & x, const Vector & dir, double & deriv) const;
  virtual double GradStopping (const Vector & x) const;

  void AddFunction(MinFunction & fun);
  
  const MinFunction & Function(int i) const;
  MinFunction & Function(int i);  
};
  

class PointFunction1 : public MinFunction
{
  Mesh::T_POINTS & points;
  const NgArray<PointIndices<3>> & faces;
  const MeshingParameters & mp;
  double h;
public:
  PointFunction1 (Mesh::T_POINTS & apoints, 
		  const NgArray<PointIndices<3>> & afaces,
		  const MeshingParameters & amp,
		  double ah);
  
  virtual double Func (const Vector & x) const;
  virtual double FuncDeriv (const Vector & x, const Vector & dir, double & deriv) const;
  virtual double FuncGrad (const Vector & x, Vector & g) const;
  virtual double GradStopping (const Vector & x) const;
};

class JacobianPointFunction : public MinFunction
{
public:
  Mesh::T_POINTS & points;
  const Array<Element, ElementIndex> & elements;
  TABLE<INDEX> elementsonpoint;
  PointIndex actpind;

  bool onplane;
  Vec<3> nv;
  
public:
  JacobianPointFunction (Mesh::T_POINTS & apoints, 
			 const Array<Element, ElementIndex> & aelements);
  virtual ~JacobianPointFunction () { ; }
  virtual void SetPointIndex (PointIndex aactpind);
  virtual double Func (const Vector & x) const;
  virtual double FuncGrad (const Vector & x, Vector & g) const;
  virtual double FuncDeriv (const Vector & x, const Vector & dir, double & deriv) const;

  inline void SetNV(const Vec<3> & anv) {nv = anv; onplane = true;}
  inline void UnSetNV(void) {onplane = false;}
};

} // namespace netgen
#endif
