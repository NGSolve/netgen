#ifndef FILE_RULER3
#define FILE_RULER3

namespace netgen
{

/**
  3D element generation rule.
 */
class fourpoints { public: RulePointIndex i1, i2, i3, i4; fourpoints() { } };

class vnetrule
{
private:
  /// rule is applicable for quality classes above this value
  int quality;
  /// name of rule
  char * name;
  /// point coordinates in reference position
  Array<Point<3>, RulePointIndex> points;
  /// old and new faces in reference numbering
  Array<RuleElement2d> faces;
  /// additional edges of rule
  Array<twoint> edges;

  /// points of freezone in reference coordinates
  Array<Point<3>> freezone;
  /// points of freezone in reference coordinates if tolcalss to infty
  Array<Point<3>> freezonelimit;
  /// point index, if point equal to mappoint, otherwise 0
  Array<int> freezonepi;
  /// faces of each convex part of freezone
  Array<Array<threeint>*> freefaces;
  /// set of points of each convex part of freezone
  Array<Array<int>*> freesets;
  /// points of transformed freezone
  Array<Point<3>> transfreezone;
  /// edges of each convex part of freezone
  Array<Array<twoint>*> freeedges;

  /// face numbers to be deleted
  Array<int> delfaces;
  /// elements to be generated
  Array<RuleElement> elements;
  /// tolerances for points and faces (used ??)
  Array<double, RulePointIndex> tolerances;
  Array<double> linetolerances;
  /// transformation matrix 
  DenseMatrix oldutonewu;
  /// transformation matrix: deviation old point to dev. freezone
  DenseMatrix * oldutofreezone;
  /** transformation matrix: deviation old point to dev. freezone, 
    quality class to infinity */
  DenseMatrix * oldutofreezonelimit;

  // can be deleted:
  // BaseMatrix *outf, *outfl;

  /**
    a point is outside of convex part of freezone, 
    iff mat * (point, 1) >= 0 for each component (correct ?)
    */
  Array<DenseMatrix*> freefaceinequ;
  /// 
  Array<fourpoints> orientations;
  /**
    flags specified in rule-description file:
    t .. test rule
    */
  Array<char> flags;

  /**
    topological distance of face to base element
    non-connected: > 100  (??) 
    */
  Array<int> fnearness;
  Array<int, RulePointIndex> pnearness;
  int maxpnearness;

  /// number of old points in rule
  int noldp;
  /// number of new poitns in rule
  int noldf;
  /// box containing free-zone
public:  
  // double fzminx, fzmaxx, fzminy, fzmaxy, fzminz, fzmaxz;
  Box3d fzbox;

public:
  
  ///
  vnetrule ();
  ///
  ~vnetrule ();
  ///
  int GetNP () const { return points.Size(); }
  ///
  int GetNF () const { return faces.Size(); }
  ///
  int GetNE () const { return elements.Size(); }
  ///
  int GetNO () const { return orientations.Size(); }
  ///
  int GetNEd () const { return edges.Size(); }
  ///
  int GetNOldP () const { return noldp; }
  ///
  int GetNOldF () const { return noldf; }
  ///
  int GetNDelF () const { return delfaces.Size(); }
  ///
  int GetQuality () const { return quality; }
  ///
  int GetFNearness (int fi) const { return fnearness[fi-1]; }
  ///
  int GetPNearness (RulePointIndex pi) const { return pnearness[pi]; }
  ///
  int GetMaxPNearness () const { return maxpnearness; }


  ///
  const Point<3> & GetPoint (RulePointIndex i) const { return points[i]; }
  ///
  const RuleElement2d & GetFace (int i) const { return faces[i-1]; }
  ///
  const RuleElement & GetElement (int i) const { return elements[i-1]; }
  ///
  const twoint & GetEdge (int i) const { return edges[i-1]; }
  ///
  int GetDelFace (int i) const { return delfaces[i-1]; }
  ///
  int IsDelFace (int fn) const;
  
  ///
  float CalcPointDist (RulePointIndex pi, const Point<3> & p) const;
  ///
  double PointDistFactor (RulePointIndex pi) const
    {
      return tolerances[pi];
    }
  ///
  void SetFreeZoneTransformation (const Vector & allp,
                                  int tolclass);
  ///
  int IsInFreeZone (const Point<3> & p) const;
  /**
    0 not in free-zone
    1 in free-zone
    -1 maybe 
   */
  int IsTriangleInFreeZone (const Point<3> & p1, const Point<3> & p2,
                            const Point<3> & p3, const Array<int> & pi, int newone);
  ///
  int IsQuadInFreeZone (const Point<3> & p1, const Point<3> & p2,
                        const Point<3> & p3, const Point<3> & p4,
                        const Array<int> & pi, int newone);
  ///
  int IsTriangleInFreeSet (const Point<3> & p1, const Point<3> & p2,
                           const Point<3> & p3, int fs, const Array<int> & pi, int newone);

  ///
  int IsQuadInFreeSet (const Point<3> & p1, const Point<3> & p2,
                       const Point<3> & p3, const Point<3> & p4,
                       int fs, const Array<int> & pi, int newone);
  
  ///
  int ConvexFreeZone () const;
  
  /// if t1 and t2 are neighbourtriangles, NTP returns the opposite Point of t1 in t2
  int NeighbourTrianglePoint (const threeint & t1, const threeint & t2) const;
  ///
  const Point<3> & GetTransFreeZone (int i) { return transfreezone[i-1]; }

  ///
  int GetNP (int fn) const
  { return faces[fn-1].GetNP(); }
  ///
  RulePointIndex GetPointNr (int fn, int endp) const
  { return faces[fn-1].PNum(endp); }
  ///
  RulePointIndex GetPointNrMod (int fn, int endp) const
  { return faces[fn-1].PNumMod(endp); }
  ///
  const fourpoints & GetOrientation (int i) { return orientations[i-1]; }

  ///
  int TestFlag (char flag) const;

  ///
  const DenseMatrix & GetOldUToNewU () const { return oldutonewu; }
  //
  //  const DenseMatrix & GetOldUToFreeZone () const { return oldutofreezone; }
  //
  //  const DenseMatrix & GetOldUToFreeZoneLimit () const 
  //    { return oldutofreezonelimit; }
  ///
  const char * Name () const { return name; }
  ///
  void LoadRule (istream & ist);

  ///
  const Array<Point<3>> & GetTransFreeZone () { return transfreezone; }
  ///
  int TestOk () const;

  ///
  friend void TestRules ();
  ///
  //  friend void Plot3DRule (const ROT3D & r, char key);
};

} // namespace netgen
#endif

