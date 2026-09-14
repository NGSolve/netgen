#ifndef FILE_NETRULE
#define FILE_NETRULE

namespace netgen
{
///
class netrule
{
private:
  ///
  typedef struct tf 
  { float f1, f2, f3; }   threefloat;
  
  class threeint 
  { 
  public: RulePointIndex i1, i2, i3; 
    threeint() { } 
    threeint(RulePointIndex ai1, RulePointIndex ai2, RulePointIndex ai3) 
    { i1 = ai1; i2 = ai2; i3 = ai3; } 
  };


  ///
  int quality;
  ///
  string name;
  ///
  Array<Point<2>, RulePointIndex> points;
  ///
  NgArray<IVec<2,RulePointIndex>> lines;
  ///
  NgArray<Point<2>> freezone, freezonelimit;
  ///
  NgArray<NgArray<Point<2>>> freezone_i;
  ///
  NgArray<Point<2>> transfreezone;

  ///
  NgArray<int> dellines;
  ///
  NgArray<RuleElement2d> elements;
  ///
  Array<threefloat, RulePointIndex> tolerances;
  ///
  NgArray<threefloat> linetolerances;
  ///
  NgArray<threeint> orientations;
  ///
  DenseMatrix oldutonewu, oldutofreearea, oldutofreearealimit;
  ///
  NgArray<DenseMatrix> oldutofreearea_i;
  ///
  MatrixFixWidth<3> freesetinequ;

  ///
  NgArray<Vec<2>> linevecs;

  ///
  int noldp, noldl;
  ///
  float fzminx, fzmaxx, fzminy, fzmaxy;

  /// topological distance of line to base element
  NgArray<int> lnearness;

public:

  ///
  netrule ();
  ///
  ~netrule();

  ///
  int GetNP () const { return points.Size(); }
  ///
  int GetNL () const { return lines.Size(); }
  ///
  int GetNE () const { return elements.Size(); }
  ///
  int GetNOldP () const { return noldp; }
  ///
  int GetNOldL () const { return noldl; }
  ///
  int GetNDelL () const { return dellines.Size(); }
  ///
  int GetNOrientations () const { return orientations.Size(); }
  ///
  int GetQuality () const { return quality; }
  ///
  int GetLNearness (int li) const { return lnearness[li-1]; }

  ///
  const Point<2>& GetPoint (RulePointIndex i) const { return points[i]; }
  ///
  const IVec<2,RulePointIndex> & GetLine (int i) const { return lines[i-1]; }
  ///
  const RuleElement2d & GetElement (int i) const { return elements[i-1]; }
  ///
  const threeint & GetOrientation (int i) const { return orientations[i-1]; }
  ///
  int GetDelLine (int i) const { return dellines[i-1]; }
  ///
  const NgArray<int> & GetDelLines() const { return dellines; }
  ///
  void GetFreeZone (NgArray<Point<2>> & afreearea);
  ///

  double CalcPointDist (RulePointIndex pi, const Point<2> & p) const
  {
    double dx = p[0] - points[pi][0];
    double dy = p[1] - points[pi][1];
    const threefloat * tfp = &tolerances[pi];
    return tfp->f1 * dx * dx + tfp->f2 * dx * dy + tfp->f3 * dy * dy;
  }

  ///
  float CalcLineError (int li, const Vec<2>& v) const;

  ///
  void SetFreeZoneTransformation (const Vector & u, int tolclass);

  ///
  bool IsInFreeZone (const Point<2> & p) const
  {
    if (p[0] < fzminx || p[0] > fzmaxx ||
	p[1] < fzminy || p[1] > fzmaxy) return 0;

    for (int i = 0; i < transfreezone.Size(); i++)
      {
	if (freesetinequ(i, 0) * p[0] + 
	    freesetinequ(i, 1) * p[1] +
	    freesetinequ(i, 2) > 0) return 0;
      }
    return 1;
  }

  ///
  int IsLineInFreeZone (const Point<2> & p1, const Point<2> & p2) const
  {
    if ( (p1[0] > fzmaxx && p2[0] > fzmaxx) ||
         (p1[0] < fzminx && p2[0] < fzminx) ||
         (p1[1] > fzmaxy && p2[1] > fzmaxy) ||
         (p1[1] < fzminy && p2[1] < fzminy) ) return 0;
    return IsLineInFreeZone2 (p1, p2);
  }
  ///
  int IsLineInFreeZone2 (const Point<2> & p1, const Point<2> & p2) const;
  ///
  int ConvexFreeZone () const;
  ///
  const NgArray<Point<2>> & GetTransFreeZone () { return transfreezone; }

  ///
  RulePointIndex GetPointNr (int ln, int endp) const { return lines[ln-1][endp-1]; }

  ///
  const DenseMatrix & GetOldUToNewU () const { return oldutonewu; }
  ///
  const DenseMatrix & GetOldUToFreeArea () const { return oldutofreearea; }
  ///
  const string & Name () const { return name; }

  ///
  void LoadRule (istream & ist);
};



/** Draws 2D rules.
    Visual testing of 2D meshing rules */
extern void DrawRules ();
} // namespace netgen
#endif

