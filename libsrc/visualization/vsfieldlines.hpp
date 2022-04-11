#ifndef VSFIELDLINES_HPP_INCLUDED
#define VSFIELDLINES_HPP_INCLUDED

namespace netgen
{

class RKStepper
{
private:
  Array<double> c,b;
  TABLE<double> *a;
  int steps;
  int order;

  double tolerance;

  Array<Vec<3>> K;

  int stepcount;

  double h;
  double startt;
  double startt_bak;
  Point<3> startval;
  Point<3> startval_bak;

  bool adaptive;
  int adrun;
  Point<3> valh;

  int notrestarted;

public:

  ~RKStepper();

  RKStepper(int type = 0);

  void SetTolerance(const double tol){tolerance = tol;}

  void StartNextValCalc(const Point<3> & astartval, const double astartt, const double ah, const bool aadaptive = false);

  bool GetNextData(Point<3> & val, double & t, double & ah);

  bool FeedNextF(const Vec<3> & f);
};




class FieldLineCalc
{
private:
  const Mesh & mesh;

  typedef std::function<bool (int elnr, const double *, Vec<3> &)> VectorFunction;

  const VectorFunction & func;
  RKStepper stepper;

  Array<double> values;
  Array<Point<3>> pstart, pend;

  double maxlength;

  int maxpoints;

  int direction;

  Point3d pmin, pmax;
  double rad;

  double critical_value;

  bool randomized;

  double thickness;

public:
  FieldLineCalc(const Mesh & amesh, const VectorFunction & afunc,
		const double rel_length, const int amaxpoints = -1,
		const double rel_thickness = -1, const double rel_tolerance = -1, const int rk_type = 0, const int adirection = 0);

  void SetCriticalValue(const double val) { critical_value = val; }

  void Randomized(void) { randomized = true; }
  void NotRandomized(void) { randomized = false; }

  void Calc(const Point<3> & startpoint, Array<Point<3>> & points, Array<double> & vals, Array<bool> & drawelems, Array<int> & dirstart);

  void GenerateFieldLines(Array<Point<3>> & potential_startpoints, const int numlines);

  const auto & GetPStart() const { return pstart; }
  const auto & GetPEnd() const { return pend; }
  const auto & GetValues() const { return values; }
  const auto GetThickness() const { return thickness; }
};

} // namespace netgen

#endif // VSFIELDLINES_HPP_INCLUDED
