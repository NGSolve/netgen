#ifndef FILE_IMPROVE2
#define FILE_IMPROVE2


class Neighbour
{
  int nr[3];
  int orient[3];

public:
  Neighbour () { ; }

  void SetNr (int side, int anr) { nr[side] = anr; }
  int GetNr (int side) { return nr[side]; }

  void SetOrientation (int side, int aorient) { orient[side] = aorient; }
  int GetOrientation (int side) { return orient[side]; }
};

///
class MeshOptimize2d
{
  int faceindex;
  int improveedges;
  double metricweight;
  int writestatus;

public:
  ///
  MeshOptimize2d ();
  virtual ~MeshOptimize2d() { ; }
  ///
  void ImproveMesh (Mesh & mesh2d, const MeshingParameters & mp);
  void ImproveMeshJacobian (Mesh & mesh2d, const MeshingParameters & mp);
  void ImproveVolumeMesh (Mesh & mesh);
  void ProjectBoundaryPoints(NgArray<int> & surfaceindex, 
			     const NgArray<Point<3>* > & from, NgArray<Point<3>* > & dest);

  bool EdgeSwapping (Mesh & mesh, const int usemetric, Array<Neighbour> &neighbors, Array<bool> &swapped,
    const SurfaceElementIndex t1, const int edge, const int t, Array<int,PointIndex> &pdef, const bool check_only=false);
  void EdgeSwapping (Mesh & mesh, int usemetric);
  void CombineImprove (Mesh & mesh);
  void SplitImprove (Mesh & mesh);

  void GenericImprove (Mesh & mesh);


  void SetFaceIndex (int fi) { faceindex = fi; }
  void SetImproveEdges (int ie) { improveedges = ie; }
  void SetMetricWeight (double mw) { metricweight = mw; }
  void SetWriteStatus (int ws) { writestatus = ws; }



  ///
  virtual void ProjectPoint (INDEX /* surfind */, Point<3> & /* p */) const { };

  /// project point, use gi as initial value, and compute new gi
  virtual int ProjectPointGI (INDEX surfind, Point<3> & p, PointGeomInfo & gi) const 
  { ProjectPoint (surfind, p); return CalcPointGeomInfo (surfind, gi, p); }

  ///
  virtual void ProjectPoint2 (INDEX /* surfind */, INDEX /* surfind2 */, Point<3> & /* p */) const { };

  /// liefert zu einem 3d-Punkt die geominfo (Dreieck) und liefert 1, wenn erfolgreich, 
  /// 0, wenn nicht (Punkt ausserhalb von chart)
  virtual int CalcPointGeomInfo(PointGeomInfo& gi, const Point<3> & /*p3*/) const
    { gi.trignum = 1; return 1;};

  virtual int CalcPointGeomInfo(int /* surfind */, PointGeomInfo& gi, const Point<3> & p3) const
    { return CalcPointGeomInfo (gi, p3); }

  ///
  virtual void GetNormalVector(INDEX surfind, const Point<3>  & p, PointGeomInfo & gi, Vec<3> & n) const;
  virtual void GetNormalVector(INDEX surfind, const Point<3> & p, Vec<3> & n) const;

  void CheckMeshApproximation (Mesh & mesh);


  ///
  friend class Opti2SurfaceMinFunction;
  ///
  friend class Opti2EdgeMinFunction;
  ///
  friend double Opti2FunctionValueGrad (const Vector & x, Vector & grad);
  ///
  friend double Opti2EdgeFunctionValueGrad (const Vector & x, Vector & grad);



};


extern void CalcTriangleBadness (double x2, double x3, double y3, 
				 double metricweight,
				 double h, double & badness, 
				 double & g1x, double & g1y);




extern double CalcTriangleBadness (const Point<3> & p1, 
				   const Point<3> & p2, 
				   const Point<3> & p3,
				   double metricweight,
				   double h);

extern double CalcTriangleBadness (const Point<3> & p1, 
				   const Point<3> & p2, 
				   const Point<3> & p3,
				   const Vec<3> & n,
				   double metricweight,
				   double h);

#endif


