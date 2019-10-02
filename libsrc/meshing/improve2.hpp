#ifndef FILE_IMPROVE2
#define FILE_IMPROVE2



///
class MeshOptimize2d
{
  int faceindex = 0;
  int improveedges = 0;
  double metricweight = 0.;
  int writestatus = 1;
  Mesh& mesh;
  const NetgenGeometry& geo;
public:
  ///
  MeshOptimize2d(Mesh& amesh) : mesh(amesh), geo(*mesh.GetGeometry())
  {}
  virtual ~MeshOptimize2d() { ; }
  ///
  void ImproveMesh (const MeshingParameters & mp);
  void ImproveMeshJacobian (const MeshingParameters & mp);
  void ImproveVolumeMesh ();
  void ProjectBoundaryPoints(NgArray<int> & surfaceindex, 
			     const NgArray<Point<3>* > & from, NgArray<Point<3>* > & dest);

  void EdgeSwapping (int usemetric);
  void CombineImprove ();
  void SplitImprove ();

  void GenericImprove ();


  void SetFaceIndex (int fi) { faceindex = fi; }
  void SetImproveEdges (int ie) { improveedges = ie; }
  void SetMetricWeight (double mw) { metricweight = mw; }
  void SetWriteStatus (int ws) { writestatus = ws; }


  /// liefert zu einem 3d-Punkt die geominfo (Dreieck) und liefert 1, wenn erfolgreich, 
  /// 0, wenn nicht (Punkt ausserhalb von chart)
  ///

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


