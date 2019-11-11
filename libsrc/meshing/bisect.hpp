#ifndef BISECT
#define BISECT

class BisectionOptions
{
public:
  const char * outfilename;
  const char * mlfilename;
  const char * refinementfilename;
  const char * femcode;
  int maxlevel;
  int usemarkedelements;
  bool refine_hp;
  bool refine_p;
  TaskManager task_manager = &DummyTaskManager;
  Tracer tracer = &DummyTracer;
  DLL_HEADER BisectionOptions ();
};

class ZRefinementOptions
{
public:
  int minref;
  DLL_HEADER ZRefinementOptions();
};



DLL_HEADER extern void BisectTetsCopyMesh (Mesh &, const NetgenGeometry *,
				BisectionOptions & opt);

DLL_HEADER extern void ZRefinement (Mesh &, const class NetgenGeometry *,
			 ZRefinementOptions & opt);





class DLL_HEADER Refinement
{
 const NetgenGeometry& geo;

public:
 Refinement (const NetgenGeometry& ageo) : geo(ageo) {}
 virtual ~Refinement () {}
  
  void Refine (Mesh & mesh) const;
  void Refine (Mesh & mesh);
  void Bisect (Mesh & mesh, class BisectionOptions & opt, NgArray<double> * quality_loss = NULL) const;

  void MakeSecondOrder (Mesh & mesh) const;
  void MakeSecondOrder (Mesh & mesh);

  void ValidateSecondOrder (Mesh & mesh);
  void ValidateRefinedMesh (Mesh & mesh, 
			    NgArray<INDEX_2> & parents);
  
  virtual void LocalizeEdgePoints(Mesh & /* mesh */) const {;}
};

#endif
