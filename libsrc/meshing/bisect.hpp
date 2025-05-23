#ifndef NETGEN_BISECT_HPP
#define NETGEN_BISECT_HPP

#include <mydefs.hpp>
#include <general/parthreads.hpp>
#include "basegeom.hpp"
#include "meshclass.hpp"

namespace netgen
{

class BisectionOptions
{
public:
  const char * outfilename;
  const char * mlfilename;
  const char * refinementfilename;
  const char * femcode;
  int maxlevel;
  int usemarkedelements;
  bool refine_hp = false;
  bool refine_p = false;
  bool onlyonce = false;
  NgTaskManager task_manager = &DummyTaskManager;
  NgTracer tracer = &DummyTracer;
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

} // namespace netgen

#endif // NETGEN_BISECT_HPP
