#ifndef FILE_SOLID
#define FILE_SOLID

/**************************************************************************/
/* File:   solid.hh                                                       */
/* Author: Joachim Schoeberl                                              */
/* Date:   1. Dez. 95                                                     */
/**************************************************************************/

#include <functional>

namespace netgen 
{

  
  /*
    
  Constructive Solid Model (csg)
  
  */
  



  class Solid;
  
  class SolidIterator
  {
  public:
    SolidIterator () { ; }
    virtual ~SolidIterator () { ; }
    virtual void Do (Solid * sol) = 0;
  };


  inline INSOLID_TYPE Intersection (INSOLID_TYPE ina, INSOLID_TYPE inb)
  {
    if (ina == IS_INSIDE && inb == IS_INSIDE) return IS_INSIDE;
    if (ina == IS_OUTSIDE || inb == IS_OUTSIDE) return IS_OUTSIDE;
    return DOES_INTERSECT;
  }
  
  inline INSOLID_TYPE Union (INSOLID_TYPE ina, INSOLID_TYPE inb)
  {
    if (ina == IS_INSIDE || inb == IS_INSIDE) return IS_INSIDE;
    if (ina == IS_OUTSIDE && inb == IS_OUTSIDE) return IS_OUTSIDE;
    return DOES_INTERSECT;
  }

  inline INSOLID_TYPE Complement (INSOLID_TYPE in)
  {
    if (in == IS_INSIDE) return IS_OUTSIDE;
    if (in == IS_OUTSIDE) return IS_INSIDE;
    return DOES_INTERSECT;
  }

  class Solid
  {
  public:
  
    typedef enum optyp1 { TERM, TERM_REF, SECTION, UNION, SUB, ROOT /*, DUMMY */ } optyp;
  
  private:
    char * name;
    Primitive * prim;
    Solid * s1, * s2;
  
    optyp op;
    bool visited;
    double maxh;
    int num_surfs;

    // static int cntnames;

  public:
    Solid (Primitive * aprim);
    Solid (optyp aop, Solid * as1, Solid * as2 = NULL);
    // default constructor for archive
    Solid () {}
    ~Solid ();

    void DoArchive(Archive& archive)
    {
      archive & name & prim & s1 & s2 & visited & maxh & num_surfs;
      if(archive.Output())
        archive << int(op);
      else
        {
          int iop;
          archive & iop;
          op = optyp(iop);
        }
    }
    const char * Name () const { return name; }
    void SetName (const char * aname);

    Solid * Copy (class CSGeometry & geom) const;
    void Transform (Transformation<3> & trans);

  
    void IterateSolid (SolidIterator & it, bool only_once = 0);

  
    void Boundaries (const Point<3> & p, NgArray<int> & bounds) const;
    int NumPrimitives () const;
    void GetSurfaceIndices (NgArray<int> & surfind) const;
    void GetSurfaceIndices (IndexSet & iset) const;

    void GetTangentialSurfaceIndices (const Point<3> & p, NgArray<int> & surfids, double eps) const;
    void GetTangentialSurfaceIndices2 (const Point<3> & p, const Vec<3> & v, NgArray<int> & surfids, double eps) const;
    void GetTangentialSurfaceIndices3 (const Point<3> & p, const Vec<3> & v, const Vec<3> & v2, NgArray<int> & surfids, double eps) const;

    void ForEachSurface (const std::function<void(Surface*,bool)> & lambda, bool inv = false) const;

    Primitive * GetPrimitive ()
    { return (op == TERM || op == TERM_REF) ? prim : NULL; }
    const Primitive * GetPrimitive () const
    { return (op == TERM || op == TERM_REF) ? prim : NULL; }

    Solid * S1() { return s1; }
    Solid * S2() { return s2; }

    // geometric tests

    INSOLID_TYPE PointInSolid (const Point<3> & p, double eps) const;
    INSOLID_TYPE VecInSolid (const Point<3> & p, const Vec<3> & v, double eps) const;

    // checks if lim s->0 lim t->0  p + t(v1 + s v2) in solid
    INSOLID_TYPE VecInSolid2 (const Point<3> & p, const Vec<3> & v1,
                              const Vec<3> & v2, double eps) const;

    
    bool IsIn (const Point<3> & p, double eps = 1e-6) const;
    bool IsStrictIn (const Point<3> & p, double eps = 1e-6) const;
    bool VectorIn (const Point<3> & p, const Vec<3> & v, double eps = 1e-6) const;
    bool VectorStrictIn (const Point<3> & p, const Vec<3> & v, double eps = 1e-6) const;
  
    bool VectorIn2 (const Point<3> & p, const Vec<3> & v1, const Vec<3> & v2,
		    double eps) const;
    /*
    bool VectorIn2Rec (const Point<3> & p, const Vec<3> & v1, const Vec<3> & v2,
		       double eps) const;
    */
    bool VectorStrictIn2 (const Point<3> & p, const Vec<3> & v1, const Vec<3> & v2,
                          double eps) const;

    /// compute localization in point p
    unique_ptr<Solid> TangentialSolid (const Point<3> & p, NgArray<int> & surfids, double eps) const;

    /// compute localization in point p tangential to vector t
    unique_ptr<Solid> TangentialSolid2 (const Point<3> & p, const Vec<3> & t,
                                        NgArray<int> & surfids, double eps) const;

    /** compute localization in point p, with second order approximation to edge
	p + s t + s*s/2 t2 **/
    unique_ptr<Solid> TangentialSolid3 (const Point<3> & p, const Vec<3> & t, const Vec<3> & t2, 
                                        NgArray<int> & surfids, double eps) const;



    /** tangential solid, which follows the edge
	p + s t + s*s/2 t2
	with second order, and the neighbouring face
	p + s t + s*s/2 t2 + r m
	with first order
    **/
    unique_ptr<Solid> TangentialEdgeSolid (const Point<3> & p, const Vec<3> & t, const Vec<3> & t2, 
                                           const Vec<3> & m, 
                                           NgArray<int> & surfids, double eps) const;


    void CalcOnePrimitiveSpecialPoints (const Box<3> & box, NgArray<Point<3> > & pts) const;

    ///
    int Edge (const Point<3> & p, const Vec<3> & v, double eps) const;
    ///
    int OnFace (const Point<3> & p, const Vec<3> & v, double eps) const;
    ///
    void Print (ostream & str) const;
    ///
    void CalcSurfaceInverse ();
    ///
    Solid * GetReducedSolid (const BoxSphere<3> & box) const;
  

    void SetMaxH (double amaxh)
    { maxh = amaxh; }
    double GetMaxH () const
    { return maxh; }

    void GetSolidData (ostream & ost, int first = 1) const;
    static Solid * CreateSolid (istream & ist, const SymbolTable<Solid*> & solids);


    static shared_ptr<BlockAllocator> ball;
    void * operator new(size_t /* s */) 
    {
      return ball->Alloc();
    }

    void operator delete (void * p)
    {
      ball->Free (p);
    }


  protected:
    ///

    void RecBoundaries (const Point<3> & p, NgArray<int> & bounds, 
			int & in, int & strin) const;
    ///
    void RecTangentialSolid (const Point<3> & p, Solid *& tansol, NgArray<int> & surfids, 
                             bool & in, bool & strin, double eps) const;

    void RecTangentialSolid2 (const Point<3> & p, const Vec<3> & vec, 
			      Solid *& tansol, NgArray<int> & surfids, 
			      bool & in, bool & strin, double eps) const;
    ///
    void RecTangentialSolid3 (const Point<3> & p, const Vec<3> & vec,const Vec<3> & vec2, 
			      Solid *& tansol, NgArray<int> & surfids, 
			      bool & in, bool & strin, double eps) const;
    ///
    void RecTangentialEdgeSolid (const Point<3> & p, const Vec<3> & t, const Vec<3> & t2, 
				 const Vec<3> & m, 
				 Solid *& tansol, NgArray<int> & surfids, 
				 bool & in, bool & strin, double eps) const;

    ///
    void RecEdge (const Point<3> & p, const Vec<3> & v,
		  bool & in, bool & strin, int & faces, double eps) const;
    ///
    void CalcSurfaceInverseRec (int inv);
    ///
    Solid * RecGetReducedSolid (const BoxSphere<3> & box, INSOLID_TYPE & in) const;
    ///
    void RecGetSurfaceIndices (NgArray<int> & surfind) const;
    void RecGetTangentialSurfaceIndices (const Point<3> & p, NgArray<int> & surfids, double eps) const;
    void RecGetTangentialSurfaceIndices2 (const Point<3> & p, const Vec<3> & v, NgArray<int> & surfids, double eps) const;
    void RecGetTangentialSurfaceIndices3 (const Point<3> & p, const Vec<3> & v, const Vec<3> & v2, 
					  NgArray<int> & surfids, double eps) const;
    void RecGetTangentialEdgeSurfaceIndices (const Point<3> & p, const Vec<3> & v, const Vec<3> & v2, const Vec<3> & m,
					     NgArray<int> & surfids, double eps) const;
    void RecGetSurfaceIndices (IndexSet & iset) const;

    void RecCalcOnePrimitiveSpecialPoints (NgArray<Point<3> > & pts) const;

    friend class SolidIterator;
    friend class ClearVisitedIt;
    friend class RemoveDummyIterator;
    friend class CSGeometry;
  };


  inline ostream & operator<< (ostream & ost, const Solid & sol)
  {
    sol.Print (ost);
    return ost;
  }






  class ReducePrimitiveIterator : public SolidIterator
  {
    BoxSphere<3> box;
  public:
    ReducePrimitiveIterator (const BoxSphere<3> & abox)
      : SolidIterator(), box(abox) { ; }
    virtual ~ReducePrimitiveIterator () { ; }
    virtual void Do (Solid * sol)
    {
      if (sol -> GetPrimitive())
	sol -> GetPrimitive() -> Reduce (box);
    }
  };


  class UnReducePrimitiveIterator : public SolidIterator
  {
  public:
    UnReducePrimitiveIterator () { ; }
    virtual ~UnReducePrimitiveIterator () { ; }
    virtual void Do (Solid * sol)
    {
      if (sol -> GetPrimitive())
	sol -> GetPrimitive() -> UnReduce ();
    }
  };

}

#endif
