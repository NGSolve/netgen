#ifndef FILE_SPLINESURFACE
#define FILE_SPLINESURFACE


namespace netgen
{
  class SplineSurface  : public  OneSurfacePrimitive
  {
  protected:
    NgArray<GeomPoint<3>> geompoints;
    NgArray<shared_ptr<SplineSeg<3>>> splines;
    NgArray<string> bcnames;
    NgArray<double> maxh;
    shared_ptr<OneSurfacePrimitive> baseprimitive;
    shared_ptr<NgArray<shared_ptr<OneSurfacePrimitive>>> cuts;
    shared_ptr<NgArray<shared_ptr<OneSurfacePrimitive>>> all_cuts;
    
  public:
    SplineSurface(shared_ptr<OneSurfacePrimitive> abaseprimitive, shared_ptr<NgArray<shared_ptr<OneSurfacePrimitive>>> acuts) :
      OneSurfacePrimitive(), baseprimitive(abaseprimitive), cuts(acuts)
    { ; }
    // default constructor for archive
    SplineSurface() {}
    virtual ~SplineSurface() { ; }
    
    const auto & GetSplines() const { return splines; }
    int GetNSplines() const { return splines.Size(); }
    const NgArray<GeomPoint<3>>& GetPoints() const { return geompoints; }
    string GetSplineType(const int i) const { return splines[i]->GetType(); }
    SplineSeg<3> & GetSpline(const int i) { return *splines[i]; }
    const SplineSeg<3> & GetSpline(const int i) const { return *splines[i]; }
    int GetNP() const { return geompoints.Size(); }
    const GeomPoint<3> & GetPoint(int i) const { return geompoints[i]; }
    string GetBCName(int i) const { return bcnames[i]; }
    string GetBCNameOf(Point<3> p1, Point<3> p2) const;
    
    DLL_HEADER void AppendPoint(const Point<3> & p, const double reffac = 1., const bool hpref=false);
    void AppendSegment(shared_ptr<SplineSeg<3>> spline, string & bcname, double amaxh = -1);

    const shared_ptr<NgArray<shared_ptr<OneSurfacePrimitive>>> CreateCuttingSurfaces();
    const shared_ptr<NgArray<shared_ptr<OneSurfacePrimitive>>> GetCuts() const { return all_cuts; }
    const shared_ptr<OneSurfacePrimitive> GetBase() const { return baseprimitive; }
    
    virtual void Project (Point<3> & p3d) const { baseprimitive->Project(p3d); }
    virtual double CalcFunctionValue (const Point<3> & point) const
    { return baseprimitive->CalcFunctionValue (point); }
    virtual void CalcGradient (const Point<3> & point, Vec<3> & grad) const
    { baseprimitive->CalcGradient (point,grad); }
    virtual double HesseNorm () const
    { return baseprimitive->HesseNorm(); }
    virtual Point<3> GetSurfacePoint () const
    { return baseprimitive->GetSurfacePoint(); }
    virtual void CalcSpecialPoints(NgArray<Point<3>> & pts) const
    { baseprimitive->CalcSpecialPoints(pts); }

    virtual INSOLID_TYPE BoxInSolid(const BoxSphere<3> & box) const
    { return baseprimitive->BoxInSolid(box); }

    virtual void DoArchive(Archive& ar)
    {
      ar & geompoints & splines & bcnames & maxh & baseprimitive & cuts & all_cuts;
    }
    
    /*
    virtual void Project (Point<3> & p3d) const;
    virtual double CalcFunctionValue (const Point<3> & point) const;
    virtual void CalcGradient (const Point<3> & point, Vec<3> & grad) const;
    virtual double HesseNorm () const;
    virtual Point<3> GetSurfacePoint () const;
    virtual void CalcSpecialPoints(NgArray<Point<3>> & pts) const;

    virtual INSOLID_TYPE BoxInSolid(const BoxSphere<3> & box) const;

    */
    
    virtual void Print (ostream & str) const;
    
    
  };



}

#endif
