
#ifndef FILE_IDENTIFY
#define FILE_IDENTIFY

/**************************************************************************/
/* File:   identify.hh                                                    */
/* Author: Joachim Schoeberl                                              */
/* Date:   1. Aug. 99                                                    */
/**************************************************************************/


namespace netgen
{

  /**
     Identify surfaces for periodic b.c. or
     thin domains
  */


  class SpecialPoint;
  class Identification
  {
  protected:
    const CSGeometry & geom;
    // identified faces, index sorted
    INDEX_2_HASHTABLE<int> identfaces;
    int nr;

  public:
    DLL_HEADER Identification (int anr, const CSGeometry & ageom);
    DLL_HEADER virtual ~Identification ();
    DLL_HEADER virtual void Print (ostream & ost) const = 0;
    DLL_HEADER virtual void GetData (ostream & ost) const = 0;

    /// obsolete
    //  virtual void IdentifySpecialPoints (NgArray<class SpecialPoint> & points);

    /// can identify both special points (fixed direction)
    /// (identified points, same tangent)
    virtual int Identifiable (const SpecialPoint & sp1, const SpecialPoint & sp2,
			      const TABLE<int> & specpoint2solid,			  
			      const TABLE<int> & specpoint2surface) const;
    ///
    virtual int Identifiable (const Point<3> & p1, const Point<3> & sp2) const;
    /// is it possible to identify sp1 with some other ?
    virtual int IdentifiableCandidate (const SpecialPoint & sp1) const;
  
    /// are points (if connected) by a short edge (direction anyhow) ?
    virtual int ShortEdge (const SpecialPoint & sp1, const SpecialPoint & sp2) const;

    /// add entries in mesh identification tables
    virtual void IdentifyPoints (class Mesh & mesh);

    /// add entries to identified faces (based on segment infos)
    virtual void IdentifyFaces (class Mesh & mesh);

    /// get point on other surface, add entry in mesh identifications
    virtual PointIndex GetIdentifiedPoint (class Mesh & mesh, PointIndex pi1);

    /// copy surfaces, or fill rectangles
    virtual void BuildSurfaceElements (NgArray<class Segment> & segs,
				       class Mesh & mesh,
				       const Surface * surf);

    /// insert volume elements in thin layers
    virtual void BuildVolumeElements (NgArray<class Element2d> & surfels,
				      class Mesh & mesh);

    /// get list of identified faces
    virtual void GetIdentifiedFaces (NgArray<INDEX_2> & idfaces) const;

    friend ostream & operator<< (ostream & ost, Identification & ident);
  };


  class PeriodicIdentification : public Identification
  {
    const Surface * s1;
    const Surface * s2;
    Transformation<3> trafo; // from s1 to s2
    Transformation<3> inv_trafo; // from s2 to s1
  public:
    PeriodicIdentification (int anr,
			    const CSGeometry & ageom,
			    const Surface * as1,
			    const Surface * as2,
                            Transformation<3> atrafo = Vec<3>(0,0,0));
    virtual ~PeriodicIdentification () override;
    virtual void Print (ostream & ost) const override;
    virtual void GetData (ostream & ost) const override;


    //  virtual void IdentifySpecialPoints (NgArray<class SpecialPoint> & points);
    virtual int Identifiable (const SpecialPoint & sp1, const SpecialPoint & sp2,
			      const TABLE<int> & specpoint2solid,
			      const TABLE<int> & specpoint2surface) const override;

    virtual int Identifiable (const Point<3> & p1, const Point<3> & sp2) const override;
    virtual PointIndex GetIdentifiedPoint (class Mesh & mesh, PointIndex pi1) override;
    virtual void IdentifyPoints (class Mesh & mesh) override;
    virtual void IdentifyFaces (class Mesh & mesh) override;
    virtual void BuildSurfaceElements (NgArray<class Segment> & segs,
				       class Mesh & mesh,
				       const Surface * surf) override;
  };


  ///
  class TopLevelObject;
  class CloseSurfaceIdentification : public Identification
  {
    const Surface * s1;
    const Surface * s2;
    const TopLevelObject * domain;
    ///
    int dom_nr;
    /// number of refinement levels (in Z-refinement)
    int ref_levels;
    /// number of refinement levels for layer next to s1 (in Z-refinement)
    int ref_levels_s1;
    /// number of refinement levels for layer next to s2 (in Z-refinement)
    int ref_levels_s2;
    ///
    double eps_n;
    Array<double> slices;
    /// used only for domain-local identification:
    NgArray<int> domain_surfaces;
    ///
    bool dom_surf_valid;

    ///
    Vec<3> direction;
    ///
    bool usedirection;
  public:
    CloseSurfaceIdentification (int anr, 
				const CSGeometry & ageom,
				const Surface * as1,
				const Surface * as2,
				const TopLevelObject * adomain,
				const Flags & flags);
    virtual ~CloseSurfaceIdentification ();

    virtual void Print (ostream & ost) const;
    virtual void GetData (ostream & ost) const;


    //  virtual void IdentifySpecialPoints (NgArray<class SpecialPoint> & points);
    virtual int Identifiable (const SpecialPoint & sp1, const SpecialPoint & sp2,
			      const TABLE<int> & specpoint2solid,
			      const TABLE<int> & specpoint2surface) const;
    virtual int Identifiable (const Point<3> & p1, const Point<3> & sp2) const;
    virtual int IdentifiableCandidate (const SpecialPoint & sp1) const;
    virtual int ShortEdge (const SpecialPoint & sp1, const SpecialPoint & sp2) const;
    virtual PointIndex GetIdentifiedPoint (class Mesh & mesh, PointIndex pi1);
    const Array<double> & GetSlices () const { return slices; }
    virtual void IdentifyPoints (class Mesh & mesh);
    virtual void IdentifyFaces (class Mesh & mesh);
    virtual void BuildSurfaceElements (NgArray<class Segment> & segs,
				       class Mesh & mesh,
				       const Surface * surf);
    void BuildSurfaceElements2 (NgArray<class Segment> & segs,
				class Mesh & mesh,
				const Surface * surf);

    virtual void BuildVolumeElements (NgArray<class Element2d> & surfels,
				      class Mesh & mesh);

    int RefLevels () const { return ref_levels; }
    int RefLevels1 () const { return ref_levels_s1; }
    int RefLevels2 () const { return ref_levels_s2; }

    bool IsSkewIdentification(void) const {return usedirection;}
    const Vec<3> & GetDirection(void) const {return direction;}

    const Surface & GetSurface1(void) const
    { return *s1;}
    const Surface & GetSurface2(void) const
    { return *s2;}
  };


  class CloseEdgesIdentification : public Identification
  {
    const Surface * facet;
    const Surface * s1;
    const Surface * s2;
  public:
    CloseEdgesIdentification (int anr,
			      const CSGeometry & ageom,
			      const Surface * afacet,
			      const Surface * as1,
			      const Surface * as2);
    virtual ~CloseEdgesIdentification ();
    virtual void Print (ostream & ost) const;
    virtual void GetData (ostream & ost) const;

    //  virtual void IdentifySpecialPoints (NgArray<class SpecialPoint> & points);
    virtual int Identifiable (const SpecialPoint & sp1, const SpecialPoint & sp2,
			      const TABLE<int> & specpoint2solid,
			      const TABLE<int> & specpoint2surface) const;


    virtual void IdentifyPoints (class Mesh & mesh);
    virtual void BuildSurfaceElements (NgArray<class Segment> & segs,
				       class Mesh & mesh,
				       const Surface * surf);
  };

}

#endif
