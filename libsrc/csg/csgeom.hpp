#ifndef FILE_CSGEOM
#define FILE_CSGEOM

/**************************************************************************/
/* File:   csgeom.hh                                                      */
/* Author: Joachim Schoeberl                                              */
/* Date:   27. Nov. 97                                                    */
/**************************************************************************/

namespace netgen
{

  /**
     Constructive Solid Geometry
  */


  class TriangleApproximation;
  class TATriangle;

  /**
     A top level object is an entity to be meshed.
     I can be either a solid, or one surface patch of a solid.
  */
  class DLL_HEADER TopLevelObject
  {
    Solid * solid;
    Surface * surface;

    double red, blue, green;
    bool visible, transp;
    double maxh;
    string material;
    int layer;
    int bc;     // for surface patches, only
    string bcname;

  public:
    TopLevelObject (Solid * asolid,
		    Surface * asurface = NULL);
    // default constructor for archive
    TopLevelObject() {}

    void DoArchive(Archive& archive)
    {
      archive & solid & surface & red & blue & green & visible & transp & maxh
        & material & layer & bc & bcname;
    }
    const Solid * GetSolid() const { return solid; }
    Solid * GetSolid() { return solid; }

    const Surface * GetSurface () const { return surface; }
    Surface  * GetSurface () { return surface; }

    void GetData (ostream & ost);
    void SetData (istream & ist);

    void SetMaxH (double amaxh) { maxh = amaxh; } 
    double GetMaxH () const { return maxh; }

    void SetRGB (double ared, double agreen, double ablue)
    {
      red = ared;
      green = agreen;
      blue = ablue;
    }

    double GetRed () const { return red; }
    double GetGreen () const { return green; }
    double GetBlue () const { return blue; }

    void SetTransparent (bool atransp) 
    { transp = atransp; }
    bool GetTransparent () const { return transp; }

    void SetVisible (bool avisible)
    { visible = avisible; }
    bool GetVisible () const { return visible; }

    const string GetMaterial () const { return material; }
    void SetMaterial (const string & mat) { material = mat; }

    int GetLayer () const { return layer; }
    void SetLayer (int alayer) { layer = alayer; }

    void SetBCProp (int abc) { bc = abc; }
    int GetBCProp () const { return bc; }

    void SetBCName (string abc) { bcname = abc; }
    const string GetBCName () const { return bcname; }
  };





  /**
     CSGeometry has the whole geometric information
  */
  class DLL_HEADER CSGeometry : public NetgenGeometry
  {
  private:
    /// all surfaces
    SymbolTable<Surface*> surfaces;

  public:
    /// primitive of surface
    NgArray<const Primitive*> surf2prim;

  private:
    NgArray<Surface*> delete_them;

    /// all named solids
    SymbolTable<Solid*> solids;

    /// all 2d splinecurves
    SymbolTable<shared_ptr<SplineGeometry<2>>> splinecurves2d;
    /// all 3d splinecurves
    SymbolTable<shared_ptr<SplineGeometry<3>>> splinecurves3d;

    /// all top level objects: solids and surfaces
    NgArray<TopLevelObject*> toplevelobjects;

  public:
    /// additional points specified by user
    class UserPoint : public Point<3>
    {
      int index;
      string name;
    public:
      UserPoint() = default;
      UserPoint (Point<3> p, int _index) : Point<3>(p), index(_index) { ; }
      UserPoint (Point<3> p, const string & _name) : Point<3>(p), index(-1), name(_name) { ; } 
      int GetIndex() const { return index; }
      const string & GetName() const { return name; } 
      void DoArchive(Archive& archive)
      {
        archive & index & name;
        Point<3>::DoArchive(archive);
      }
    };
    
  private:
    // NgArray<Point<3> > userpoints;
    NgArray<UserPoint> userpoints;
    NgArray<double> userpoints_ref_factor;

    mutable NgArray<Point<3> > identpoints;

    /// triangular approximation of top level objects
    NgArray<TriangleApproximation*> triapprox;

    /// increment, if geometry is changed
    static int changeval;
  
    /// bounding box of geometry
    Box<3> boundingbox;

    /// bounding box, if not set by input file
    static Box<3> default_boundingbox;

    /// identic surfaces are stored by pair of indizes, val = inverse
    INDEX_2_HASHTABLE<int> identicsurfaces;
    NgArray<int> isidenticto;
    /// identification of boundaries (periodic, thin domains, ...)

    double ideps;

    /// filename of inputfile
    string filename;

    /// store splinesurfaces, such that added ones do not get deleted before geometry does
    NgArray<shared_ptr<SplineSurface>> spline_surfaces;

    shared_ptr<BlockAllocator> solid_ball = Solid::ball;
    
  public:
    CSGeometry ();
    CSGeometry (const string & afilename);
    virtual ~CSGeometry ();

    void Clean ();

    virtual void Save (const filesystem::path & filename) const override;
    void Save (ostream & ost) const;
    void Load (istream & ist);

    void SaveSurfaces (ostream & out) const;
    void LoadSurfaces (istream & in);

    virtual void SaveToMeshFile (ostream & ost) const override;

    PointGeomInfo ProjectPoint(INDEX surfind, Point<3> & p) const override;
    bool ProjectPointGI (int surfind, Point<3> & p, PointGeomInfo & gi) const override;
    void ProjectPointEdge(INDEX surfind, INDEX surfind2, Point<3> & p,
                          EdgePointGeomInfo* gi = nullptr) const override;
    Vec<3> GetNormal(int surfind, const Point<3> & p, const PointGeomInfo* gi = nullptr) const override;

    void PointBetween(const Point<3> & p1, const Point<3> & p2,
                      double secpoint, int surfi,
                      const PointGeomInfo & gi1,
                      const PointGeomInfo & gi2,
                      Point<3> & newp, PointGeomInfo & newgi) const override;

    void PointBetweenEdge(const Point<3> & p1, const Point<3> & p2, double secpoint,
                      int surfi1, int surfi2,
                      const EdgePointGeomInfo & ap1,
                      const EdgePointGeomInfo & ap2,
                      Point<3> & newp, EdgePointGeomInfo & newgi) const override;

    Vec<3> GetTangent (const Point<3> & p, int surfi1, int surfi2,
                       const EdgePointGeomInfo & ap1) const override;

    int GetChangeVal() { return changeval; }
    void Change() { changeval++; }

    void AddSurface (Surface * surf);
    void AddSurface (char * name, Surface * surf);
    void AddSurfaces (Primitive * prim);

    int GetNSurf () const { return surfaces.Size(); }
    const Surface * GetSurface (const char * name) const;
    const Surface * GetSurface (int i) const
    { return surfaces[i]; }

    void SetSolid (const char * name, Solid * sol);
    const Solid * GetSolid (const char * name) const;
    const Solid * GetSolid (const string & name) const;
    int GetNSolids () const { return solids.Size(); }
    const Solid * GetSolid (int i) const { return solids[i]; }
    const SymbolTable<Solid*> & GetSolids () const { return solids; }


    void SetSplineCurve (const char * name, shared_ptr<SplineGeometry<2>> spl);
    void SetSplineCurve (const char * name, shared_ptr<SplineGeometry<3>> spl);
    shared_ptr<SplineGeometry<2>> GetSplineCurve2d (const string & name) const;
    shared_ptr<SplineGeometry<3>> GetSplineCurve3d (const string & name) const;

    void DoArchive(Archive& archive) override;
    

    void SetFlags (const char * solidname, const Flags & flags);


    int GetNTopLevelObjects () const
    { return toplevelobjects.Size(); }
    int SetTopLevelObject (Solid * sol, Surface * surf = NULL);
    void GetTopLevelObject (int nr, Solid *& sol, Surface *& surf)
    {
      sol = toplevelobjects[nr]->GetSolid();
      surf = toplevelobjects[nr]->GetSurface();
    }
    void GetTopLevelObject (int nr, const Solid *& sol, const Surface *& surf) const
    {
      sol = toplevelobjects[nr]->GetSolid();
      surf = toplevelobjects[nr]->GetSurface();
    }

    TopLevelObject * GetTopLevelObject (const Solid * sol, const Surface * surf = NULL);
    TopLevelObject * GetTopLevelObject (int nr) const
    { return toplevelobjects[nr]; }
    // const TopLevelObject * GetTopLevelObject (int nr) const
    // { return toplevelobjects[nr]; }
    void RemoveTopLevelObject (Solid * sol, Surface * surf = NULL); 


    void AddUserPoint (const Point<3> & p, double ref_factor = 0)      
    { userpoints.Append (UserPoint(p,userpoints.Size()+1)); userpoints_ref_factor.Append (ref_factor); }
    void AddUserPoint (const UserPoint up, double ref_factor = 0)
    { userpoints.Append (up); userpoints_ref_factor.Append (ref_factor); }
    int GetNUserPoints () const
    { return userpoints.Size(); }
    const UserPoint & GetUserPoint (int nr) const
    { return userpoints[nr]; }
    double GetUserPointRefFactor (int nr) const
    { return userpoints_ref_factor[nr]; }
  
    void AddIdentPoint (const Point<3> & p) const
    { identpoints.Append(p);}
    int GetNIdentPoints (void) const
    { return identpoints.Size();}
    const Point<3> & GetIdentPoint(int nr) const
    { return identpoints[nr]; }
    void DeleteIdentPoints(void) const
    { identpoints.DeleteAll();}


    // quick implementations:
    NgArray<SingularFace*> singfaces;
    NgArray<SingularEdge*> singedges;
    NgArray<SingularPoint*> singpoints;
    NgArray<Identification*> identifications;

    int GetNIdentifications (void) const { return identifications.Size(); }
    void AddIdentification (Identification * ident);


    ///
    void CalcTriangleApproximation(double detail, double facets);

    ///
    void FindIdenticSurfaces (double eps);
    ///
    void GetSurfaceIndices (const Solid * sol, 
			    const BoxSphere<3> & box, 
			    NgArray<int> & locsurf) const;
    ///
    void GetIndependentSurfaceIndices (const Solid * sol, 
				       const BoxSphere<3> & box, 
				       NgArray<int> & locsurf) const;
    ///
    /*
    void GetIndependentSurfaceIndices (const Solid * sol, 
				       const Point<3> & p, Vec<3> & v,
				       NgArray<int> & locsurf) const;
    */
    ///
    void GetIndependentSurfaceIndices (NgArray<int> & locsurf) const;

    ///
    int GetSurfaceClassRepresentant (int si) const
    { return isidenticto[si]; }

    ///
    const TriangleApproximation * GetTriApprox (int msnr)
    {
      if (msnr < triapprox.Size())
	return triapprox[msnr];
      return 0;
    }
  

    void IterateAllSolids (SolidIterator & it, bool only_once = false) const;

    void RefineTriangleApprox (Solid * locsol, 
			       int surfind,
			       const BoxSphere<3> & box, 
			       double detail,
			       const TATriangle & tria, 
			       TriangleApproximation & tams,
			       IndexSet & iset,
			       int level);

    const Box<3> & BoundingBox () const { return boundingbox; }

    void SetBoundingBox (const Box<3> & abox)
    {
      boundingbox = abox;
    }


    static void SetDefaultBoundingBox (const Box<3> & abox)
    {
      default_boundingbox = abox;
    }

    double MaxSize () const;

    void SetIdEps(double eps){ideps = eps;}
    double GetIdEps(void) const {return ideps;}

    class BCModification {
    public:
      int si;
      int tlonr;
      int bcnr;
      string * bcname;
    };

    NgArray<BCModification> bcmodifications;


    map<tuple<Surface*,Surface*>, string> named_edges;
      

    
    virtual int GenerateMesh (shared_ptr<Mesh> & mesh, MeshingParameters & mparam) override;

    void AddSplineSurface (shared_ptr<SplineSurface> ss) { spline_surfaces.Append(ss); }
  };


  


}

#endif

