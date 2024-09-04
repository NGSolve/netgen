#ifndef FILE_MVDRAW
#define FILE_MVDRAW


namespace netgen
{


  class VisualScene
  {
  protected:
    static NGGUI_API Point3d center;
    static NGGUI_API double rad;

    static double lookatmat[16];
    static double transmat[16];
    static double rotmat[16];
    static double centermat[16];
	
    static NGGUI_API double transformationmat[16];

    GLdouble clipplane[4];

    int changeval;
    static NGGUI_API GLdouble backcolor;

    static int NGGUI_API selface;
    static int selelement;
    static PointIndex NGGUI_API selpoint;
    static PointIndex selpoint2;
    static int locpi;
    static int NGGUI_API seledge;

    static optional<Point<3>> marker;

    static int subdivision_timestamp;
    static int subdivisions;
  public:
    static int viewport[4];
    static GLuint coltexname;
    static int ntexcols;
    int invcolor;


  public:
    NGGUI_API VisualScene ();
    NGGUI_API virtual ~VisualScene();

    NGGUI_API virtual void BuildScene (int zoomall = 0);
    NGGUI_API virtual void DrawScene ();
  
    NGGUI_API void CalcTransformationMatrices();
    NGGUI_API void StandardRotation (const char * dir);
    NGGUI_API void ArbitraryRotation (const NgArray<double> & alpha, const NgArray<Vec3d> & vec);
    NGGUI_API void ArbitraryRotation (const double alpha, const Vec3d & vec);

    NGGUI_API virtual void MouseMove(int oldx, int oldy,
                                      int newx, int newy,
                                      char mode);

    NGGUI_API void LookAt (const Point<3> & cam, const Point<3> & obj,
                 const Point<3> & camup);

    NGGUI_API void SetClippingPlane ();

    NGGUI_API virtual void MouseDblClick (int px, int py);

    NGGUI_API void SetLight ();
    static void SetBackGroundColor (double col)
    { backcolor = col; }

    NGGUI_API void CreateTexture (int ncols, int linear, double alpha, int typ);
    NGGUI_API void DrawColorBar (double minval, double maxval, int logscale = 0, bool linear = 1, string format="%8.3e", string unit="");
    NGGUI_API void DrawTitle (string title);
    NGGUI_API void DrawCoordinateCross ();
    NGGUI_API void DrawMarker();
    NGGUI_API void DrawNetgenLogo ();
    NGGUI_API void SetOpenGlColor(double val, double valmin, double valmax, int logscale = 0);
	

#ifdef PARALLELGL
    NGGUI_API void InitParallelGL ();
    NGGUI_API void Broadcast ();
#endif 
  };


  NGGUI_API extern void MyOpenGLText (const char * text);
  NGGUI_API extern int MyOpenGLTextWidth ();
  NGGUI_API extern void Set_OpenGLText_Callback ( void (*fun) (const char * text), int width );
  NGGUI_API extern VisualScene visual_scene_cross;
  NGGUI_API extern VisualScene *visual_scene;










  class VisualSceneSurfaceMeshing : public VisualScene
  {
    double scalex = 1., scaley = 1., shiftx = 0., shifty = 0.;
  public:
    shared_ptr<NgArray<Point<3>>> locpointsptr;
    shared_ptr<NgArray<INDEX_2>> loclinesptr;
    shared_ptr<NgArray<Point<2>>> plainpointsptr;
    int oldnl;
    bool clearptr;
    VisualSceneSurfaceMeshing ();
    virtual ~VisualSceneSurfaceMeshing ();

    void BuildScene (int zoomall = 0) override;
    void DrawScene () override;
    NGGUI_API void MouseMove(int oldx, int oldy, int newx, int newy,
                   char mode) override;
  };

  NGGUI_API extern VisualSceneSurfaceMeshing vssurfacemeshing;






  class VisualSceneMesh : public VisualScene
  {
    int filledlist = 0;
    int linelist = 0;
    int edgelist = 0;
    int pointnumberlist = 0;

    int tetlist = 0;
    int prismlist = 0;
    int pyramidlist = 0;
    int hexlist = 0;

    int badellist = 0;
    int identifiedlist = 0;
    int domainsurflist = 0;

    int vstimestamp = -1;
    int filledtimestamp = -1;
    int linetimestamp = -1;
    int edgetimestamp = -1;
    // int pointnumbertimestamp = -1;

    int tettimestamp = -1;
    int prismtimestamp = -1;
    int pyramidtimestamp = -1;
    int hextimestamp = -1;

    // int badeltimestamp = -1;
    // int identifiedtimestamp = -1;
    // int domainsurftimestamp = -1;

    struct {
      unsigned texture = -1;
      int width = 0;
      int height = 0;
      int size = 0;
    } colors;

    struct {
      unsigned framebuffer = 0;
      unsigned render_buffers[2];
      unsigned width = 0;
      unsigned height = 0;
      unsigned x = 0;
      unsigned y = 0;
      int list = 0;
      int list_timestamp = -1;
      double projmat[16];
      int viewport[4];
    } select;

#ifdef PARALLELGL
    NgArray<int> par_linelists;
    NgArray<int> par_filledlists;
#endif

    MouseEventHandler * user_me_handler;

    NgLock *lock;

    //  int selface, selelement;
    //  int selpoint, selpoint2, locpi;
    //  int seledge;

    double minh, maxh; // for meshsize coloring

    // weak_ptr<Mesh> wp_mesh;

  public:
    NGGUI_API VisualSceneMesh ();
	NGGUI_API virtual ~VisualSceneMesh ();

	NGGUI_API virtual void BuildScene (int zoomall = 0);
	NGGUI_API virtual void DrawScene ();
	NGGUI_API virtual void MouseDblClick (int px, int py);

    // void SetMesh (shared_ptr<Mesh> mesh) { wp_mesh = mesh; }
    // shared_ptr<Mesh> GetMesh () { return shared_ptr<Mesh>(wp_mesh); }
    shared_ptr<Mesh> GetMesh () const { return shared_ptr<Mesh>(global_mesh); }

    void SetMouseEventHandler (MouseEventHandler * handler)
    { user_me_handler = handler; }


	NGGUI_API int SelectedFace () const
    { return selface; }
	NGGUI_API void SetSelectedFace (int asf);
    //    { selface = asf; selecttimestamp = GetTimeStamp(); }

	NGGUI_API int SelectedEdge () const
    { return seledge; }
	NGGUI_API int SelectedElement () const
    { return selelement; }
	NGGUI_API int SelectedPoint () const
    { return selpoint; }
    void BuildFilledList (bool select);
    void BuildColorTexture();
    void SelectCenter(int zoomall);
    // private:
    void BuildLineList();
    void BuildEdgeList();
    void BuildPointNumberList();

    void BuildTetList(const BitArray & shownode);
    void BuildPrismList(const BitArray & shownode);
    void BuildPyramidList(const BitArray & shownode);
    void BuildHexList(const BitArray & shownode);

    void BuildBadelList();
    void BuildIdentifiedList();
    void BuildDomainSurfList();

    bool SelectSurfaceElement (int px, int py, Point<3> &p, bool select_on_clipping_plane);
    bool Unproject(int px, int py, Point<3> &p);
    ngcore::IVec<2> Project(Point<3> p);
  };

  NGGUI_API extern VisualSceneMesh vsmesh;
  

  class NGGUI_API VisualSceneSpecPoints : public VisualScene
  {
  public:
    VisualSceneSpecPoints ();
    virtual ~VisualSceneSpecPoints ();

    virtual void BuildScene (int zoomall = 0);
    virtual void DrawScene ();

    double len;
  };



  
  


  // extern struct Tcl_Interp * hinterp;


  extern void AddVisualizationScene (const string & name, 
                                     VisualScene * vs);


  void MouseDblClickSelect (const int px, const int py,
                            const GLdouble * clipplane, const GLdouble backcolor,
                            const double * transformationmat,
                            const Point3d & center,
                            const double rad,
                            const int displaylist,
                            int & selelement, int & selface, int & seledge, PointIndex & selpoint,
                            PointIndex & selpoint2, int & locpi);


  NGGUI_API std::vector<unsigned char> Snapshot( int w, int h );
}


#endif

