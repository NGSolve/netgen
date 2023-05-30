#ifndef FILE_STLGEOM
#define FILE_STLGEOM

/**************************************************************************/
/* File:   stlgeom.hpp                                                    */
/* Author: Joachim Schoeberl                                              */
/* Author2: Johannes Gerstmayr                                            */
/* Date:   26. Jul. 99                                                    */
/**************************************************************************/

/**
   STL Geometry


   Terminology:
   
   Point ... coordinates of STL triangles
   Triangle  (short Trig)  STL triangle
   TopEdge .... edge in topology, boundary of STL triangles (many)
   Edge .... Edges which will occur in the mesh (confirmed edges, less)
*/


#include <meshing.hpp>

#include "stltopology.hpp"
#include "stltool.hpp"
#include "stlline.hpp"
 


namespace netgen
{
  /*
  inline int IsInArray(int n, const NgArray<int>& ia)
  {
    return ia.Contains(n); 
  }

  inline bool AddIfNotExists(NgArray<int>& list, int x)
  {
    if (list.Contains(x)) return false;
    list.Append(x);
    return true;
  }
  */
  
// extern DLL_HEADER MeshingParameters mparam;
  








  class STLEdgeDataList
  {
    NgArray<int> storedstatus;
    STLTopology & geom;
  public:
  
    STLEdgeDataList(STLTopology & ageom);
    ~STLEdgeDataList();

    void Store ();
    void Restore ();

    void SetSize(int /* size */) { };
    void Clear() { };
    int Size() const { return geom.GetNTE(); }
    const STLTopEdge & Get(int i) const { return geom.GetTopEdge(i); }
    STLTopEdge & Elem(int i) { return geom.GetTopEdge(i); }

    int GetNEPP(int pn) const {return geom.NTopEdgesPerPoint(pn); }
    int GetEdgePP(int pn, int vi) const {return geom.TopEdgePerPoint(pn, vi);};

    //void AddEdgePP(int pn, int vn) { } ;

    void ResetAll();
    void ChangeStatus(int status1, int status2);

    int GetEdgeNum(int np1, int np2) const
    { return geom.GetTopEdgeNum (np1, np2); }

    int GetNConfEdges() const;

    void Write(ofstream& of) const;
    void Read(ifstream& ifs);

    void BuildLineWithEdge(int ep1, int ep2, NgArray<twoint>& line);
    void BuildClusterWithEdge(int ep1, int ep2, NgArray<twoint>& line);

    int GetNEPPStat(int p, int status) const;
    int GetNConfCandEPP(int p) const;
  };






  class DLL_HEADER STLGeometry : public NetgenGeometry, public STLTopology
  {
    // edges to be meshed:
    NgArray<STLEdge> edges;
    //edges per point
    TABLE<int> edgesperpoint;

    // line: a connection of edges
    NgArray<STLLine*> lines;
    NgArray<int> lineendpoints; //per geometrypoint, 1 = is endpoint; 0 = no endpoint,

    NgArray<Vec3d> normals; //normals belong to points!

    NgArray<twoint> externaledges;

    int undoexternaledges;
    NgArray<twoint> storedexternaledges;

    unique_ptr<STLEdgeDataList> edgedata;
    //  STLEdgeDataList edgedata_store;
    int calcedgedataanglesnew;

    int edgedatastored;



    int facecnt; 
    //meshpoint is only set, if an edge is at this point!!!

    NgArray<int> vicinity; //is one, if a triangle belongs to vicinity (eg. of selecttrig)
    NgArray<int> markedtrigs; //is one, if a triangle belongs to marked triangles (calcdirtystrigs)
    NgArray<Point3d> markedsegs; //every pointpair is a segment!!!  
    NgArray<twoint> selectedmultiedge;


    //spiralpoints:
    NgArray<int> spiralpoints;
    //
    Array<unique_ptr<STLChart>, ChartId> atlas;
    //marks all already charted trigs with chartnumber
    Array<ChartId, STLTrigId> chartmark; 
    //outerchartspertrig, ascending sorted
    TABLE<int> outerchartspertrig;


    //for meshing and project:
    NgArray<int> meshcharttrigs; //per trig: 1=belong to chart, 0 not
    mutable int meshchart;

    NgArray<int> ha_points;  // help array, np long, filled with 0 


    // sharp geometric edges not declared as edges
    // (not considered for spiral check)
    INDEX_2_HASHTABLE<int> * smoothedges;


    //transformation:
    mutable Vec<3> meshtrignv;
    Vec<3> ex, ey, ez;
    Point<3> p1;

  public:
    int edgesfound;
    int surfacemeshed;
    int surfaceoptimized;
    int volumemeshed;

    int trigsconverted; //when STLTriangles exist -> 1

    //for selecting nodes
    //int selecttrig, nodeofseltrig;

    //only for testing;
    NgArray<STLLine*> meshlines;
    NgArray<Point3d> meshpoints;

    double area;
  public:
    STLGeometry();
    virtual ~STLGeometry();

    void DoArchive(Archive& ar) override
    {
      STLTopology::DoArchive(ar);
    }

    void Clear();

    virtual void Save (const filesystem::path & filename) const override;

    bool CalcPointGeomInfo(int surfind, PointGeomInfo& gi, const Point<3> & p3) const override;
    PointGeomInfo ProjectPoint(INDEX surfind, Point<3> & p) const override;
    bool ProjectPointGI (int surfind, Point<3> & p, PointGeomInfo & gi) const override;
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



	void STLInfo(double* data);
    //stldoctor:
	void SmoothNormals(const STLParameters& stlparam);
	void MarkNonSmoothNormals(const STLParameters& stlparam);

	void CalcEdgeData();
	void CalcEdgeDataAngles();

    const STLEdgeDataList& EdgeDataList() const {return *edgedata;}

	void UndoEdgeChange();
	void StoreEdgeData();
	void RestoreEdgeData();

    //void ClearSelectedMultiEdge() {selectedmultiedge.SetSize(0);}
    //void AddSelectedMultiEdge(twoint ep) {selectedmultiedge.Append(ep);}
    //int SelectedMultiEdgeSize() {return selectedmultiedge.Size();}
    const NgArray<twoint>& SelectedMultiEdge() {return selectedmultiedge;}
    twoint GetNearestSelectedDefinedEdge();
    void BuildSelectedMultiEdge(twoint ep);
    void BuildSelectedEdge(twoint ep);
    void BuildSelectedCluster(twoint ep);

	void ImportEdges();
	void AddEdges(const NgArray<Point<3> >& eps);
	void ExportEdges();
	void LoadEdgeData(const filesystem::path & file);
	void SaveEdgeData(const filesystem::path & file);
    //  void SetEdgeAtSelected(int mode);
  

	void STLDoctorConfirmEdge();
	void STLDoctorCandidateEdge();
	void STLDoctorExcludeEdge();
	void STLDoctorUndefinedEdge();

	void STLDoctorSetAllUndefinedEdges();
	void STLDoctorEraseCandidateEdges();
	void STLDoctorConfirmCandidateEdges();
	void STLDoctorConfirmedToCandidateEdges();

	void STLDoctorDirtyEdgesToCandidates();
	void STLDoctorLongLinesToCandidates();

	void UndoExternalEdges();
	void StoreExternalEdges();
	void RestoreExternalEdges();

	void ImportExternalEdges(const char * filename);  // Flame edges, JS
    //  void LoadExternalEdges();

	void BuildExternalEdgesFromEdges();
	void SaveExternalEdges();
	void AddExternalEdgeAtSelected();
	void AddClosedLinesToExternalEdges();
	void AddLongLinesToExternalEdges();
	void AddAllNotSingleLinesToExternalEdges();
	void STLDoctorBuildEdges(const STLParameters& stlparam);
	void AddExternalEdgesFromGeomLine();
	void DeleteDirtyExternalEdges();
	void DeleteExternalEdgeAtSelected();
	void DeleteExternalEdgeInVicinity();
    void AddExternalEdge(int p1, int p2);
    void DeleteExternalEdge(int p1, int p2);
    int IsExternalEdge(int p1, int p2);
    int NOExternalEdges() const {return externaledges.Size();}
    twoint GetExternalEdge(int i) const {return externaledges.Get(i);}

	void DestroyDirtyTrigs();
	void CalcNormalsFromGeometry();
	void MoveSelectedPointToMiddle();
	void NeighbourAnglesOfSelectedTrig();
	void PrintSelectInfo();
	void ShowSelectedTrigChartnum();
	void ShowSelectedTrigCoords();
	void SmoothGeometry ();


	void LoadMarkedTrigs();
	void SaveMarkedTrigs();
	void ClearMarkedSegs() {markedsegs.SetSize(0);}
    void AddMarkedSeg(const Point<3> & ap1, const Point<3> & ap2) 
    {
      markedsegs.Append(ap1);markedsegs.Append(ap2);
    }

    void GetMarkedSeg(int i, Point<3> & ap1, Point<3> & ap2) 
    {
      ap1=markedsegs.Get(i*2-1); 
      ap2=markedsegs.Get(i*2);
    }
    int GetNMarkedSegs() {return markedsegs.Size()/2;}
	void CalcVicinity(int starttrig);
	void GetVicinity(int starttrig, int size, NgArray<int>& vic);

	int Vicinity(int trig) const;

	void InitMarkedTrigs();
	void MarkDirtyTrigs(const STLParameters& stlparam);
	void SmoothDirtyTrigs(const STLParameters& stlparam);
	void GeomSmoothRevertedTrigs(const STLParameters& stlparam);
	void MarkRevertedTrigs(const STLParameters& stlparam);
	double CalcTrigBadness(int i);
	int IsMarkedTrig(int trig) const;
	void SetMarkedTrig(int trig, int num);
	void MarkTopErrorTrigs ();

    //Selected triangle
	void SetSelectTrig(int trig);
	int GetSelectTrig() const;
	void SetNodeOfSelTrig(int n);
	int GetNodeOfSelTrig() const;


    int AddNormal(const Vec3d& n) { normals.Append(n); return normals.Size(); }
    const Vec3d & GetNormal(int nr) const {return normals.Get(nr);}
    void SetNormal(int nr, const Vec3d& n) {normals.Elem(nr) = n;}

    int AddEdge(const STLEdge& v) { edges.Append(v); return edges.Size(); }
    int AddEdge(int p1, int p2);

    STLEdge GetEdge(int nr) {return edges.Get(nr);}
    int GetNE() {return edges.Size();}

    double Area();

    double GetAngle(int t1, int t2);
    double GetGeomAngle(int t1, int t2);
    //if triangles t1 and t2 touch, return 1 and in p1, p2 the touching points
    //int TrigsTouch(int t1, int t2, int& p1, int& p2);


  
    ///

    ///ReadTriangle->STLTriangle, initialise some important variables, always after load!!!
    virtual void InitSTLGeometry (const NgArray<STLReadTriangle> & readtrigs) override;
    virtual void TopologyChanged() override; //do some things, if topology changed!
    int CheckGeometryOverlapping();

    //get NO edges per point
    int GetEPPSize() const {return edgesperpoint.Size();};
    int GetNEPP(int pn) 
    {
      if (edgesperpoint.Size() == 0) {BuildEdgesPerPoint();}
      return edgesperpoint.EntrySize(pn);
    };
    int GetEdgePP(int pn, int vi)
    {
      if (edgesperpoint.Size() == 0) {BuildEdgesPerPoint();}
      return edgesperpoint.Get(pn,vi);
    };
    void AddEdgePP(int pn, int vn) {edgesperpoint.Add1(pn,vn);};
    //von 2 punkten ermitteln, ob sie eine Kante sind
    int IsEdge(int p1, int p2);
    int IsEdgeNum(int p1, int p2);

    ///Build EdgeSegments
    void ClearEdges();
    void BuildEdges(const STLParameters& stlparam);
    void BuildEdgesPerPoint();
    void UseExternalEdges();


    void FindEdgesFromAngles(const STLParameters& stlparam);
    void CalcFaceNums();
    int GetNOBodys();
    int GetNOFaces() {return facecnt;}
    void LinkEdges(const STLParameters& stlparam);

    void AddConeAndSpiralEdges(const STLParameters& stlparam);
    void AddFaceEdges(); //each face should have at least one starting edge (outherwise it won't be meshed)

    void GetDirtyChartTrigs(int chartnum, STLChart& chart, const Array<ChartId, STLTrigId>& outercharttrigs, 
			    NgArray<ChartId>& chartpointchecked, NgArray<int>& dirtytrigs);

    void ClearSpiralPoints();
    void SetSpiralPoint(int pn) {spiralpoints.Elem(pn) = 1;};
    int GetSpiralPoint(int pn) const {return spiralpoints.Get(pn);};

    void GetSortedTrianglesAroundPoint(STLPointId p, STLTrigId starttrig, Array<STLTrigId>& trigs);

    // smooth edges: sharp geometric edges not declared as edges
    void BuildSmoothEdges ();
    bool IsSmoothEdge (int pi1, int pi2) const;


    //make charts with regions of a max. angle
    void MakeAtlas(class Mesh & mesh, const MeshingParameters& mparam, const STLParameters& stlparam);

    //outerchartspertrig, sorted!
    int GetOCPTSize() const {return outerchartspertrig.Size();};
    int GetNOCPT(int tn) const {return outerchartspertrig.EntrySize(tn);};
    int GetOCPT(int tn, int vi) const {return outerchartspertrig.Get(tn,vi);};
    void SetOCPT(int tn, int vi, int ocn) {outerchartspertrig.Set(tn,vi,ocn);};
    void AddOCPT(int tn, int ocn) {outerchartspertrig.Add1(tn, ocn);};
    int TrigIsInOC(int tn, int ocn) const;
 
    //get chart number of a trig or 0 if unmarked
    ChartId GetChartNr(STLTrigId i) const;
    ChartId GetMarker(STLTrigId i) const  { return chartmark[i]; }
    void SetMarker(STLTrigId nr, ChartId m);
    size_t GetNOCharts() const { return atlas.Size(); }
    //get a chart from atlas
    const STLChart& GetChart(ChartId nr) const { return *atlas[nr];};
    STLChart & GetChart(ChartId nr) { return *atlas[nr];};
    int AtlasMade() const;
  
    void GetInnerChartLimes(NgArray<twoint>& limes, ChartId chartnum);

    //FOR MESHING
    int GetMeshChartNr () { return meshchart; }
    void GetMeshChartBoundary (NgArray<Point<2>> & points,
			       NgArray<Point<3>> & points3d,
			       NgArray<INDEX_2> & lines, double h);


    Point<3> PointBetween(const Point<3> & p1, int t1, const Point<3> & p2, int t2);

    //select triangles in meshcharttrigs of actual (defined by trig) whole chart
    void PrepareSurfaceMeshing();
    //
    void DefineTangentialPlane(const Point<3> & ap1, const Point<3> & ap2, int trig);
    //
    void SelectChartOfTriangle (int trignum) const;
    //
    void SelectChartOfPoint (const Point<3> & p);
    //
    const Vec<3> & GetChartNormalVector () const { return meshtrignv; }

    // list of trigs
    void ToPlane (const Point<3> & locpoint, int * trigs, Point<2> & plainpoint, 
		  double h, int& zone, int checkchart);
    //return 0, wenn alles OK, 1 sonst
    int FromPlane (const Point<2> & plainpoint, Point<3> & locpoint, double h);
  
    //get nearest point in actual chart and return any triangle where it lies on
    int ProjectNearest(Point<3> & p3d) const;
    //project point with normal nv from last define tangential plane

    int LastTrig() const;
    int Project(Point<3> & p3d) const;
    int ProjectOnWholeSurface (Point<3> & p3d) const;

    int GetNLines() const {return lines.Size();}
    int AddLine(STLLine* line) { lines.Append(line); return lines.Size(); }
    STLLine* GetLine(int nr) const {return lines.Get(nr);}
    int GetLineP(int lnr, int pnr) const {return lines.Get(lnr)->PNum(pnr);}
    int GetLineNP(int nr) const {return lines.Get(nr)->NP();}

    void SetLineEndPoint(int pn);
    int IsLineEndPoint(int pn);
    int LineEndPointsSet() const {return lineendpoints.Size() == GetNP();}
    void ClearLineEndPoints();

    void RestrictLocalH(class Mesh & mesh, double gh, const STLParameters& stlparam, const MeshingParameters& mparam);
    void RestrictLocalHCurv(class Mesh & mesh, double gh, const STLParameters& stlparam);
    void RestrictHChartDistOneChart(ChartId chartnum, NgArray<int>& acttrigs, class Mesh & mesh, 
				    double gh, double fact, double minh, const STLParameters& stlparam);

    friend class MeshingSTLSurface;

    int GenerateMesh (shared_ptr<Mesh> & mesh, MeshingParameters & mparam) override;
    
    // Add additional Point to chart to close the surface and write the resulting stl to a file
    void WriteChartToFile( ChartId chartnumber, filesystem::path filename="chart.slb" );
  };
 

#include "meshstlsurface.hpp"



extern int STLMeshingDummy (STLGeometry* stlgeometry, shared_ptr<Mesh> & mesh, const MeshingParameters & mparam,
                            const STLParameters& stlpar);


}
#endif
