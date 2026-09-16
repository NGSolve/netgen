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
  inline int IsInArray(int n, const Array<int>& ia)
  {
    return ia.Contains(n); 
  }

  inline bool AddIfNotExists(Array<int>& list, int x)
  {
    if (list.Contains(x)) return false;
    list.Append(x);
    return true;
  }
  */
  
// extern DLL_HEADER MeshingParameters mparam;
  








  class STLEdgeDataList
  {
    Array<int> storedstatus;
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

    void BuildLineWithEdge(int ep1, int ep2, Array<IVec<2>>& line);
    void BuildClusterWithEdge(int ep1, int ep2, Array<IVec<2>>& line);

    int GetNEPPStat(int p, int status) const;
    int GetNConfCandEPP(int p) const;
  };






  class DLL_HEADER STLGeometry : public NetgenGeometry, public STLTopology
  {
    // edges to be meshed:
    Array<STLEdge> edges;
    //edges per point
    TABLE<int> edgesperpoint;

    // line: a connection of edges
    Array<STLLine*> lines;
    Array<int> lineendpoints; //per geometrypoint, 1 = is endpoint; 0 = no endpoint,

    Array<Vec<3>> normals; //normals belong to points!

    Array<IVec<2>> externaledges;

    int undoexternaledges;
    Array<IVec<2>> storedexternaledges;

    unique_ptr<STLEdgeDataList> edgedata;
    //  STLEdgeDataList edgedata_store;
    int calcedgedataanglesnew;

    int edgedatastored;



    int facecnt; 
    //meshpoint is only set, if an edge is at this point!!!

    Array<int> vicinity; //is one, if a triangle belongs to vicinity (eg. of selecttrig)
    Array<int> markedtrigs; //is one, if a triangle belongs to marked triangles (calcdirtystrigs)
    Array<Point<3>> markedsegs; //every pointpair is a segment!!!  
    Array<IVec<2>> selectedmultiedge;


    //spiralpoints:
    Array<int> spiralpoints;
    //
    Array<unique_ptr<STLChart>, ChartId> atlas;
    //marks all already charted trigs with chartnumber
    Array<ChartId, STLTrigId> chartmark; 
    //outerchartspertrig, ascending sorted
    TABLE<int> outerchartspertrig;


    //for meshing and project:
    Array<int> meshcharttrigs; //per trig: 1=belong to chart, 0 not
    mutable int meshchart;

    Array<int> ha_points;  // help array, np long, filled with 0 


    // sharp geometric edges not declared as edges
    // (not considered for spiral check)
    unique_ptr<ClosedHashTable<IVec<2>, int>> smoothedges;


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
    Array<STLLine*> meshlines;
    Array<Point<3>> meshpoints;

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
                          Point<3> & newp, EdgePointGeomInfo & newgi,
                          int edgenr) const override;



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
    //void AddSelectedMultiEdge(IVec<2> ep) {selectedmultiedge.Append(ep);}
    //int SelectedMultiEdgeSize() {return selectedmultiedge.Size();}
    const Array<IVec<2>>& SelectedMultiEdge() {return selectedmultiedge;}
    IVec<2> GetNearestSelectedDefinedEdge();
    void BuildSelectedMultiEdge(IVec<2> ep);
    void BuildSelectedEdge(IVec<2> ep);
    void BuildSelectedCluster(IVec<2> ep);

        void ImportEdges();
        void AddEdges(const Array<Point<3> >& eps);
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
    IVec<2> GetExternalEdge(int i) const {return externaledges[i-1];}

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
      ap1=markedsegs[i*2-2]; 
      ap2=markedsegs[i*2-1];
    }
    int GetNMarkedSegs() {return markedsegs.Size()/2;}
        void CalcVicinity(int starttrig);
        void GetVicinity(int starttrig, int size, Array<int>& vic);

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


    int AddNormal(const Vec<3>& n) { normals.Append(n); return normals.Size(); }
    const Vec<3> & GetNormal(int nr) const {return normals[nr-1];}
    void SetNormal(int nr, const Vec<3>& n) {normals[nr-1] = n;}

    int AddEdge(const STLEdge& v) { edges.Append(v); return edges.Size(); }
    int AddEdge(int p1, int p2);

    STLEdge GetEdge(int nr) {return edges[nr-1];}
    int GetNE() {return edges.Size();}

    double Area();

    double GetAngle(int t1, int t2);
    double GetGeomAngle(int t1, int t2);
    //if triangles t1 and t2 touch, return 1 and in p1, p2 the touching points
    //int TrigsTouch(int t1, int t2, int& p1, int& p2);


  
    ///

    ///ReadTriangle->STLTriangle, initialise some important variables, always after load!!!
    virtual void InitSTLGeometry (const Array<STLReadTriangle> & readtrigs) override;
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
                            Array<ChartId>& chartpointchecked, Array<int>& dirtytrigs);

    void ClearSpiralPoints();
    void SetSpiralPoint(int pn) {spiralpoints[pn-1] = 1;};
    int GetSpiralPoint(int pn) const {return spiralpoints[pn-1];};

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
  
    void GetInnerChartLimes(Array<IVec<2>>& limes, ChartId chartnum);

    //FOR MESHING
    int GetMeshChartNr () { return meshchart; }
    void GetMeshChartBoundary (Array<Point<2>> & points,
                               Array<Point<3>> & points3d,
                               Array<IVec<2>> & lines, double h);


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
    STLLine* GetLine(int nr) const {return lines[nr-1];}
    int GetLineP(int lnr, int pnr) const {return lines[lnr-1]->PNum(pnr);}
    int GetLineNP(int nr) const {return lines[nr-1]->NP();}

    void SetLineEndPoint(int pn);
    int IsLineEndPoint(int pn);
    int LineEndPointsSet() const {return lineendpoints.Size() == GetNP();}
    void ClearLineEndPoints();

    void RestrictLocalH(class Mesh & mesh, double gh, const STLParameters& stlparam, const MeshingParameters& mparam);
    void RestrictLocalHCurv(class Mesh & mesh, double gh, const STLParameters& stlparam);
    void RestrictHChartDistOneChart(ChartId chartnum, Array<int>& acttrigs, class Mesh & mesh, 
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
