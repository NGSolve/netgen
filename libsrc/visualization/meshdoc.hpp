namespace netgen

{

class VisualSceneMeshDoctor : public VisualScene
{
  int filledlist;
  int outlinelist;
  int edgelist;

  int selelement, locpi;
  PointIndex selpoint, selpoint2;

  // for edgemarking:
  Array<int,PointIndex> edgedist;
  int markedgedist;
  

public:
  NGGUI_API VisualSceneMeshDoctor ();
  NGGUI_API virtual ~VisualSceneMeshDoctor ();

  NGGUI_API virtual void BuildScene (int zoomall = 0);
  NGGUI_API virtual void DrawScene ();
  NGGUI_API virtual void MouseDblClick (int px, int py);

  NGGUI_API void SetMarkEdgeDist (int dist);
  NGGUI_API void ClickElement (int elnr);
  NGGUI_API void UpdateTables ();
  NGGUI_API int IsSegmentMarked (int segnr) const;
};

class MeshDoctorParameters 
{
public:
  int active;
};


NGGUI_API extern MeshDoctorParameters meshdoctor;

}
