#ifndef LOCALH
#define LOCALH

/**************************************************************************/
/* File:   localh.hh                                                      */
/* Author: Joachim Schoeberl                                              */
/* Date:   29. Jan. 97                                                    */
/**************************************************************************/


namespace netgen
{


  /// box for grading
  class GradingBox
  {
    /// xmid
    float xmid[3];
    /// half edgelength
    float h2;
    ///
    GradingBox * childs[8] = {nullptr};
    ///
    GradingBox * father = nullptr;
    ///
    double hopt;
    ///
  public:

    struct 
    {
      /*
      unsigned int cutboundary:1;
      unsigned int isinner:1;
      unsigned int oldcell:1;
      unsigned int pinner:1;
      */
      bool cutboundary;
      bool isinner;
      bool oldcell;
      bool pinner;
    } flags;

    ///
    GradingBox (const double * ax1, const double * ax2);
    /// default constructor for Archive
    GradingBox() = default;
    ///
    void DeleteChilds();
    ///

    void DoArchive(Archive& ar);

    Point<3> PMid() const { return Point<3> (xmid[0], xmid[1], xmid[2]); }
    double H2() const { return h2; }
    double HOpt() const { return hopt; }

    bool HasChilds() const
    {
      for (int i = 0; i < 8; i++)
        if (childs[i]) return true;
      return false;
    }
    
    friend class LocalH;

    static BlockAllocator ball;
    void * operator new(size_t);
    void operator delete (void *);
  };




  /**
     Control of 3D mesh grading
  */
  class LocalH 
  {
    ///
    GradingBox * root;
    ///
    double grading;
    ///
    Array<GradingBox*> boxes;
    ///
    Box<3> boundingbox;
    /// octree or quadtree
    int dimension;
  public:
    ///
    DLL_HEADER LocalH (Point<3> pmin, Point<3> pmax, double grading, int adimension = 3);
    ///
    LocalH (const Box<3> & box, double grading, int adimension = 3)
      : LocalH (box.PMin(), box.PMax(), grading, adimension) { ; }
    /// Default ctor for archive
    LocalH() = default;

    DLL_HEADER ~LocalH();
    ///
    DLL_HEADER unique_ptr<LocalH> Copy();
    DLL_HEADER unique_ptr<LocalH> Copy( const Box<3> & bbox );
    ///
    DLL_HEADER void Delete();
    ///
    DLL_HEADER void DoArchive(Archive& ar);
    ///
    void SetGrading (double agrading) { grading = agrading; }
    ///
    DLL_HEADER void SetH (Point<3> x, double h);
    ///
    DLL_HEADER double GetH (Point<3> x) const;
    /// minimal h in box (pmin, pmax)
    DLL_HEADER double GetMinH (Point<3> pmin, Point<3> pmax) const;

    /// mark boxes intersecting with boundary-box
    // void CutBoundary (const Point3d & pmin, const Point3d & pmax)
    // { CutBoundaryRec (pmin, pmax, root); }

    void CutBoundary (const Box<3> & box)
    { CutBoundaryRec (box.PMin(), box.PMax(), root); }

    GradingBox * Find(Point<3> p) const;
  
    /// find inner boxes
    void FindInnerBoxes (const class AdFront3 & adfront,
			 int (*testinner)(const Point3d & p1));

    void FindInnerBoxes (const class AdFront2 & adfront,
			 int (*testinner)(const Point<2> & p1));


    /// clears all flags 
    void ClearFlags ()
    { ClearFlagsRec(root); }

    void ClearRootFlags ();

    /// widen refinement zone
    void WidenRefinement ();

    /// get points in inner elements
    void GetInnerPoints (NgArray<Point<3> > & points) const;
    void GetInnerPointsRec (const GradingBox * box, NgArray<Point<3> > & points) const;

    /// get points in outer closure
    void GetOuterPoints (NgArray<Point<3> > & points);

    ///
    void Convexify ();
    ///
    int GetNBoxes () { return boxes.Size(); } 
    const Box<3> & GetBoundingBox () const
    { return boundingbox; }
    ///
    void PrintMemInfo (ostream & ost) const;
  private:
    /// 
    double GetMinHRec (const Point3d & pmin, const Point3d & pmax,
		       const GradingBox * box) const;
    ///
    void CutBoundaryRec (const Point3d & pmin, const Point3d & pmax,
			 GradingBox * box);

    ///
    void FindInnerBoxesRec ( int (*inner)(const Point3d & p),
			     GradingBox * box);

    ///
    void FindInnerBoxesRec2 (GradingBox * box,
			     const class AdFront3 & adfront,
			     NgArray<Box3d> & faceboxes,
			     NgArray<int> & finds, int nfinbox);



    void FindInnerBoxesRec ( int (*inner)(const Point<2> & p),
			     GradingBox * box);

    ///
    void FindInnerBoxesRec2 (GradingBox * box,
			     const class AdFront2 & adfront,
			     FlatArray<Box<2>> faceboxes,
			     FlatArray<int> finds); // , int nfinbox);



    ///
    void SetInnerBoxesRec (GradingBox * box);

    ///
    void ClearFlagsRec (GradingBox * box);
  
    ///
    void ConvexifyRec (GradingBox * box);

    unique_ptr<LocalH> CopyRec( const Box<3> & bbox, GradingBox * current );

    friend ostream & operator<< (ostream & ost, const LocalH & loch);
  };




  inline ostream & operator<< (ostream & ost, const GradingBox & box)
  {
    ost << "gradbox, pmid = " << box.PMid() << ", h2 = " << box.H2() 
	<< " cutbound = " << box.flags.cutboundary << " isinner = " << box.flags.isinner 
	<< endl;
    return ost;
  }

  inline ostream & operator<< (ostream & ost, const LocalH & loch)
  {
    for (int i = 0; i < loch.boxes.Size(); i++)
      ost << "box[" << i << "] = " << *(loch.boxes[i]);
    return ost;
  }

}

#endif
