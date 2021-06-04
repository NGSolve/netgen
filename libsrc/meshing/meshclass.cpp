#include <mystdlib.h>
#include <atomic>
#include <set>
#include "meshing.hpp"

#ifdef NG_PYTHON
// must be included to instantiate Archive::Shallow(NetgenGeometry&)
#include <core/python_ngcore.hpp>
#endif

namespace netgen
{

  static mutex buildsearchtree_mutex;

  Mesh :: Mesh ()
    : topology(*this), surfarea(*this)
  {
    boundaryedges = nullptr;
    surfelementht = nullptr; 
    segmentht = nullptr;

    lochfunc = nullptr;
    // mglevels = 1;
    elementsearchtree = nullptr;
    elementsearchtreets = NextTimeStamp();
    majortimestamp = timestamp = NextTimeStamp();
    hglob = 1e10;
    hmin = 0;
    numvertices = -1;
    dimension = 3;

    curvedelems = make_unique<CurvedElements> (*this);
    clusters = make_unique<AnisotropicClusters> (*this);
    ident = make_unique<Identifications> (*this);

    hpelements = NULL;
    coarsemesh = NULL;

    ps_startelement = 0;

    geomtype = NO_GEOM;

    bcnames.SetSize(0);
    cd2names.SetSize(0);

    // this->comm = netgen :: ng_comm;
#ifdef PARALLEL
    paralleltop = make_unique<ParallelMeshTopology> (*this);
#endif
  }


  Mesh :: ~Mesh()
  {
    // delete lochfunc;
    // delete boundaryedges;
    // delete surfelementht;
    // delete segmentht;
    // delete curvedelems;
    // delete clusters;
    // delete ident;
    // delete elementsearchtree;
    // delete coarsemesh;
    // delete hpelements;

    for (int i = 0; i < materials.Size(); i++)
      delete materials[i];
    for(int i = 0; i < userdata_int.Size(); i++)
      delete userdata_int[i];
    for(int i = 0; i < userdata_double.Size(); i++)
      delete userdata_double[i];

    for (int i = 0; i < bcnames.Size(); i++ )
      delete bcnames[i];

    for (int i = 0; i < cd2names.Size(); i++)
      delete cd2names[i];

    // #ifdef PARALLEL
    // delete paralleltop;
    // #endif
  }

  void Mesh :: SetCommunicator(NgMPI_Comm acomm)
  {
    this->comm = acomm;
  }

  Mesh & Mesh :: operator= (const Mesh & mesh2)
  {
    geometry = mesh2.geometry;
    dimension = mesh2.dimension;
    points = mesh2.points;
    segments = mesh2.segments;
    surfelements = mesh2.surfelements;
    volelements = mesh2.volelements;
    lockedpoints = mesh2.lockedpoints;
    facedecoding = mesh2.facedecoding;
    dimension = mesh2.dimension;


    materials.SetSize( mesh2.materials.Size() );
    for ( int i = 0; i < mesh2.materials.Size(); i++ )
      if ( mesh2.materials[i] ) materials[i] = new string ( *mesh2.materials[i] );
      else materials[i] = 0;

    std::map<const string*, string*> bcmap;
    bcnames.SetSize( mesh2.bcnames.Size() );
    for ( int i = 0; i < mesh2.bcnames.Size(); i++ )
    {
      if ( mesh2.bcnames[i] ) bcnames[i] = new string ( *mesh2.bcnames[i] );
      else bcnames[i] = 0;
      bcmap[mesh2.bcnames[i]] = bcnames[i];
    }

    // Remap string* members in FaceDescriptor to new mesh
    for (auto & f : facedecoding)
      f.SetBCName( bcmap[&f.GetBCName()] );


    cd2names.SetSize(mesh2.cd2names.Size());
    for (int i=0; i < mesh2.cd2names.Size(); i++)
      if (mesh2.cd2names[i]) cd2names[i] = new string(*mesh2.cd2names[i]);
      else cd2names[i] = 0;

    cd3names.SetSize(mesh2.cd3names.Size());
    for (int i=0; i < mesh2.cd3names.Size(); i++)
      if (mesh2.cd3names[i]) cd3names[i] = new string(*mesh2.cd3names[i]);
      else cd3names[i] = 0;

    numvertices = mesh2.numvertices;

    return *this;
  }


  void Mesh :: DeleteMesh()
  {
    NgLock lock(mutex);
    lock.Lock();
    points.SetSize(0);
    segments.SetSize(0);
    surfelements.SetSize(0);
    volelements.SetSize(0);
    lockedpoints.SetSize(0);
    // surfacesonnode.SetSize(0);

    // delete boundaryedges;
    boundaryedges = nullptr;
    segmentht = nullptr;
    surfelementht = nullptr;

    openelements.SetSize(0);
    facedecoding.SetSize(0);

    ident = make_unique<Identifications> (*this);
    topology = MeshTopology (*this);
    curvedelems = make_unique<CurvedElements> (*this);
    clusters = make_unique<AnisotropicClusters> (*this);

    for ( int i = 0; i < bcnames.Size(); i++ )
      if ( bcnames[i] ) delete bcnames[i];
    for (int i= 0; i< cd2names.Size(); i++)
      if (cd2names[i]) delete cd2names[i];

#ifdef PARALLEL
    paralleltop = make_unique<ParallelMeshTopology> (*this);
#endif

    lock.UnLock();

    timestamp = NextTimeStamp();
  }


  void Mesh :: ClearSurfaceElements()
  { 
    surfelements.SetSize(0);
    /*
    for (int i = 0; i < facedecoding.Size(); i++)
      facedecoding[i].firstelement = -1;
    */
    for (auto & fd : facedecoding)
      fd.firstelement = -1;
    
    timestamp = NextTimeStamp();
  }



  PointIndex Mesh :: AddPoint (const Point3d & p, int layer)
  { 
    return AddPoint (p, layer, INNERPOINT);
  }

  PointIndex Mesh :: AddPoint (const Point3d & p, int layer, POINTTYPE type)
  { 

    // PointIndex pi = points.End();
    PointIndex pi = *points.Range().end();
    if (points.Size() == points.AllocSize())
      {
        NgLock lock(mutex);
        lock.Lock();
        points.Append ( MeshPoint (p, layer, type) ); 
        lock.UnLock();
      }
    else
      {
        points.Append ( MeshPoint (p, layer, type) ); 
      }

    timestamp = NextTimeStamp();

    return pi;
  }


  SegmentIndex Mesh :: AddSegment (const Segment & s)
  { 
    NgLock lock(mutex);	
    lock.Lock();
    timestamp = NextTimeStamp();

    int maxn = max2 (s[0], s[1]);
    maxn += 1-PointIndex::BASE;

    /*
      if (maxn > ptyps.Size())
      {
      int maxo = ptyps.Size();
      ptyps.SetSize (maxn);
      for (int i = maxo; i < maxn; i++)
      ptyps[i] = INNERPOINT;
      }

      if (ptyps[s[0]] > EDGEPOINT) ptyps[s[0]] = EDGEPOINT;
      if (ptyps[s[1]] > EDGEPOINT) ptyps[s[1]] = EDGEPOINT;
    */

    if (maxn <= points.Size())
      {
        if (points[s[0]].Type() > EDGEPOINT)
          points[s[0]].SetType (EDGEPOINT);
        if (points[s[1]].Type() > EDGEPOINT)
          points[s[1]].SetType (EDGEPOINT);
      }
    /*
      else
      {
      cerr << "edge points nrs > points.Size" << endl;
      }
    */

    SegmentIndex si = segments.Size();
    segments.Append (s); 

    lock.UnLock();
    return si;
  }

  SurfaceElementIndex Mesh :: AddSurfaceElement (const Element2d & el)
  {     
    timestamp = NextTimeStamp();

    PointIndex maxn = el[0];
    for (int i = 1; i < el.GetNP(); i++)
      if (el[i] > maxn) maxn = el[i];

    /*
    maxn += 1-PointIndex::BASE;
    if (maxn <= points.Size())
      {
        for (int i = 0; i < el.GetNP(); i++)
          if (points[el[i]].Type() > SURFACEPOINT)
            points[el[i]].SetType(SURFACEPOINT);
      }
    */
    // if (maxn < points.End())
    if (maxn < *points.Range().end())
      for (PointIndex pi : el.PNums())
        if (points[pi].Type() > SURFACEPOINT)
          points[pi].SetType(SURFACEPOINT);

    
    SurfaceElementIndex si = surfelements.Size();
    if (surfelements.AllocSize() == surfelements.Size())
      {
        NgLock lock(mutex);
        lock.Lock();
        surfelements.Append (el);
        lock.UnLock();
      }
    else
      {
        surfelements.Append (el);        
      }

    if (el.index<=0 || el.index > facedecoding.Size())
      cerr << "has no facedecoding: fd.size = " << facedecoding.Size() << ", ind = " << el.index << endl;

    surfelements.Last().next = facedecoding[el.index-1].firstelement;
    facedecoding[el.index-1].firstelement = si;

    if (SurfaceArea().Valid())
      SurfaceArea().Add (el);

    return si;
  }

  void Mesh :: SetSurfaceElement (SurfaceElementIndex sei, const Element2d & el)
  {
    int maxn = el[0];
    for (int i = 1; i < el.GetNP(); i++)
      if (el[i] > maxn) maxn = el[i];

    maxn += 1-PointIndex::BASE;

    if (maxn <= points.Size())
      {
        for (int i = 0; i < el.GetNP(); i++)
          if (points[el[i]].Type() > SURFACEPOINT)
            points[el[i]].SetType(SURFACEPOINT);
      }

    surfelements[sei] = el;
    if (el.index > facedecoding.Size())
      cerr << "has no facedecoding: fd.size = " << facedecoding.Size() << ", ind = " << el.index << endl;

    // add lock-free to list ... slow, call RebuildSurfaceElementLists later
    /*
    surfelements[sei].next = facedecoding[el.index-1].firstelement;
    auto & head = reinterpret_cast<atomic<SurfaceElementIndex>&> (facedecoding[el.index-1].firstelement);
    while (!head.compare_exchange_weak (surfelements[sei].next, sei))
      ;
    */

    /*
    if (SurfaceArea().Valid())
      SurfaceArea().Add (el);
    */
  }


  ElementIndex Mesh :: AddVolumeElement (const Element & el)
  { 
    /*
    int maxn = el[0];
    for (int i = 1; i < el.GetNP(); i++)
      if (el[i] > maxn) maxn = el[i];

    maxn += 1-PointIndex::BASE;
    */
    
    /*
      if (maxn > ptyps.Size())
      {
      int maxo = ptyps.Size();
      ptyps.SetSize (maxn);
      for (i = maxo+PointIndex::BASE; 
      i < maxn+PointIndex::BASE; i++)
      ptyps[i] = INNERPOINT;
      }
    */
    /*
      if (maxn > points.Size())
      {
      cerr << "add vol element before point" << endl;
      }
    */

    int ve = volelements.Size();

    if (volelements.Size() == volelements.AllocSize())
      {
        NgLock lock(mutex);
        lock.Lock();
        volelements.Append (el);
        lock.UnLock();
      }
    else
      {
        volelements.Append (el);
      }
    volelements.Last().flags.illegal_valid = 0;

    // while (volelements.Size() > eltyps.Size())
    // eltyps.Append (FREEELEMENT);

    timestamp = NextTimeStamp();

    return ve;
  }

  void Mesh :: SetVolumeElement (ElementIndex ei, const Element & el)
  {
    /*
    int maxn = el[0];
    for (int i = 1; i < el.GetNP(); i++)
      if (el[i] > maxn) maxn = el[i];

    maxn += 1-PointIndex::BASE;
    */

    volelements[ei]  = el;
    volelements[ei].flags.illegal_valid = 0;
  }





  void Mesh :: Save (const string & filename) const
  {
    if (filename.find(".vol.bin") != string::npos)
    {
        BinaryOutArchive in(filename);
        in & const_cast<Mesh&>(*this);
        return;
    }

    ostream * outfile;
    if (filename.find(".vol.gz")!=string::npos)
      outfile = new ogzstream(filename.c_str());
    else if (filename.find(".vol")!=string::npos)
      outfile = new ofstream(filename.c_str());
    else
      outfile = new ogzstream((filename+".vol.gz").c_str());

    Save(*outfile);
    delete outfile;
  }



  void Mesh :: Save (ostream & outfile) const
  {
    int i, j;

    double scale = 1;  // globflags.GetNumFlag ("scale", 1);
    int inverttets = 0;  // globflags.GetDefineFlag ("inverttets");
    int invertsurf = 0;  // globflags.GetDefineFlag ("invertsurfacemesh");

    outfile << "# Generated by NETGEN " << GetLibraryVersion("netgen") << endl << endl;


    outfile << "mesh3d" << "\n";

    outfile << "dimension\n" << GetDimension() << "\n";

    outfile << "geomtype\n" << int(geomtype) << "\n";


    outfile << "\n";
    outfile << "# surfnr    bcnr   domin  domout      np      p1      p2      p3"
            << "\n";


    switch (geomtype)
      {
      case GEOM_STL:
        outfile << "surfaceelementsgi" << "\n";
        break;
      case GEOM_OCC: case GEOM_ACIS:
        outfile << "surfaceelementsuv" << "\n";
        break;
      default:
        outfile << "surfaceelements" << "\n";
      }

    outfile << GetNSE() << "\n";

    for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
      {
        if ((*this)[sei].GetIndex())
          {
            outfile << " " << GetFaceDescriptor((*this)[sei].GetIndex ()).SurfNr()+1;
            outfile << " " << GetFaceDescriptor((*this)[sei].GetIndex ()).BCProperty();
            outfile << " " << GetFaceDescriptor((*this)[sei].GetIndex ()).DomainIn();
            outfile << " " << GetFaceDescriptor((*this)[sei].GetIndex ()).DomainOut();
          }
        else
          outfile << " 0 0 0";

        Element2d sel = (*this)[sei];
        if (invertsurf)
          sel.Invert();

        outfile << " " << sel.GetNP();
        for (j = 0; j < sel.GetNP(); j++)
          outfile << " " << sel[j];

        switch (geomtype)
          {
          case GEOM_STL:
            for (j = 1; j <= sel.GetNP(); j++)
              outfile << " " << sel.GeomInfoPi(j).trignum;
            break;
          case GEOM_OCC: case GEOM_ACIS:
            for (j = 1; j <= sel.GetNP(); j++)
              {
                outfile << " " << sel.GeomInfoPi(j).u;
                outfile << " " << sel.GeomInfoPi(j).v;
              }
            break;
          default:
            ; 
          }
        outfile << "\n";
      }

    outfile << "\n" << "\n";
    outfile << "#  matnr      np      p1      p2      p3      p4" << "\n";
    outfile << "volumeelements" << "\n";
    outfile << GetNE() << "\n";

    for (ElementIndex ei = 0; ei < GetNE(); ei++)
      {
        outfile << (*this)[ei].GetIndex();
        outfile << " " << (*this)[ei].GetNP();

        Element el = (*this)[ei];
        if (inverttets) el.Invert();

        for (j = 0; j < el.GetNP(); j++)
	  outfile << " " << el[j];
        outfile << "\n";
      }


    outfile << "\n" << "\n";
    //     outfile << "   surf1   surf2      p1      p2" << "\n";
    outfile << "# surfid  0   p1   p2   trignum1    trignum2   domin/surfnr1    domout/surfnr2   ednr1   dist1   ednr2   dist2 \n";
    outfile << "edgesegmentsgi2" << "\n";
    outfile << GetNSeg() << "\n";

    for (i = 1; i <= GetNSeg(); i++)
      {
        const Segment & seg = LineSegment (i);
        outfile.width(8);
        outfile << seg.si; // 2D: bc number, 3D: wievielte Kante
        outfile.width(8);
        outfile << 0;
        outfile.width(8);
        outfile << seg[0];
        outfile.width(8);
        outfile << seg[1];
        outfile << " ";
        outfile.width(8);
        outfile << seg.geominfo[0].trignum;  // stl dreiecke
        outfile << " ";
        outfile.width(8);
        outfile << seg.geominfo[1].trignum; // << endl;  // stl dreieck

        if (dimension == 3)
          {
            outfile << " ";
            outfile.width(8);
            outfile << seg.surfnr1+1;
            outfile << " ";
            outfile.width(8);
            outfile << seg.surfnr2+1;
          }
        else
          {
            outfile << " ";
            outfile.width(8);
            outfile << seg.domin;
            outfile << " ";
            outfile.width(8);
            outfile << seg.domout;
          }

        outfile << " ";
        outfile.width(8);
        outfile << seg.edgenr;
        outfile << " ";
        outfile.width(12);
        outfile.precision(16);
        outfile << seg.epgeominfo[0].dist;  // splineparameter (2D)
        outfile << " ";
        outfile.width(8);
        outfile.precision(16);
        outfile << seg.epgeominfo[1].edgenr;  // geometry dependent
        outfile << " ";
        outfile.width(12);
        outfile << seg.epgeominfo[1].dist;

        outfile << "\n";
      }


    outfile << "\n" << "\n";
    outfile << "#          X             Y             Z" << "\n";
    outfile << "points" << "\n";
    outfile << GetNP() << "\n";
    outfile.precision(16);
    outfile.setf (ios::fixed, ios::floatfield);
    outfile.setf (ios::showpoint);

    PointIndex pi;
    for (pi = PointIndex::BASE; 
         pi < GetNP()+PointIndex::BASE; pi++)
      {
        outfile.width(22);
        outfile << (*this)[pi](0)/scale << "  ";
        outfile.width(22);
        outfile << (*this)[pi](1)/scale << "  ";
        outfile.width(22);
        outfile << (*this)[pi](2)/scale << "\n";
      }

    outfile << "\n" << "\n";
    outfile << "#          pnum             index" << "\n";
    outfile << "pointelements" << "\n";
    outfile << pointelements.Size() << "\n";

    for (i = 0; i < pointelements.Size(); i++)
      {
        outfile.width(8);
        outfile << pointelements[i].pnum << "  ";
        outfile.width(8);
        outfile << pointelements[i].index << "\n";
      }

    if (ident -> GetMaxNr() > 0)
      {
        outfile << "identifications\n";
        NgArray<INDEX_2> identpairs;
        int cnt = 0;
        for (i = 1; i <= ident -> GetMaxNr(); i++)
          {
            ident -> GetPairs (i, identpairs);
            cnt += identpairs.Size();
          }
        outfile << cnt << "\n";
        for (i = 1; i <= ident -> GetMaxNr(); i++)
          {
            ident -> GetPairs (i, identpairs);
            for (j = 1; j <= identpairs.Size(); j++)
              {
                outfile.width (8);
                outfile << identpairs.Get(j).I1();
                outfile.width (8);
                outfile << identpairs.Get(j).I2();
                outfile.width (8);
                outfile << i << "\n";
              }
          }

        outfile << "identificationtypes\n";
        outfile << ident -> GetMaxNr() << "\n";
        for (i = 1; i <= ident -> GetMaxNr(); i++)
          {
            int type = ident -> GetType(i);
            outfile << " " << type;
          }
        outfile << "\n";
      }

    int cntmat = 0;
    for (i = 1; i <= materials.Size(); i++)
      if (materials.Get(i) && materials.Get(i)->length())
        cntmat++;

    if (cntmat)
      {
        outfile << "materials" << endl;
        outfile << cntmat << endl;
        for (i = 1; i <= materials.Size(); i++)
          if (materials.Get(i) && materials.Get(i)->length())
            outfile << i << " " << *materials.Get(i) << endl;
      }


    int cntbcnames = 0;
    for ( int ii = 0; ii < bcnames.Size(); ii++ )
      if ( bcnames[ii] ) cntbcnames++;

    if ( cntbcnames )
      {
        outfile << "\n\nbcnames" << endl << bcnames.Size() << endl;
        for ( i = 0; i < bcnames.Size(); i++ )
          outfile << i+1 << "\t" << GetBCName(i) << endl;
        outfile << endl << endl;
      }
    int cntcd2names = 0;
    for (int ii = 0; ii<cd2names.Size(); ii++)
      if(cd2names[ii]) cntcd2names++;

    if(cntcd2names)
      {
	outfile << "\n\ncd2names" << endl << cd2names.Size() << endl;
	for (i=0; i<cd2names.Size(); i++)
	  outfile << i+1 << "\t" << GetCD2Name(i) << endl;
	outfile << endl << endl;
      }

    int cntcd3names = 0;
    for (int ii = 0; ii<cd3names.Size(); ii++)
      if(cd3names[ii]) cntcd3names++;

    if(cntcd3names)
      {
	outfile << "\n\ncd3names" << endl << cd3names.Size() << endl;
	for (i=0; i<cd3names.Size(); i++)
	  outfile << i+1 << "\t" << GetCD3Name(i) << endl;
	outfile << endl << endl;
      }

    /*
      if ( GetDimension() == 2 )
      {
      for (i = 1; i <= GetNSeg(); i++)
      {
      const Segment & seg = LineSegment (i);
      if ( ! bcprops.Contains(seg.si) && seg.GetBCName() != "" )
      {
      bcprops.Append(seg.si);
      cntbcnames++;
      }
      }
      }
      else
      {
      for (sei = 0; sei < GetNSE(); sei++)
      {
      if ((*this)[sei].GetIndex())
      {
      int bcp = GetFaceDescriptor((*this)[sei].GetIndex ()).BCProperty();
      string name = GetFaceDescriptor((*this)[sei].GetIndex ()).BCName();
      if ( !bcprops.Contains(bcp) &&
      name != "" )
      {
      bcprops.Append(bcp);
      cntbcnames++;
      }
      }
      }
      }

      bcprops.SetSize(0);
      if ( cntbcnames )
      {
      outfile << "\nbcnames" << endl << cntbcnames << endl;
      if ( GetDimension() == 2 )
      {
      for (i = 1; i <= GetNSeg(); i++)
      {
      const Segment & seg = LineSegment (i);
      if ( ! bcprops.Contains(seg.si) && seg.GetBCName() != "" )
      {
      bcprops.Append(seg.si);
      outfile << seg.si << "\t" << seg.GetBCName() << endl;
      }
      }
      }
      else
      {
      for (sei = 0; sei < GetNSE(); sei++)
      {
      if ((*this)[sei].GetIndex())
      {
      int bcp = GetFaceDescriptor((*this)[sei].GetIndex ()).BCProperty();
      string name = GetFaceDescriptor((*this)[sei].GetIndex ()).BCName();
      if ( !bcprops.Contains(bcp) &&
      name != "" )
      {
      bcprops.Append(bcp);
      outfile << bcp << "\t" << name << endl;
      }
      }
      }
      }
      outfile << endl << endl;
      }
    */

    int cnt_sing = 0;
    // for (PointIndex pi = points.Begin(); pi < points.End(); pi++)
    // if ((*this)[pi].Singularity()>=1.) cnt_sing++;
    for (auto & p : points)
      if (p.Singularity() >= 1.) cnt_sing++;
      
    if (cnt_sing)
      {
        outfile << "singular_points" << endl << cnt_sing << endl;
        // for (PointIndex pi = points.Begin(); pi < points.End(); pi++)
        for (PointIndex pi : points.Range())
          if ((*this)[pi].Singularity()>=1.) 
            outfile << int(pi) << "\t" << (*this)[pi].Singularity() << endl;
      }

    cnt_sing = 0;
    for (SegmentIndex si = 0; si < GetNSeg(); si++)
      if ( segments[si].singedge_left ) cnt_sing++;
    if (cnt_sing)
      {
        outfile << "singular_edge_left" << endl << cnt_sing << endl;
        for (SegmentIndex si = 0; si < GetNSeg(); si++)
          if ( segments[si].singedge_left )
            outfile << int(si) << "\t" << segments[si].singedge_left << endl;
      }

    cnt_sing = 0;
    for (SegmentIndex si = 0; si < GetNSeg(); si++)
      if ( segments[si].singedge_right ) cnt_sing++;
    if (cnt_sing)
      {
        outfile << "singular_edge_right" << endl << cnt_sing << endl;
        for (SegmentIndex si = 0; si < GetNSeg(); si++)
          if ( segments[si].singedge_right  )
            outfile << int(si) << "\t" << segments[si].singedge_right << endl;
      }


    cnt_sing = 0;
    for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
      if ( GetFaceDescriptor ((*this)[sei].GetIndex()).domin_singular) 
        cnt_sing++;

    if (cnt_sing)
      {
        outfile << "singular_face_inside" << endl << cnt_sing << endl;
        for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
          if ( GetFaceDescriptor ((*this)[sei].GetIndex()).domin_singular) 
            outfile << int(sei)  << "\t" << 
              GetFaceDescriptor ((*this)[sei].GetIndex()).domin_singular  << endl;
      }

    cnt_sing = 0;
    for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
      if ( GetFaceDescriptor ((*this)[sei].GetIndex()).domout_singular) cnt_sing++;
    if (cnt_sing)
      {
        outfile << "singular_face_outside" << endl << cnt_sing << endl;
        for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
          if ( GetFaceDescriptor ((*this)[sei].GetIndex()).domout_singular) 
            outfile << int(sei) << "\t" 
                    << GetFaceDescriptor ((*this)[sei].GetIndex()).domout_singular << endl;
      }


    // Philippose - 09/07/2009
    // Add mesh face colours to Netgen Vol file format
    // The colours are saved in RGB triplets
    int cnt_facedesc = GetNFD();
    if (cnt_facedesc)
    {
       outfile << endl << endl << "#   Surfnr     Red     Green     Blue" << endl;
       outfile << "face_colours" << endl << cnt_facedesc << endl;

       outfile.precision(8);
       outfile.setf(ios::fixed, ios::floatfield);
       outfile.setf(ios::showpoint);

       for(i = 1; i <= cnt_facedesc; i++)
       {
          outfile.width(8);
          outfile << GetFaceDescriptor(i).SurfNr()+1 << " ";
          outfile.width(12);
          outfile << GetFaceDescriptor(i).SurfColour()[0] << " ";
          outfile.width(12);
          outfile << GetFaceDescriptor(i).SurfColour()[1] << " ";
          outfile.width(12);
          outfile << GetFaceDescriptor(i).SurfColour()[2];
          outfile << endl;
       }
    }

    outfile << endl << endl << "endmesh" << endl << endl;
    if (geometry)
      geometry -> SaveToMeshFile (outfile);
  }



  void Mesh :: Load (const string & filename)
  {
    cout << "filename = " << filename << endl;

    if (filename.find(".vol.bin") != string::npos)
    {
        BinaryInArchive in(filename);
        in & (*this);
        return;
    }

    istream * infile = NULL;

    if (filename.find(".vol.gz") != string::npos)
      infile = new igzstream (filename.c_str());
    else
      infile = new ifstream (filename.c_str());

    // ifstream infile(filename.c_str());
    if (! (infile -> good()) )
      throw NgException ("mesh file not found");

    Load(*infile);
    delete infile;
  }



  // Reads mandatory integer and optional string token from input stream
  // used for parsing bcnames, cd2names etc.
  void ReadNumberAndName( istream & infile, int & i, string & s )
  {
    string line;
    std::istringstream iline;

    bool empty_line = true;

    while(empty_line && infile)
      {
        std::getline(infile, line);
        iline = std::istringstream{line};
        iline >> i;

        if(iline)
            empty_line = false;

        iline >> s;
      }

    if(!infile)
        throw Exception("Reached end of file while parsing");
  }

  void Mesh :: Load (istream & infile)
  {
    static Timer timer("Mesh::Load"); RegionTimer rt(timer);
    if (! (infile.good()) )
      {
        cout << "cannot load mesh" << endl;
        throw NgException ("mesh file not found");
      }

    int rank = GetCommunicator().Rank();
    int ntasks = GetCommunicator().Size();
    
    char str[100];
    int i, n;

    double scale = 1;  // globflags.GetNumFlag ("scale", 1);
    int inverttets = 0;  // globflags.GetDefineFlag ("inverttets");
    int invertsurf = 0;  // globflags.GetDefineFlag ("invertsurfacemesh");


    facedecoding.SetSize(0);

    bool endmesh = false;
    

    while (infile.good() && !endmesh)
      {
        infile >> str;
        if (strcmp (str, "dimension") == 0)
          {
            infile >> dimension;
          }

        if (strcmp (str, "geomtype") == 0)
          {
            int hi;
            infile >> hi;
            geomtype = GEOM_TYPE(hi);
          }


        if (strcmp (str, "surfaceelements") == 0 || strcmp (str, "surfaceelementsgi")==0 || strcmp (str, "surfaceelementsuv") == 0)
          {
            static Timer t1("read surface elements"); RegionTimer rt1(t1);
            infile >> n;
            PrintMessage (3, n, " surface elements");

	    bool geominfo = strcmp (str, "surfaceelementsgi") == 0;
	    bool uv = strcmp (str, "surfaceelementsuv") == 0;


            for (i = 1; i <= n; i++)
              {
                int surfnr, bcp, domin, domout, nep, faceind = 0;

                infile >> surfnr >> bcp >> domin >> domout;
                surfnr--;

		bool invert_el = false;
		/*
		if (domin == 0) 
		  {
		    invert_el = true;
		    Swap (domin, domout);
		  }
		*/
		
                for (int j = 1; j <= facedecoding.Size(); j++)
                  if (GetFaceDescriptor(j).SurfNr() == surfnr &&
                      GetFaceDescriptor(j).BCProperty() == bcp &&
                      GetFaceDescriptor(j).DomainIn() == domin &&
                      GetFaceDescriptor(j).DomainOut() == domout)
                    faceind = j;

		// if (facedecoding.Size()) faceind = 1;   // for timing 

                if (!faceind)
                  {
                    faceind = AddFaceDescriptor (FaceDescriptor(surfnr, domin, domout, 0));
                    GetFaceDescriptor(faceind).SetBCProperty (bcp);
                  }

                infile >> nep;
                if (!nep) nep = 3;

                Element2d tri(nep);
                tri.SetIndex(faceind);

                for (int j = 1; j <= nep; j++)
                  infile >> tri.PNum(j);

                if (geominfo)
                  for (int j = 1; j <= nep; j++)
                    infile >> tri.GeomInfoPi(j).trignum;

                if (uv)
                  for (int j = 1; j <= nep; j++)
                    infile >> tri.GeomInfoPi(j).u >> tri.GeomInfoPi(j).v;
		
                if (invertsurf) tri.Invert();
		if (invert_el) tri.Invert();

		AddSurfaceElement (tri);
              }
          }

        if (strcmp (str, "volumeelements") == 0)
          {
            static Timer t1("read volume elements"); RegionTimer rt1(t1);
            infile >> n;
            PrintMessage (3, n, " volume elements");
            for (i = 1; i <= n; i++)
              {
                Element el(TET);
                int hi, nep;
                infile >> hi;
                if (hi == 0) hi = 1;
                el.SetIndex(hi);
                infile >> nep;
                el.SetNP(nep);
                el.SetCurved (nep != 4);
                for (int j = 0; j < nep; j++)
                  infile >> (int&)(el[j]);

                if (inverttets)
                  el.Invert();

		AddVolumeElement (el);
              }
          }


        if (strcmp (str, "edgesegments") == 0)
          {
            static Timer t1("read edge segments"); RegionTimer rt1(t1);
            infile >> n;
            for (i = 1; i <= n; i++)
              {
                Segment seg;
                int hi;
                infile >> seg.si >> hi >> seg[0] >> seg[1];
                AddSegment (seg);
              }
          }



        if (strcmp (str, "edgesegmentsgi") == 0)
          {
            static Timer t1("read edge segmentsgi"); RegionTimer rt1(t1);
            infile >> n;
            for (i = 1; i <= n; i++)
              {
                Segment seg;
                int hi;
                infile >> seg.si >> hi >> seg[0] >> seg[1]
                       >> seg.geominfo[0].trignum
                       >> seg.geominfo[1].trignum;
                AddSegment (seg);
              }
          }

        if (strcmp (str, "edgesegmentsgi2") == 0)
          {
            static Timer t1("read edge segmentsgi2"); RegionTimer rt1(t1);
            int a; 
            infile >> a;
            n=a; 

            PrintMessage (3, n, " curve elements");

            for (i = 1; i <= n; i++)
              {
                Segment seg;
                int hi;
                infile >> seg.si >> hi >> seg[0] >> seg[1]
                       >> seg.geominfo[0].trignum
                       >> seg.geominfo[1].trignum
                       >> seg.surfnr1 >> seg.surfnr2
                       >> seg.edgenr
                       >> seg.epgeominfo[0].dist
                       >> seg.epgeominfo[1].edgenr
                       >> seg.epgeominfo[1].dist;

                seg.epgeominfo[0].edgenr = seg.epgeominfo[1].edgenr;

                seg.domin = seg.surfnr1;
                seg.domout = seg.surfnr2;

                seg.surfnr1--;
                seg.surfnr2--;

                AddSegment (seg);
              }
          }

        if (strcmp (str, "points") == 0)
          {
            static Timer t1("read points"); RegionTimer rt1(t1);
            infile >> n;
            PrintMessage (3, n, " points");
            for (i = 1; i <= n; i++)
              {
                Point3d p;
                infile >> p.X() >> p.Y() >> p.Z();
                p.X() *= scale;
                p.Y() *= scale;
                p.Z() *= scale;
                AddPoint (p);
              }
	    PrintMessage (3, n, " points done");
          }

        if (strcmp (str, "pointelements") == 0)
          {
            static Timer t1("read point elements"); RegionTimer rt1(t1);
            infile >> n;
            PrintMessage (3, n, " pointelements");
            for (i = 1; i <= n; i++)
              {
                Element0d el;
                infile >> el.pnum >> el.index;
                pointelements.Append (el);
              }
	    PrintMessage (3, n, " pointelements done");
          }

        if (strcmp (str, "identifications") == 0)
          {
            infile >> n;
            PrintMessage (3, n, " identifications");
            for (i = 1; i <= n; i++)
              {
                PointIndex pi1, pi2;
                int ind;
                infile >> pi1 >> pi2 >> ind;
                ident -> Add (pi1, pi2, ind);
              }
          }

        if (strcmp (str, "identificationtypes") == 0)
          {
            infile >> n;
            PrintMessage (3, n, " identificationtypes");
            for (i = 1; i <= n; i++)
              {
                int type;
                infile >> type;
                ident -> SetType(i,Identifications::ID_TYPE(type));
              }
          }

        if (strcmp (str, "materials") == 0)
          {
            infile >> n;
            for ( auto i : Range(n) )
              {
                int nr;
                string mat;
                ReadNumberAndName( infile, nr, mat );
                SetMaterial (nr, mat.c_str());
              }
          }

        if ( strcmp (str, "bcnames" ) == 0 )
          {
            infile >> n;
            Array<int> bcnrs(n);
            SetNBCNames(n);
            for ( auto i : Range(n) )
              {
                string nextbcname;
                ReadNumberAndName( infile, bcnrs[i], nextbcname );
                bcnames[bcnrs[i]-1] = new string(nextbcname);
              }

            if ( GetDimension() == 3 )
              {
                for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
                  {
                    if ((*this)[sei].GetIndex())
                      {
                        int bcp = GetFaceDescriptor((*this)[sei].GetIndex ()).BCProperty();
                        if ( bcp <= n )
                          GetFaceDescriptor((*this)[sei].GetIndex ()).SetBCName(bcnames[bcp-1]);
                        else
                          GetFaceDescriptor((*this)[sei].GetIndex ()).SetBCName(0);

                      }
                  }

              }
          }

	if ( strcmp (str, "cd2names" ) == 0)
	  {
	    infile >> n;
	    Array<int> cd2nrs(n);
	    SetNCD2Names(n);
            for ( auto i : Range(n) )
              {
                string nextcd2name;
                ReadNumberAndName( infile, cd2nrs[i], nextcd2name );
                cd2names[cd2nrs[i]-1] = new string(nextcd2name);
              }
	    if (GetDimension() < 2)
	      {
		throw NgException("co dim 2 elements not implemented for dimension < 2");
	      }
	  }

        if ( strcmp (str, "cd3names" ) == 0)
	  {
	    infile >> n;
	    Array<int> cd3nrs(n);
	    SetNCD3Names(n);
	    for( auto i : Range(n) )
	      {
		string nextcd3name;
                ReadNumberAndName( infile, cd3nrs[i], nextcd3name );
		infile >> cd3nrs[i-1] >> nextcd3name;
		cd3names[cd3nrs[i-1]-1] = new string(nextcd3name);
	      }
	    if (GetDimension() < 3)
	      {
		throw NgException("co dim 3 elements not implemented for dimension < 3");
	      }
	  }

        if (strcmp (str, "singular_points") == 0)
          {
            infile >> n;
            for (i = 1; i <= n; i++)
              {
                PointIndex pi;
                double s; 
                infile >> pi;
                infile >> s; 
                (*this)[pi].Singularity (s);
              }
          }

        if (strcmp (str, "singular_edge_left") == 0)
          {
            infile >> n;
            for (i = 1; i <= n; i++)
              {
                SegmentIndex si;
                double s; 
                infile >> si;
                infile >> s; 
                (*this)[si].singedge_left = s;
              }
          }
        if (strcmp (str, "singular_edge_right") == 0)
          {
            infile >> n;
            for (i = 1; i <= n; i++)
              {
                SegmentIndex si;
                double s; 
                infile >> si;
                infile >> s; 
                (*this)[si].singedge_right = s;
              }
          }

        if (strcmp (str, "singular_face_inside") == 0)
          {
            infile >> n;
            for (i = 1; i <= n; i++)
              {
                SurfaceElementIndex sei;
                double s; 
                infile >> sei;
                infile >> s; 
                GetFaceDescriptor((*this)[sei].GetIndex()).domin_singular = s;
              }
          }

        if (strcmp (str, "singular_face_outside") == 0)
          {
            infile >> n;
            for (i = 1; i <= n; i++)
              {
                SurfaceElementIndex sei;
                double s; 
                infile >> sei;
                infile >> s; 
                GetFaceDescriptor((*this)[sei].GetIndex()).domout_singular = s;
              }
          }

        // Philippose - 09/07/2009
        // Add mesh face colours to Netgen Vol file format
        // The colours are read in as RGB triplets
        if (strcmp (str, "face_colours") == 0)
        {
           int cnt_facedesc = GetNFD();
           infile >> n;
           if(n == cnt_facedesc)
           {
              for(i = 1; i <= n; i++)
              {
                 int surfnr = 0;
                 Vec3d surfcolour(0.0,1.0,0.0);

                 infile >> surfnr 
                        >> surfcolour.X() 
                        >> surfcolour.Y() 
                        >> surfcolour.Z();

                 surfnr--;

                 if(surfnr > 0) 
                 {
                    for(int facedesc = 1; facedesc <= cnt_facedesc; facedesc++)
                    {
                       if(surfnr == GetFaceDescriptor(facedesc).SurfNr())
                       {
                          GetFaceDescriptor(facedesc).SetSurfColour(surfcolour);
                       }
                    }
                 }
              }
           }
        }

        if (strcmp (str, "endmesh") == 0)
          endmesh = true;



        strcpy (str, "");
      }




    CalcSurfacesOfNode ();
 
    if (ntasks == 1) // sequential run only
      {
	topology.Update();
	clusters -> Update();
      }

    SetNextMajorTimeStamp();
    //  PrintMemInfo (cout);
  }


  void Mesh :: DoArchive (Archive & archive)
  {
    static Timer t("Mesh::Archive"); RegionTimer r(t);

#ifdef PARALLEL
    auto comm = GetCommunicator();
    if (archive.IsParallel() && comm.Size() > 1)
      { // parallel pickling supported only for output archives
        if (comm.Rank() == 0)
          archive & dimension;

        auto rank = comm.Rank();
        
        auto & partop = GetParallelTopology();
        
        // global enumration of points:
        // not used now, but will be needed for refined meshes
        // GridFunciton pickling is not compatible, now
        // should go to paralleltopology
        
        
        
        // merge points
        Array<PointIndex, PointIndex> globnum(points.Size());
        PointIndex maxglob = -1;
        for (auto pi : Range(points))
          {
            globnum[pi] = partop.GetGlobalPNum(pi);
            // globnum[pi] = global_pnums[pi];
            maxglob = max(globnum[pi], maxglob);
          }
        
        maxglob = comm.AllReduce (maxglob, MPI_MAX);
        int numglob = maxglob+1-PointIndex::BASE;
        if (comm.Rank() > 0)
          {
            comm.Send (globnum, 0, 200);
            comm.Send (points, 0, 200);
          }
        else
          {
            Array<PointIndex, PointIndex> globnumi;
            Array<MeshPoint, PointIndex> pointsi;
            Array<MeshPoint, PointIndex> globpoints(numglob);
            for (int j = 1; j < comm.Size(); j++)
              {
                comm.Recv (globnumi, j, 200);
                comm.Recv (pointsi, j, 200);
                for (auto i : Range(globnumi))
                  globpoints[globnumi[i]] = pointsi[i];
              }
            archive & globpoints;
          }

        
        // sending surface elements
        auto copy_el2d  (surfelements);
        for (auto & el : copy_el2d)
          for (auto & pi : el.PNums())
            pi = globnum[pi];

        if (comm.Rank() > 0)
          comm.Send(copy_el2d, 0, 200);
        else
          {
            Array<Element2d, SurfaceElementIndex> el2di;
            for (int j = 1; j < comm.Size(); j++)
              {
                comm.Recv(el2di, j, 200);
                for (auto & el : el2di)
                  copy_el2d += el;
              }
            archive & copy_el2d;
          }


        // sending volume elements
        auto copy_el3d  (volelements);
        for (auto & el : copy_el3d)
          for (auto & pi : el.PNums())
            pi = globnum[pi];

        if (comm.Rank() > 0)
          comm.Send(copy_el3d, 0, 200);
        else
          {
            Array<Element, ElementIndex> el3di;
            for (int j = 1; j < comm.Size(); j++)
              {
                comm.Recv(el3di, j, 200);
                for (auto & el : el3di)
                  copy_el3d += el;
              }
            archive & copy_el3d;
          }


        // sending 1D elements
        auto copy_el1d  (segments);
        for (auto & el : copy_el1d)
          for (auto & pi : el.pnums)
            if (pi != PointIndex(PointIndex::INVALID))
              pi = globnum[pi];

        if (comm.Rank() > 0)
          comm.Send(copy_el1d, 0, 200);
        else
          {
            Array<Segment, SegmentIndex> el1di;
            for (int j = 1; j < comm.Size(); j++)
              {
                comm.Recv(el1di, j, 200);
                for (auto & el : el1di)
                  copy_el1d += el;
              }
            archive & copy_el1d;
          }

        if (comm.Rank() == 0)
          {
            archive & facedecoding;
            archive & materials & bcnames & cd2names & cd3names;
            auto mynv = numglob;
            archive & mynv;   // numvertices;
            archive & *ident;
            
            archive.Shallow(geometry);
            archive & *curvedelems;
          }
        
        if (comm.Rank() == 0)
          return;
      }
#endif
    
    
    archive & dimension;
    archive & points;
    archive & surfelements;
    archive & volelements;
    archive & segments;
    archive & facedecoding;
    archive & materials & bcnames & cd2names & cd3names;
    archive & numvertices;

    archive & *ident;

    archive.Shallow(geometry);
    archive & *curvedelems;
    
    if (archive.Input())
      {
	int rank = GetCommunicator().Rank();
	int ntasks = GetCommunicator().Size();
	
        RebuildSurfaceElementLists();
        
        CalcSurfacesOfNode ();
        if (ntasks == 1) // sequential run only
          {
            topology.Update();
            clusters -> Update();
          }
        SetNextMajorTimeStamp();
      }
  }


  void Mesh :: Merge (const string & filename, const int surfindex_offset)
  {
    ifstream infile(filename.c_str());
    if (!infile.good())
      throw NgException ("mesh file not found");

    Merge(infile,surfindex_offset);

  }



  void Mesh :: Merge (istream & infile, const int surfindex_offset)
  {
    char str[100];
    int i, n;


    int inverttets = 0;  // globflags.GetDefineFlag ("inverttets");

    int oldnp = GetNP();
    int oldne = GetNSeg();
    int oldnd = GetNDomains();

    for(SurfaceElementIndex si = 0; si < GetNSE(); si++)
      for(int j=1; j<=(*this)[si].GetNP(); j++) (*this)[si].GeomInfoPi(j).trignum = -1;

    int max_surfnr = 0;
    for (i = 1; i <= GetNFD(); i++)
      max_surfnr = max2 (max_surfnr, GetFaceDescriptor(i).SurfNr());
    max_surfnr++;

    if(max_surfnr < surfindex_offset) max_surfnr = surfindex_offset;


    bool endmesh = false;

    while (infile.good() && !endmesh)
      {
        infile >> str;

        if (strcmp (str, "surfaceelementsgi") == 0 || strcmp (str, "surfaceelements") == 0)
          {
            infile >> n;
            PrintMessage (3, n, " surface elements");
            for (i = 1; i <= n; i++)
              {
                int j;
                int surfnr, bcp, domin, domout, nep, faceind = 0;
                infile >> surfnr >> bcp >> domin >> domout;

                surfnr--;

                if(domin > 0) domin += oldnd;
                if(domout > 0) domout += oldnd;
                surfnr += max_surfnr;


                for (j = 1; j <= facedecoding.Size(); j++)
                  if (GetFaceDescriptor(j).SurfNr() == surfnr &&
                      GetFaceDescriptor(j).BCProperty() == bcp &&
                      GetFaceDescriptor(j).DomainIn() == domin &&
                      GetFaceDescriptor(j).DomainOut() == domout)
                    faceind = j;

                if (!faceind)
                  {
                    faceind = AddFaceDescriptor (FaceDescriptor(surfnr, domin, domout, 0));
                    if(GetDimension() == 2) bcp++;
                    GetFaceDescriptor(faceind).SetBCProperty (bcp);
                  }

                infile >> nep;
                if (!nep) nep = 3;

                Element2d tri(nep);
                tri.SetIndex(faceind);

                for (j = 1; j <= nep; j++)
                  {
                    infile >> tri.PNum(j);
                    tri.PNum(j) = tri.PNum(j) + oldnp;
                  }


                if (strcmp (str, "surfaceelementsgi") == 0)
                  for (j = 1; j <= nep; j++)
                    {
                      infile >> tri.GeomInfoPi(j).trignum;
                      tri.GeomInfoPi(j).trignum = -1;
                    }

                AddSurfaceElement (tri);
              }
          }


        if (strcmp (str, "edgesegments") == 0)
          {
            infile >> n;
            for (i = 1; i <= n; i++)
              {
                Segment seg;
                int hi;
                infile >> seg.si >> hi >> seg[0] >> seg[1];
                seg[0] = seg[0] + oldnp;
                seg[1] = seg[1] + oldnp;
                AddSegment (seg);
              }
          }



        if (strcmp (str, "edgesegmentsgi") == 0)
          {
            infile >> n;
            for (i = 1; i <= n; i++)
              {
                Segment seg;
                int hi;
                infile >> seg.si >> hi >> seg[0] >> seg[1]
                       >> seg.geominfo[0].trignum
                       >> seg.geominfo[1].trignum;
                seg[0] = seg[0] + oldnp;
                seg[1] = seg[1] + oldnp;
                AddSegment (seg);
              }
          }
        if (strcmp (str, "edgesegmentsgi2") == 0)
          {
            infile >> n;
            PrintMessage (3, n, " curve elements");

            for (i = 1; i <= n; i++)
              {
                Segment seg;
                int hi;
                infile >> seg.si >> hi >> seg[0] >> seg[1]
                       >> seg.geominfo[0].trignum
                       >> seg.geominfo[1].trignum
                       >> seg.surfnr1 >> seg.surfnr2
                       >> seg.edgenr
                       >> seg.epgeominfo[0].dist
                       >> seg.epgeominfo[1].edgenr
                       >> seg.epgeominfo[1].dist;
                seg.epgeominfo[0].edgenr = seg.epgeominfo[1].edgenr;

                seg.surfnr1--;
                seg.surfnr2--;

                if(seg.surfnr1 >= 0)  seg.surfnr1 = seg.surfnr1 + max_surfnr;
                if(seg.surfnr2 >= 0)  seg.surfnr2 = seg.surfnr2 + max_surfnr;
                seg[0] = seg[0] +oldnp;
                seg[1] = seg[1] +oldnp;
		*testout << "old edgenr: " << seg.edgenr << endl;
                seg.edgenr = seg.edgenr + oldne;
		*testout << "new edgenr: " << seg.edgenr << endl;
                seg.epgeominfo[1].edgenr = seg.epgeominfo[1].edgenr + oldne;

                AddSegment (seg);
              }
          }

        if (strcmp (str, "volumeelements") == 0)
          {
            infile >> n;
            PrintMessage (3, n, " volume elements");
            for (i = 1; i <= n; i++)
              {
                Element el(TET);
                int hi, nep;
                infile >> hi;
                if (hi == 0) hi = 1;
                el.SetIndex(hi+oldnd);
                infile >> nep;
                el.SetNP(nep);

                for (int j = 0; j < nep; j++)
                  {
                    infile >> (int&)(el[j]);
                    el[j] = el[j]+oldnp;
                  }

                if (inverttets)
                  el.Invert();

                AddVolumeElement (el);
              }
          }


        if (strcmp (str, "points") == 0)
          {
            infile >> n;
            PrintMessage (3, n, " points");
            for (i = 1; i <= n; i++)
              {
                Point3d p;
                infile >> p.X() >> p.Y() >> p.Z();
                AddPoint (p);
              }
          }


        if (strcmp (str, "endmesh") == 0)
          {
            endmesh = true;
          }


        if (strcmp (str, "materials") == 0)
          {
            infile >> n;
            for (i = 1; i <= n; i++)
              {
                int nr;
                string mat;
                infile >> nr >> mat;
                SetMaterial (nr+oldnd, mat.c_str());
              }
          }


        strcpy (str, "");
      }

    CalcSurfacesOfNode ();

    topology.Update();
    clusters -> Update();

    SetNextMajorTimeStamp();
  }










  bool Mesh :: TestOk () const
  {
    for (ElementIndex ei = 0; ei < volelements.Size(); ei++)
      {
        for (int j = 0; j < 4; j++)
          if ( (*this)[ei][j] <= PointIndex::BASE-1)
            {
              (*testout) << "El " << ei << " has 0 nodes: ";
              for (int k = 0; k < 4; k++)
                (*testout) << (*this)[ei][k];
              break;
            }
      }
    CheckMesh3D (*this);
    return 1;
  }

  void Mesh :: SetAllocSize(int nnodes, int nsegs, int nsel, int nel)
  {
    points.SetAllocSize(nnodes);
    segments.SetAllocSize(nsegs);
    surfelements.SetAllocSize(nsel);
    volelements.SetAllocSize(nel);
  }

  void Mesh :: BuildBoundaryEdges(bool rebuild)
  {
    static Timer t("Mesh::BuildBoundaryEdges"); RegionTimer reg(t);
    
    if(!rebuild && boundaryedges)
      return;

    boundaryedges = make_unique<INDEX_2_CLOSED_HASHTABLE<int>>
      (3 * (GetNSE() + GetNOpenElements()) + GetNSeg() + 1);


    for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
      {
        const Element2d & sel = surfelements[sei];
        if (sel.IsDeleted()) continue;

        // int si = sel.GetIndex();

        if (sel.GetNP() <= 4)
          for (int j = 0; j < sel.GetNP(); j++)
            {
              INDEX_2 i2;
              i2.I1() = sel.PNumMod(j+1);
              i2.I2() = sel.PNumMod(j+2);
              i2.Sort();
              boundaryedges->Set (i2, 1);
            }
        else if (sel.GetType()==TRIG6)
          {
            for (int j = 0; j < 3; j++)
              {
                INDEX_2 i2;
                i2.I1() = sel[j];
                i2.I2() = sel[(j+1)%3];
                i2.Sort();
                boundaryedges->Set (i2, 1);
              }
          }
        else 
          cerr << "illegal element for buildboundaryedges" << endl;
      }


    for (int i = 0; i < openelements.Size(); i++)
      {
        const Element2d & sel = openelements[i];
        for (int j = 0; j < sel.GetNP(); j++)
          {
            INDEX_2 i2;
            i2.I1() = sel.PNumMod(j+1);
            i2.I2() = sel.PNumMod(j+2);
            i2.Sort();
            boundaryedges->Set (i2, 1);

            points[sel[j]].SetType(FIXEDPOINT);
          }
      }

    for (int i = 0; i < GetNSeg(); i++)
      {
        const Segment & seg = segments[i];
        INDEX_2 i2(seg[0], seg[1]);
        i2.Sort();

        boundaryedges -> Set (i2, 2);
        //segmentht -> Set (i2, i);
      }


  }

  void Mesh :: CalcSurfacesOfNode ()
  {
    static Timer t("Mesh::CalcSurfacesOfNode"); RegionTimer reg (t);
    static Timer tn2se("Mesh::CalcSurfacesOfNode - surf on node");     
    static Timer tht("Mesh::CalcSurfacesOfNode - surfelementht"); 
    // surfacesonnode.SetSize (GetNP());
    TABLE<int,PointIndex::BASE> surfacesonnode(GetNP());

    // delete boundaryedges;
    // boundaryedges = NULL;
    boundaryedges = nullptr;

    // delete surfelementht;
    // surfelementht = nullptr;
    surfelementht = nullptr;
    // delete segmentht;

    /*
      surfelementht = new INDEX_3_HASHTABLE<int> (GetNSE()/4 + 1);
      segmentht = new INDEX_2_HASHTABLE<int> (GetNSeg() + 1);
    */

    if (dimension == 3)
      surfelementht = make_unique<INDEX_3_CLOSED_HASHTABLE<int>> (3*GetNSE() + 1);
    segmentht = make_unique<INDEX_2_CLOSED_HASHTABLE<int>> (3*GetNSeg() + 1);

    tn2se.Start();
    if (dimension == 3)
      /*
    for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
      {
        const Element2d & sel = surfelements[sei];
      */
      for (const Element2d & sel : surfelements)
        {
        if (sel.IsDeleted()) continue;

        int si = sel.GetIndex();

        /*
        for (int j = 0; j < sel.GetNP(); j++)
          {
            PointIndex pi = sel[j];
        */
        for (PointIndex pi : sel.PNums())
          {
            if (!surfacesonnode[pi].Contains(si))
              surfacesonnode.Add (pi, si);
            /*
            bool found = 0;
            for (int k = 0; k < surfacesonnode[pi].Size(); k++)
              if (surfacesonnode[pi][k] == si)
                {
                  found = 1;
                  break;
                }

            if (!found)
              surfacesonnode.Add (pi, si);
            */
          }
      }
    /*
      for (sei = 0; sei < GetNSE(); sei++)
      {
      const Element2d & sel = surfelements[sei];
      if (sel.IsDeleted()) continue;

      INDEX_3 i3;
      i3.I1() = sel.PNum(1);
      i3.I2() = sel.PNum(2);
      i3.I3() = sel.PNum(3);
      i3.Sort();
      surfelementht -> PrepareSet (i3);
      }

      surfelementht -> AllocateElements();
    */
    tn2se.Stop();
    
    tht.Start();
    if (dimension==3)
    for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
      {
        const Element2d & sel = surfelements[sei];
        if (sel.IsDeleted()) continue;

        INDEX_3 i3;
        i3.I1() = sel.PNum(1);
        i3.I2() = sel.PNum(2);
        i3.I3() = sel.PNum(3);
        i3.Sort();
        surfelementht -> Set (i3, sei);   // war das wichtig ???    sel.GetIndex());
      }
    tht.Stop();
    
    // int np = GetNP();

    if (dimension == 3)
      {
        static Timer t("Mesh::CalcSurfacesOfNode, pointloop"); RegionTimer reg (t);            
        /*
        for (PointIndex pi = points.Begin(); pi < points.End(); pi++)
          points[pi].SetType (INNERPOINT);
        */
        for (auto & p : points)
          p.SetType (INNERPOINT);
        
        if (GetNFD() == 0) 
          {
            for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
              {
                const Element2d & sel = surfelements[sei];
                if (sel.IsDeleted()) continue;
                for (int j = 0;  j < sel.GetNP(); j++)
                  {
                    PointIndex pi = SurfaceElement(sei)[j];
                    points[pi].SetType(FIXEDPOINT);
                  }
              }
          }
        else
          {
            for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
              {
                const Element2d & sel = surfelements[sei];
                if (sel.IsDeleted()) continue;
                for (int j = 0; j < sel.GetNP(); j++)
                  {
                    PointIndex pi = sel[j];
                    int ns = surfacesonnode[pi].Size();
                    if (ns == 1)
                      points[pi].SetType(SURFACEPOINT);
                    if (ns == 2)
                      points[pi].SetType(EDGEPOINT);
                    if (ns >= 3)
                      points[pi].SetType(FIXEDPOINT);
                  }      
              }
          }
      }

    /*
    for (int i = 0; i < segments.Size(); i++)
      {
        const Segment & seg = segments[i];
    */
    for (const Segment & seg : segments)
      {
        for (int j = 1; j <= 2; j++)
          {
            PointIndex hi = (j == 1) ? seg[0] : seg[1];
            if (points[hi].Type() == INNERPOINT ||
                points[hi].Type() == SURFACEPOINT)
              points[hi].SetType(EDGEPOINT);
          }
      }
    
    for (int i = 0; i < lockedpoints.Size(); i++)
      points[lockedpoints[i]].SetType(FIXEDPOINT);


    /*
      for (i = 0; i < openelements.Size(); i++)
      {
      const Element2d & sel = openelements[i];
      for (j = 0; j < sel.GetNP(); j++)
      {
      INDEX_2 i2;
      i2.I1() = sel.PNumMod(j+1);
      i2.I2() = sel.PNumMod(j+2);
      i2.Sort();
      boundaryedges->Set (i2, 1);

      points[sel[j]].SetType(FIXEDPOINT);
      }
      }
    */

    // eltyps.SetSize (GetNE());
    // eltyps = FREEELEMENT;

    for (int i = 0; i < GetNSeg(); i++)
      {
        const Segment & seg = segments[i];
        INDEX_2 i2(seg[0], seg[1]);
        i2.Sort();

        //boundaryedges -> Set (i2, 2);
        segmentht -> Set (i2, i);
      }
  }

  // NgBitArray base is PointIndex::BASE ... 
  void Mesh :: FixPoints (const NgBitArray & fixpoints)
  {
    if (fixpoints.Size() != GetNP())
      {
        cerr << "Mesh::FixPoints: sizes don't fit" << endl;
        return;
      }
    int np = GetNP();
    /*
    for (int i = 1; i <= np; i++)
      if (fixpoints.Test(i))
        {
          points.Elem(i).SetType (FIXEDPOINT);
        }
    */
    for (PointIndex pi : points.Range())
      if (fixpoints.Test(pi))
        points[pi].SetType(FIXEDPOINT);
  }


  void Mesh :: FindOpenElements (int dom)
  {
    static Timer t("Mesh::FindOpenElements"); RegionTimer reg (t);
    static Timer t_table("Mesh::FindOpenElements - build table"); 
    static Timer t_pointloop("Mesh::FindOpenElements - pointloop"); 

    int np = GetNP();
    int ne = GetNE();
    int nse = GetNSE();
    
    t_table.Start();

    auto elsonpoint = ngcore::CreateSortedTable<ElementIndex, PointIndex>( volelements.Range(),
           [&](auto & table, ElementIndex ei)
           {
             const Element & el = (*this)[ei];
             if (dom == 0 || dom == el.GetIndex())
               {
                 if (el.GetNP() == 4)
                   {
                     INDEX_4 i4(el[0], el[1], el[2], el[3]);
                     i4.Sort();
                     table.Add (PointIndex(i4.I1()), ei);
                     table.Add (PointIndex(i4.I2()), ei);
                   }
                 else
                   {
                     for (PointIndex pi : el.PNums())
                       table.Add(pi, ei);
                   }
               }
           }, GetNP());


    NgArray<int,PointIndex::BASE> numonpoint(np);
    /*
    numonpoint = 0;
    for (ElementIndex ei = 0; ei < ne; ei++)
      {
        const Element & el = (*this)[ei];
        if (dom == 0 || dom == el.GetIndex())
          {
            if (el.GetNP() == 4)
              {
                INDEX_4 i4(el[0], el[1], el[2], el[3]);
                i4.Sort();
                numonpoint[i4.I1()]++;
                numonpoint[i4.I2()]++;
              }
            else
              for (int j = 0; j < el.GetNP(); j++)
                numonpoint[el[j]]++;
          }
      }

    TABLE<ElementIndex,PointIndex::BASE> elsonpoint(numonpoint);
    for (ElementIndex ei = 0; ei < ne; ei++)
      {
        const Element & el = (*this)[ei];
        if (dom == 0 || dom == el.GetIndex())
          {
            if (el.GetNP() == 4)
              {
                INDEX_4 i4(el[0], el[1], el[2], el[3]);
                i4.Sort();
                elsonpoint.Add (i4.I1(), ei);
                elsonpoint.Add (i4.I2(), ei);
              }
            else
              for (int j = 0; j < el.GetNP(); j++)
                elsonpoint.Add (el[j], ei);
          }
      }
    */
    t_table.Stop();


    NgArray<bool, 1> hasface(GetNFD());

    for (int i = 1; i <= GetNFD(); i++)
      {
        int domin = GetFaceDescriptor(i).DomainIn();
        int domout = GetFaceDescriptor(i).DomainOut();
        hasface[i] = 
          ( dom == 0 && (domin != 0 || domout != 0) ) ||
          ( dom != 0 && (domin == dom || domout == dom) );
      }

    numonpoint = 0;
    for (SurfaceElementIndex sii = 0; sii < nse; sii++)
      {
        int ind = surfelements[sii].GetIndex();
        /*
          if (
          GetFaceDescriptor(ind).DomainIn() && 
          (dom == 0 || dom == GetFaceDescriptor(ind).DomainIn())
          ||
          GetFaceDescriptor(ind).DomainOut() && 
          (dom == 0 || dom == GetFaceDescriptor(ind).DomainOut())
          )
        */
        if (hasface[ind])
          {
            /*
              Element2d hel = surfelements[i];
              hel.NormalizeNumbering();	  
              numonpoint[hel[0]]++;
            */
            const Element2d & hel = surfelements[sii];
            int mini = 0;
            for (int j = 1; j < hel.GetNP(); j++)
              if (hel[j] < hel[mini])
                mini = j;
            numonpoint[hel[mini]]++;
          }
      }

    TABLE<SurfaceElementIndex,PointIndex::BASE> selsonpoint(numonpoint);
    for (SurfaceElementIndex sii = 0; sii < nse; sii++)
      {
        int ind = surfelements[sii].GetIndex();

        /*
          if (
          GetFaceDescriptor(ind).DomainIn() && 
          (dom == 0 || dom == GetFaceDescriptor(ind).DomainIn())
          ||
          GetFaceDescriptor(ind).DomainOut() && 
          (dom == 0 || dom == GetFaceDescriptor(ind).DomainOut())
          )
        */
        if (hasface[ind])
          {
            /*
              Element2d hel = surfelements[i];
              hel.NormalizeNumbering();	  
              selsonpoint.Add (hel[0], i);
            */
            const Element2d & hel = surfelements[sii];
            int mini = 0;
            for (int j = 1; j < hel.GetNP(); j++)
              if (hel[j] < hel[mini])
                mini = j;
            selsonpoint.Add (hel[mini], sii);
          }
      }


    // PointIndex pi;
    // SurfaceElementIndex sei;
    // Element2d hel;

    struct tval { int index; PointIndex p4; };
    openelements.SetSize(0);
    
    t_pointloop.Start();

    /*
    INDEX_3_CLOSED_HASHTABLE<tval> faceht(100);
    
    for (PointIndex pi : points.Range())
      if (selsonpoint[pi].Size()+elsonpoint[pi].Size())
        {
          faceht.SetSize (2 * selsonpoint[pi].Size() + 4 * elsonpoint[pi].Size());

          for (SurfaceElementIndex sei : selsonpoint[pi])
            {
              Element2d hel = SurfaceElement(sei);
              if (hel.GetType() == TRIG6) hel.SetType(TRIG);
              int ind = hel.GetIndex();	  

              if (GetFaceDescriptor(ind).DomainIn() && 
                  (dom == 0 || dom == GetFaceDescriptor(ind).DomainIn()) )
                {
                  hel.NormalizeNumbering();
                  if (hel.PNum(1) == pi)
                    {
                      INDEX_3 i3(hel[0], hel[1], hel[2]);
                      tval i2;
                      i2.index = GetFaceDescriptor(ind).DomainIn();
                      i2.p4 = (hel.GetNP() == 3)
                            ? PointIndex (PointIndex::INVALID)
                      : hel.PNum(4);
                      faceht.Set (i3, i2);
                    }
                }
              if (GetFaceDescriptor(ind).DomainOut() &&
                  (dom == 0 || dom == GetFaceDescriptor(ind).DomainOut()) )
                {
                  hel.Invert();
                  hel.NormalizeNumbering();
                  if (hel.PNum(1) == pi)
                    {
                      INDEX_3 i3(hel[0], hel[1], hel[2]);
                      tval i2;
                      i2.index = GetFaceDescriptor(ind).DomainOut();
                      i2.p4 = (hel.GetNP() == 3)
                        ? PointIndex (PointIndex::INVALID)
                        : hel.PNum(4);
                      faceht.Set (i3, i2);
                    }
                }
            }

          for (ElementIndex ei : elsonpoint[pi])
            {
              const Element & el = VolumeElement(ei);

              if (dom == 0 || el.GetIndex() == dom)
                {
                  for (int j = 1; j <= el.GetNFaces(); j++)
                    {
                      Element2d hel(TRIG);
                      el.GetFace (j, hel);
                      hel.Invert();
                      hel.NormalizeNumbering();

                      if (hel[0] == pi)
                        {
                          INDEX_3 i3(hel[0], hel[1], hel[2]);

                          if (faceht.Used (i3))
                            {
                              tval i2 = faceht.Get(i3);
                              if (i2.index == el.GetIndex())
                                {
                                  i2.index = PointIndex::BASE-1;
                                  faceht.Set (i3, i2);
                                }
                              else
                                {
                                  if (i2.index == 0)
                                    {
                                      PrintSysError ("more elements on face");
                                      (*testout)  << "more elements on face!!!" << endl;
                                      (*testout) << "el = " << el << endl;
                                      (*testout) << "hel = " << hel << endl;
                                      (*testout) << "face = " << i3 << endl;
                                      (*testout) << "points = " << endl;
                                      for (int jj = 1; jj <= 3; jj++)
                                        (*testout) << "p = " << Point(i3.I(jj)) << endl;
                                    }
                                }
                            }
                          else
                            {
                              hel.Invert();
                              hel.NormalizeNumbering();
                              INDEX_3 i3(hel[0], hel[1], hel[2]);
                              
                              tval i2;
                              i2.index = el.GetIndex();
                              i2.p4 = (hel.GetNP() == 3)
                                ? PointIndex (PointIndex::INVALID)
                                : hel[3];
                              faceht.Set (i3, i2);
                            }
                        }
                    }
                }
            }
          
          for (int i = 0; i < faceht.Size(); i++)
            if (faceht.UsedPos (i))
              {
                INDEX_3 i3;
                //INDEX_2 i2;
                tval i2;
                faceht.GetData (i, i3, i2);
                if (i2.index != PointIndex::BASE-1)
                  {
                    Element2d tri ( (i2.p4 == PointIndex::BASE-1) ? TRIG : QUAD);
                    for (int l = 0; l < 3; l++)
                      tri[l] = i3.I(l+1);
                    tri.PNum(4) = i2.p4;
                    tri.SetIndex (i2.index);
                    openelements.Append (tri);
                  }
              }
        }

    */

    size_t numtasks = 4*ngcore::TaskManager::GetNumThreads();
    Array<Array<Element2d>> thread_openelements(numtasks);
    ParallelJob
      ( [&](TaskInfo & ti)
      {
        auto myrange = points.Range().Split(ti.task_nr, ti.ntasks);
        INDEX_3_CLOSED_HASHTABLE<tval> faceht(100);        
        for (PointIndex pi : myrange)
          if (selsonpoint[pi].Size()+elsonpoint[pi].Size())
            {
              faceht.SetSize (2 * selsonpoint[pi].Size() + 4 * elsonpoint[pi].Size());

              for (SurfaceElementIndex sei : selsonpoint[pi])
                {
                  Element2d hel = SurfaceElement(sei);
                  if (hel.GetType() == TRIG6) hel.SetType(TRIG);
                  int ind = hel.GetIndex();	  

                  if (GetFaceDescriptor(ind).DomainIn() && 
                      (dom == 0 || dom == GetFaceDescriptor(ind).DomainIn()) )
                    {
                      hel.NormalizeNumbering();
                      if (hel.PNum(1) == pi)
                        {
                          INDEX_3 i3(hel[0], hel[1], hel[2]);
                          tval i2;
                          i2.index = GetFaceDescriptor(ind).DomainIn();
                          i2.p4 = (hel.GetNP() == 3)
                            ? PointIndex (PointIndex::INVALID)
                            : hel.PNum(4);
                          faceht.Set (i3, i2);
                        }
                    }
                  if (GetFaceDescriptor(ind).DomainOut() &&
                      (dom == 0 || dom == GetFaceDescriptor(ind).DomainOut()) )
                    {
                      hel.Invert();
                      hel.NormalizeNumbering();
                      if (hel.PNum(1) == pi)
                        {
                          INDEX_3 i3(hel[0], hel[1], hel[2]);
                          tval i2;
                          i2.index = GetFaceDescriptor(ind).DomainOut();
                          i2.p4 = (hel.GetNP() == 3)
                            ? PointIndex (PointIndex::INVALID)
                            : hel.PNum(4);
                          faceht.Set (i3, i2);
                        }
                    }
                }
              
              for (ElementIndex ei : elsonpoint[pi])
                {
                  const Element & el = VolumeElement(ei);
                  
                  if (dom == 0 || el.GetIndex() == dom)
                    {
                      for (int j = 1; j <= el.GetNFaces(); j++)
                        {
                          Element2d hel(TRIG);
                          el.GetFace (j, hel);
                          hel.Invert();
                          hel.NormalizeNumbering();
                          
                          if (hel[0] == pi)
                            {
                              INDEX_3 i3(hel[0], hel[1], hel[2]);
                              
                              if (faceht.Used (i3))
                                {
                                  tval i2 = faceht.Get(i3);
                                  if (i2.index == el.GetIndex())
                                    {
                                      i2.index = PointIndex::BASE-1;
                                      faceht.Set (i3, i2);
                                    }
                                  else
                                    {
                                      if (i2.index == 0)
                                        {
                                          PrintSysError ("more elements on face");
                                          (*testout)  << "more elements on face!!!" << endl;
                                          (*testout) << "el = " << el << endl;
                                          (*testout) << "hel = " << hel << endl;
                                          (*testout) << "face = " << i3 << endl;
                                          (*testout) << "points = " << endl;
                                          for (int jj = 1; jj <= 3; jj++)
                                            (*testout) << "p = " << Point(i3.I(jj)) << endl;
                                        }
                                    }
                                }
                              else
                                {
                                  hel.Invert();
                                  hel.NormalizeNumbering();
                                  INDEX_3 i3(hel[0], hel[1], hel[2]);
                                  
                                  tval i2;
                                  i2.index = el.GetIndex();
                                  i2.p4 = (hel.GetNP() == 3)
                                    ? PointIndex (PointIndex::INVALID)
                                    : hel[3];
                                  faceht.Set (i3, i2);
                                }
                            }
                        }
                    }
                }
              
              for (int i = 0; i < faceht.Size(); i++)
                if (faceht.UsedPos (i))
                  {
                    INDEX_3 i3;
                    tval i2;
                    faceht.GetData (i, i3, i2);
                    if (i2.index != PointIndex::BASE-1)
                      {
                        Element2d tri ( (i2.p4 == PointIndex::BASE-1) ? TRIG : QUAD);
                        for (int l = 0; l < 3; l++)
                          tri[l] = i3.I(l+1);
                        tri.PNum(4) = i2.p4;
                        tri.SetIndex (i2.index);
                        thread_openelements[ti.task_nr].Append (tri);
                      }
                  }
            }}, numtasks);

    for (auto & a : thread_openelements)
      for (auto & el : a)
        openelements.Append (el);
    
    t_pointloop.Stop();
    
    int cnt3 = 0;
    for (int i = 0; i < openelements.Size(); i++)
      if (openelements[i].GetNP() == 3)
        cnt3++;

    int cnt4 = openelements.Size() - cnt3;


    MyStr treequad;
    if (cnt4)
      treequad = MyStr(" (") + MyStr(cnt3) + MyStr (" + ") + 
        MyStr(cnt4) + MyStr(")");

    PrintMessage (5, openelements.Size(), treequad, " open elements");

    BuildBoundaryEdges();


    for (int i = 1; i <= openelements.Size(); i++)
      {
        const Element2d & sel = openelements.Get(i);

        if (boundaryedges)
          for (int j = 1; j <= sel.GetNP(); j++)
            {
              INDEX_2 i2;
              i2.I1() = sel.PNumMod(j);
              i2.I2() = sel.PNumMod(j+1);
              i2.Sort();
              boundaryedges->Set (i2, 1);
            }

        for (int j = 1; j <= 3; j++)
          {
            PointIndex pi = sel.PNum(j);
            // if (pi < points.End())
            if (pi < *points.Range().end())
              points[pi].SetType (FIXEDPOINT);
          }
      }



    /*
      for (i = 1; i <= GetNSeg(); i++)
      {
      const Segment & seg = LineSegment(i);
      INDEX_2 i2(seg[0], seg[1]);
      i2.Sort();

      if (!boundaryedges->Used (i2))
      cerr << "WARNING: no boundedge, but seg edge: " << i2 << endl;

      boundaryedges -> Set (i2, 2);
      segmentht -> Set (i2, i-1);
      }
    */
  }

  bool Mesh :: HasOpenQuads () const
  {
    int no = GetNOpenElements();
    for (int i = 0; i < no; i++)
      if (openelements[i].GetNP() == 4)
        return true;
    return false;
  }





  void Mesh :: FindOpenSegments (int surfnr)
  {
    // int i, j, k;

    // new version, general elements
    // hash index: pnum1-2, surfnr
    // hash data : surfel-nr (pos) or segment nr(neg)
    INDEX_3_HASHTABLE<int> faceht(4 * GetNSE()+GetNSeg()+1);   

    PrintMessage (5, "Test Opensegments");
    for (int i = 1; i <= GetNSeg(); i++)
      {
        const Segment & seg = LineSegment (i);

        if (surfnr == 0 || seg.si == surfnr)
          {
            INDEX_3 key(seg[0], seg[1], seg.si);
            int data = -i;

            if (faceht.Used (key))
              {
                cerr << "ERROR: Segment " << seg << " already used" << endl;
                (*testout) << "ERROR: Segment " << seg << " already used" << endl;
              }

            faceht.Set (key, data);
          }
      }


    /*
      // not possible with surfnr as hash-index
    for (int i = 1; i <= GetNSeg(); i++)
      {
        const Segment & seg = LineSegment (i);

        if (surfnr == 0 || seg.si == surfnr)
          {
            INDEX_2 key(seg[1], seg[0]);
            if (!faceht.Used(key))
              {
                cerr << "ERROR: Segment " << seg << " brother not used" << endl;
                (*testout) << "ERROR: Segment " << seg << " brother not used" << endl;
              }
          }
      }
    */

    
    // bool buggy = false;
    // ofstream bout("buggy.out");

    for (int i = 1; i <= GetNSE(); i++)
      {
        const Element2d & el = SurfaceElement(i);
        if (el.IsDeleted()) continue;

        if (surfnr == 0 || el.GetIndex() == surfnr)
          {
            for (int j = 1; j <= el.GetNP(); j++)
              {
                INDEX_3 seg (el.PNumMod(j), el.PNumMod(j+1), el.GetIndex());
                int data;

                if (seg.I1() < PointIndex::BASE || seg.I2() < PointIndex::BASE)
                  cerr << "seg = " << seg << endl;

                if (faceht.Used(seg))
                  {
                    faceht.Set (seg, 0);
                    /*
                    data = faceht.Get(seg);
                    
                    if (data.I1() == el.GetIndex())
                      {
                        data.I1() = 0;
                        faceht.Set (seg, data);
                      }
                    else
                      {
			// buggy = true;
                        PrintWarning ("hash table si not fitting for segment: ",
                                       seg.I1(), "-", seg.I2(), " other = ",
                                      data.I2(), ", surfnr = ", surfnr);
                      }
                    */
                  }
                else
                  {
                    Swap (seg.I1(), seg.I2());
                    // data.I1() = el.GetIndex();
                    // data.I2() = i;
                    faceht.Set (seg, i);
                  }
              }
          }
      }  

    /*
    if (buggy)
      {
	for (int i = 1; i <= GetNSeg(); i++)
	  bout << "seg" << i << " " << LineSegment(i) << endl;

	for (int i = 1; i <= GetNSE(); i++)
	  bout << "sel" << i << " " << SurfaceElement(i) << " ind = " 
	       << SurfaceElement(i).GetIndex() << endl;

	bout << "hashtable: " << endl;
	for (int j = 1; j <= faceht.GetNBags(); j++)
	  {
	    bout << "bag " << j << ":" << endl;
	    for (int k = 1; k <= faceht.GetBagSize(j); k++)
	      {
		INDEX_2 i2, data;
		faceht.GetData (j, k, i2, data);
		bout << "key = " << i2 << ", data = " << data << endl;
	      }
	  }
	exit(1);
      }
    */

    (*testout) << "open segments: " << endl;
    opensegments.SetSize(0);
    for (int i = 1; i <= faceht.GetNBags(); i++)
      for (int j = 1; j <= faceht.GetBagSize(i); j++)
        {
          INDEX_3 i2;
          int data;
          faceht.GetData (i, j, i2, data);
          if (data)  // surfnr
            {
              Segment seg;
              seg[0] = i2.I1();
              seg[1] = i2.I2();
              seg.si = i2.I3();

              // find geomdata:
              if (data > 0)
                {
                  // segment due to triangle
                  const Element2d & el = SurfaceElement (data);
                  for (int k = 1; k <= el.GetNP(); k++)
                    {
                      if (seg[0] == el.PNum(k))
                        seg.geominfo[0] = el.GeomInfoPi(k);
                      if (seg[1] == el.PNum(k))
                        seg.geominfo[1] = el.GeomInfoPi(k);
                    }

                  (*testout) << "trig seg: ";
                }
              else
                {
                  // segment due to line
                  const Segment & lseg = LineSegment (-data);
                  seg.geominfo[0] = lseg.geominfo[0];
                  seg.geominfo[1] = lseg.geominfo[1];

                  (*testout) << "line seg: ";
                }

              (*testout) << seg[0] << " - " << seg[1] 
                         << " len = " << Dist (Point(seg[0]), Point(seg[1]))
                         << endl;

              opensegments.Append (seg);
              if (seg.geominfo[0].trignum <= 0 || seg.geominfo[1].trignum <= 0)
                {
                  (*testout) << "Problem with open segment: " << seg << endl;
                }

            }
        }

    PrintMessage (3, opensegments.Size(), " open segments found");
    (*testout) << opensegments.Size() << " open segments found" << endl;

    /*
      ptyps.SetSize (GetNP());
      for (i = 1; i <= ptyps.Size(); i++)
      ptyps.Elem(i) = SURFACEPOINT;

      for (i = 1; i <= GetNSeg(); i++)
      {
      const Segment & seg = LineSegment (i);
      ptyps.Elem(seg[0]) = EDGEPOINT;
      ptyps.Elem(seg[1]) = EDGEPOINT;
      }
      for (i = 1; i <= GetNOpenSegments(); i++)
      {
      const Segment & seg = GetOpenSegment (i);
      ptyps.Elem(seg[0]) = EDGEPOINT;
      ptyps.Elem(seg[1]) = EDGEPOINT;
      }
    */
    /*
    for (int i = 1; i <= points.Size(); i++)
      points.Elem(i).SetType(SURFACEPOINT);
    */
    for (auto & p : points)
      p.SetType (SURFACEPOINT);
    
    for (int i = 1; i <= GetNSeg(); i++)
      {
        const Segment & seg = LineSegment (i);
        points[seg[0]].SetType(EDGEPOINT);
        points[seg[1]].SetType(EDGEPOINT);
      }
    for (int i = 1; i <= GetNOpenSegments(); i++)
      {
        const Segment & seg = GetOpenSegment (i);
        points[seg[0]].SetType (EDGEPOINT);
        points[seg[1]].SetType (EDGEPOINT);
      }



    /*

    for (i = 1; i <= openelements.Size(); i++)
    {
    const Element2d & sel = openelements.Get(i);

    if (boundaryedges)
    for (j = 1; j <= sel.GetNP(); j++)
    {
    INDEX_2 i2;
    i2.I1() = sel.PNumMod(j);
    i2.I2() = sel.PNumMod(j+1);
    i2.Sort();
    boundaryedges->Set (i2, 1);
    }

    for (j = 1; j <= 3; j++)
    {
    int pi = sel.PNum(j);
    if (pi <= ptyps.Size())
    ptyps.Elem(pi) = FIXEDPOINT;
    }
    }
    */
  }


  void Mesh :: RemoveOneLayerSurfaceElements ()
  {
    int np = GetNP();

    FindOpenSegments();
    NgBitArray frontpoints(np+1);  // for 0- and 1-based
    frontpoints.Clear();
    
    for (int i = 1; i <= GetNOpenSegments(); i++)
      {
        const Segment & seg = GetOpenSegment(i);
        frontpoints.Set (seg[0]);
        frontpoints.Set (seg[1]);
      }

    for (int i = 1; i <= GetNSE(); i++)
      {
        Element2d & sel = surfelements[i-1];
        bool remove = false;
        for (int j = 1; j <= sel.GetNP(); j++)
          if (frontpoints.Test(sel.PNum(j)))
            remove = true;
        if (remove)
          sel.PNum(1).Invalidate();
      }

    for (int i = surfelements.Size(); i >= 1; i--)
      {
        if (!surfelements[i-1].PNum(1).IsValid())
          {
            surfelements[i-1] = surfelements.Last();
            surfelements.DeleteLast();
          }
      }

    RebuildSurfaceElementLists ();
    /*
    for (int i = 0; i < facedecoding.Size(); i++)
      facedecoding[i].firstelement = -1;
    for (int i = surfelements.Size()-1; i >= 0; i--)
      {
        int ind = surfelements[i].GetIndex();
        surfelements[i].next = facedecoding[ind-1].firstelement;
        facedecoding[ind-1].firstelement = i;
      }
    */

    timestamp = NextTimeStamp();
    //  Compress();
  }





  void Mesh :: FreeOpenElementsEnvironment (int layers)
  {
    int i, j, k;
    PointIndex pi;
    const int large = 9999;
    NgArray<int,PointIndex::BASE> dist(GetNP());

    dist = large;

    for (int i = 1; i <= GetNOpenElements(); i++)
      {
        const Element2d & face = OpenElement(i);
        for (j = 0; j < face.GetNP(); j++)
          dist[face[j]] = 1;
      }

    for (k = 1; k <= layers; k++)
      for (i = 1; i <= GetNE(); i++)
        {
          const Element & el = VolumeElement(i);
          if (el[0] == -1 || el.IsDeleted()) continue;

          int elmin = large;
          for (j = 0; j < el.GetNP(); j++)
            if (dist[el[j]] < elmin)
              elmin = dist[el[j]];

          if (elmin < large)
            {
              for (j = 0; j < el.GetNP(); j++)
                if (dist[el[j]] > elmin+1)
                  dist[el[j]] = elmin+1;
            }
        }

    int cntfree = 0;
    for (i = 1; i <= GetNE(); i++)
      {
        Element & el = VolumeElement(i);
        if (el[0] == -1 || el.IsDeleted()) continue;

        int elmin = large;
        for (j = 0; j < el.GetNP(); j++)
          if (dist[el[j]] < elmin)
            elmin = dist[el[j]];

        el.flags.fixed = elmin > layers;
        // eltyps.Elem(i) = (elmin <= layers) ? 
        // FREEELEMENT : FIXEDELEMENT;
        if (elmin <= layers)
          cntfree++;
      }

    PrintMessage (5, "free: ", cntfree, ", fixed: ", GetNE()-cntfree);
    (*testout) << "free: " << cntfree << ", fixed: " << GetNE()-cntfree << endl;

    for (pi = PointIndex::BASE; 
         pi < GetNP()+PointIndex::BASE; pi++)
      {
        if (dist[pi] > layers+1)
          points[pi].SetType(FIXEDPOINT);
      }
  }



  void Mesh :: SetLocalH (netgen::Point<3> pmin, netgen::Point<3> pmax, double grading)
  {
    using netgen::Point;
    Point<3> c = Center (pmin, pmax);
    double d = max3 (pmax(0)-pmin(0),
                     pmax(1)-pmin(1),
                     pmax(2)-pmin(2));
    d /= 2;
    Point<3> pmin2 = c - Vec<3> (d, d, d);
    Point<3> pmax2 = c + Vec<3> (d, d, d);

    lochfunc = make_unique<LocalH> (pmin2, pmax2, grading, dimension);
  }

  void Mesh :: RestrictLocalH (const Point3d & p, double hloc)
  {
    if(hloc < hmin)
      hloc = hmin;

    //cout << "restrict h in " << p << " to " << hloc << endl;
    if (!lochfunc)
      {
        PrintWarning("RestrictLocalH called, creating mesh-size tree");

        Point3d boxmin, boxmax;
        GetBox (boxmin, boxmax);
        SetLocalH (boxmin, boxmax, 0.8);
      }

    lochfunc -> SetH (p, hloc);
  }

  void Mesh :: RestrictLocalHLine (const Point3d & p1, 
                                   const Point3d & p2,
                                   double hloc)
  {
    if(hloc < hmin)
      hloc = hmin;

    // cout << "restrict h along " << p1 << " - " << p2 << " to " << hloc << endl;
    int i;
    int steps = int (Dist (p1, p2) / hloc) + 2;
    Vec3d v(p1, p2);

    for (i = 0; i <= steps; i++)
      {
        Point3d p = p1 + (double(i)/double(steps) * v);
        RestrictLocalH (p, hloc);
      }
  }


  void Mesh :: SetMinimalH (double h)
  {
    hmin = h;
  }


  void Mesh :: SetGlobalH (double h)
  {
    hglob = h;
  }

  double Mesh :: MaxHDomain (int dom) const
  {
    if (maxhdomain.Size())
      return maxhdomain.Get(dom);
    else
      return 1e10;
  }

  void Mesh :: SetMaxHDomain (const NgArray<double> & mhd)
  {
    maxhdomain.SetSize(mhd.Size());
    for (int i = 1; i <= mhd.Size(); i++)
      maxhdomain.Elem(i) = mhd.Get(i);
  }


  double Mesh :: GetH (const Point3d & p) const
  {
    double hmin = hglob;
    if (lochfunc)
      {
        double hl = lochfunc->GetH (p);
        if (hl < hglob)
          hmin = hl;
      }
    return hmin;
  }

  double Mesh :: GetMinH (const Point3d & pmin, const Point3d & pmax)
  {
    double hmin = hglob;
    if (lochfunc)
      {
        double hl = lochfunc->GetMinH (pmin, pmax);
        if (hl < hmin)
          hmin = hl;
      }
    return hmin;
  }





  double Mesh :: AverageH (int surfnr) const
  {
    int i, j, n;
    double hi, hsum;
    double maxh = 0, minh = 1e10;

    hsum = 0;
    n = 0;
    for (i = 1; i <= GetNSE(); i++)
      {
        const Element2d & el = SurfaceElement(i);
        if (surfnr == 0 || el.GetIndex() == surfnr)
          {
            for (j = 1; j <= 3; j++)
              {
                hi = Dist (Point (el.PNumMod(j)), 
                           Point (el.PNumMod(j+1)));

                hsum += hi;

                if (hi > maxh) maxh = hi;
                if (hi < minh) minh = hi;
                n++;
              }
          }
      }

    PrintMessage (5, "minh = ", minh, " avh = ", (hsum/n), " maxh = ", maxh);
    return (hsum / n);
  }



  void Mesh :: CalcLocalH (double grading) 
  {
    static Timer t("Mesh::CalcLocalH"); RegionTimer reg(t);
    
    if (!lochfunc)
      {
        Point3d pmin, pmax;
        GetBox (pmin, pmax);
        // SetLocalH (pmin, pmax, mparam.grading);
	SetLocalH (pmin, pmax, grading);
      }

    PrintMessage (3,
                  "CalcLocalH: ", 
                  GetNP(), " Points ", 
                  GetNE(), " Elements ", 
                  GetNSE(), " Surface Elements");


    for (int i = 0; i < GetNSE(); i++)
      {
        const Element2d & el = surfelements[i];
        int j;

        if (el.GetNP() == 3)
          {
            double hel = -1;
            for (j = 1; j <= 3; j++)
              {
                const Point3d & p1 = points[el.PNumMod(j)];
                const Point3d & p2 = points[el.PNumMod(j+1)];

                /*
                  INDEX_2 i21(el.PNumMod(j), el.PNumMod(j+1));
                  INDEX_2 i22(el.PNumMod(j+1), el.PNumMod(j));
                  if (! identifiedpoints->Used (i21) &&
                  ! identifiedpoints->Used (i22) )
                */
                if (!ident -> UsedSymmetric (el.PNumMod(j),
                                             el.PNumMod(j+1)))
                  {
                    double hedge = Dist (p1, p2);
                    if (hedge > hel)
                      hel = hedge;
                    //		  lochfunc->SetH (Center (p1, p2), 2 * Dist (p1, p2));
                    //		  (*testout) << "trigseth, p1,2 = " << el.PNumMod(j) << ", " << el.PNumMod(j+1) 
                    //			     << " h = " << (2 * Dist(p1, p2)) << endl;
                  }
              }

            if (hel > 0)
              {
                const Point3d & p1 = points[el.PNum(1)];
                const Point3d & p2 = points[el.PNum(2)];
                const Point3d & p3 = points[el.PNum(3)];
                lochfunc->SetH (Center (p1, p2, p3), hel);
              }
          }
        else
          {
            {
              const Point3d & p1 = points[el.PNum(1)];
              const Point3d & p2 = points[el.PNum(2)];
              lochfunc->SetH (Center (p1, p2), 2 * Dist (p1, p2));
            }
            {
              const Point3d & p1 = points[el.PNum(3)];
              const Point3d & p2 = points[el.PNum(4)];
              lochfunc->SetH (Center (p1, p2), 2 * Dist (p1, p2));
            }
          }
      }

    for (int i = 0; i < GetNSeg(); i++)
      {
        const Segment & seg = segments[i];
        const Point3d & p1 = points[seg[0]];
        const Point3d & p2 = points[seg[1]];
        /*
          INDEX_2 i21(seg[0], seg[1]);
          INDEX_2 i22(seg[1], seg[0]);
          if (identifiedpoints)
          if (!identifiedpoints->Used (i21) && !identifiedpoints->Used (i22))
        */
        if (!ident -> UsedSymmetric (seg[0], seg[1]))
          {
            lochfunc->SetH (Center (p1, p2), Dist (p1, p2));
          }
      }
    /*
      cerr << "do vol" << endl;
      for (i = 1; i <= GetNE(); i++)
      {
      const Element & el = VolumeElement(i);
      if (el.GetType() == TET)
      {
      int j, k;
      for (j = 2; j <= 4; j++)
      for (k = 1; k < j; k++)  
      {
      const Point3d & p1 = Point (el.PNum(j));
      const Point3d & p2 = Point (el.PNum(k));
      lochfunc->SetH (Center (p1, p2), 2 * Dist (p1, p2));
      (*testout) << "set vol h to " << (2 * Dist (p1, p2)) << endl;

      }
      }
      }
    */

    /*
      const char * meshsizefilename = 
      globflags.GetStringFlag ("meshsize", NULL);
      if (meshsizefilename)
      {
      ifstream msf(meshsizefilename);
      if (msf)
      {
      int nmsp;
      msf >> nmsp;
      for (i = 1; i <= nmsp; i++)
      {
      Point3d pi;
      double hi;
      msf >> pi.X() >> pi.Y() >> pi.Z();
      msf >> hi;
      lochfunc->SetH (pi, hi);
      }
      }
      }
    */
    //  lochfunc -> Convexify();
    //  lochfunc -> PrintMemInfo (cout);
  }


  void Mesh :: CalcLocalHFromPointDistances(double grading)
  {
    PrintMessage (3, "Calculating local h from point distances");

    if (!lochfunc)
      {
        Point3d pmin, pmax;
        GetBox (pmin, pmax);

        // SetLocalH (pmin, pmax, mparam.grading);
	SetLocalH (pmin, pmax, grading);
      }

    PointIndex i,j;
    double hl;


    for (i = PointIndex::BASE; 
         i < GetNP()+PointIndex::BASE; i++)
      {
        for(j=i+1; j<GetNP()+PointIndex::BASE; j++)
          {
            const Point3d & p1 = points[i];
            const Point3d & p2 = points[j];
            hl = Dist(p1,p2);
            RestrictLocalH(p1,hl);
            RestrictLocalH(p2,hl);
            //cout << "restricted h at " << p1 << " and " << p2 << " to " << hl << endl;
          }
      }


  }


  void Mesh :: CalcLocalHFromSurfaceCurvature (double grading, double elperr) 
  {
    PrintMessage (3, "Calculating local h from surface curvature");

    if (!lochfunc)
      {
        Point3d pmin, pmax;
        GetBox (pmin, pmax);

        // SetLocalH (pmin, pmax, mparam.grading);
	SetLocalH (pmin, pmax, grading);
      }


    INDEX_2_HASHTABLE<int> edges(3 * GetNP() + 2);
    INDEX_2_HASHTABLE<int> bedges(GetNSeg() + 2);
    int i, j;

    for (i = 1; i <= GetNSeg(); i++)
      {
        const Segment & seg = LineSegment(i);
        INDEX_2 i2(seg[0], seg[1]);
        i2.Sort();
        bedges.Set (i2, 1);
      }
    for (i = 1; i <= GetNSE(); i++)
      {
        const Element2d & sel = SurfaceElement(i);
        if (!sel.PNum(1))
          continue;
        for (j = 1; j <= 3; j++)
          {
            INDEX_2 i2(sel.PNumMod(j), sel.PNumMod(j+1));
            i2.Sort();
            if (bedges.Used(i2)) continue;

            if (edges.Used(i2))
              {
                int other = edges.Get(i2);

                const Element2d & elother = SurfaceElement(other);

                int pi3 = 1;
                while ( (sel.PNum(pi3) == i2.I1()) || 
                        (sel.PNum(pi3) == i2.I2()))
                  pi3++;
                pi3 = sel.PNum(pi3);

                int pi4 = 1;
                while ( (elother.PNum(pi4) == i2.I1()) || 
                        (elother.PNum(pi4) == i2.I2()))
                  pi4++;
                pi4 = elother.PNum(pi4);

                double rad = ComputeCylinderRadius (Point (PointIndex(i2.I1())),
                                                    Point (PointIndex(i2.I2())),
                                                    Point (PointIndex(pi3)),
                                                    Point (PointIndex(pi4)));

                RestrictLocalHLine (Point(PointIndex(i2.I1())), Point(PointIndex(i2.I2())), rad/elperr);


                /*	      
                  (*testout) << "pi1,2, 3, 4 = " << i2.I1() << ", " << i2.I2() << ", " << pi3 << ", " << pi4
                  << " p1 = " << Point(i2.I1()) 
                  << ", p2 = " << Point(i2.I2()) 
                  //			 << ", p3 = " << Point(pi3) 
                  //			 << ", p4 = " << Point(pi4) 
                  << ", rad = " << rad << endl;
                */
              }
            else
              edges.Set (i2, i);
          }
      }


    // Restrict h due to line segments

    for (i = 1; i <= GetNSeg(); i++)
      {
        const Segment & seg = LineSegment(i);
        const Point3d & p1 = Point(seg[0]);
        const Point3d & p2 = Point(seg[1]);
        RestrictLocalH (Center (p1, p2),  Dist (p1, p2));
      }



    /*


    int i, j;
    int np = GetNP();
    int nseg = GetNSeg();
    int nse = GetNSE();

    NgArray<Vec3d> normals(np);
    NgBitArray linepoint(np);

    linepoint.Clear();
    for (i = 1; i <= nseg; i++)
    {
    linepoint.Set (LineSegment(i)[0]);
    linepoint.Set (LineSegment(i)[1]);
    }

    for (i = 1; i <= np; i++)
    normals.Elem(i) = Vec3d(0,0,0);

    for (i = 1; i <= nse; i++)
    {
    Element2d & el = SurfaceElement(i);
    Vec3d nf = Cross (Vec3d (Point (el.PNum(1)), Point(el.PNum(2))),
    Vec3d (Point (el.PNum(1)), Point(el.PNum(3))));
    for (j = 1; j <= 3; j++)
    normals.Elem(el.PNum(j)) += nf;
    }

    for (i = 1; i <= np; i++)
    normals.Elem(i) /= (1e-12 + normals.Elem(i).Length());

    for (i = 1; i <= nse; i++)
    {
    Element2d & el = SurfaceElement(i);
    Vec3d nf = Cross (Vec3d (Point (el.PNum(1)), Point(el.PNum(2))),
    Vec3d (Point (el.PNum(1)), Point(el.PNum(3))));
    nf /= nf.Length();
    Point3d c = Center (Point(el.PNum(1)),
    Point(el.PNum(2)),
    Point(el.PNum(3)));

    for (j = 1; j <= 3; j++)
    {
    if (!linepoint.Test (el.PNum(j)))
    {
    double dist = Dist (c, Point(el.PNum(j)));
    double dn = (nf - normals.Get(el.PNum(j))).Length();

    RestrictLocalH (Point(el.PNum(j)), dist / (dn+1e-12) /elperr);
    }
    }
    }
    */
  }


  void Mesh :: RestrictLocalH (resthtype rht, int nr, double loch)
  {
    int i;
    switch (rht)
      {
      case RESTRICTH_FACE:
        {
          for (i = 1; i <= GetNSE(); i++)
            {
              const Element2d & sel = SurfaceElement(i);
              if (sel.GetIndex() == nr)
                RestrictLocalH (RESTRICTH_SURFACEELEMENT, i, loch);
            }
          break;
        }
      case RESTRICTH_EDGE:
        {
          for (i = 1; i <= GetNSeg(); i++)
            {
              const Segment & seg = LineSegment(i);
              if (seg.edgenr == nr)
                RestrictLocalH (RESTRICTH_SEGMENT, i, loch);
            }
          break;
        }
      case RESTRICTH_POINT:
        {
          RestrictLocalH (Point (nr), loch);
          break;
        }

      case RESTRICTH_SURFACEELEMENT:
        {
          const Element2d & sel = SurfaceElement(nr);
          Point3d p = Center (Point(sel.PNum(1)),
                              Point(sel.PNum(2)),
                              Point(sel.PNum(3)));
          RestrictLocalH (p, loch);
          break;
        }
      case RESTRICTH_SEGMENT:
        {
          const Segment & seg = LineSegment(nr);
          RestrictLocalHLine (Point (seg[0]), Point(seg[1]), loch);
          break;
        }
      }
  }


  void Mesh :: LoadLocalMeshSize (const string &  meshsizefilename)
  {
    // Philippose - 10/03/2009
    // Improve error checking when loading and reading
    // the local mesh size file

    if (meshsizefilename.empty()) return;

    ifstream msf(meshsizefilename.c_str());

    // Philippose - 09/03/2009
    // Adding print message information in case the specified 
    // does not exist, or does not load successfully due to 
    // other reasons such as access rights, etc...
    if (!msf) 
      {
        PrintMessage(3, "Error loading mesh size file: ", meshsizefilename, "....","Skipping!");
        return;
      }

    PrintMessage (3, "Load local mesh-size file: ", meshsizefilename);

    int nmsp = 0;
    int nmsl = 0;

    msf >> nmsp;
    if(!msf.good())
      throw NgException ("Mesh-size file error: No points found\n");

    if(nmsp > 0)
      PrintMessage (4, "Number of mesh-size restriction points: ", nmsp);

    for (int i = 0; i < nmsp; i++)
      {
        Point3d pi;
        double hi;
        msf >> pi.X() >> pi.Y() >> pi.Z();
        msf >> hi;
        if (!msf.good())
          throw NgException ("Mesh-size file error: Number of points don't match specified list size\n");
        RestrictLocalH (pi, hi);
      }

    msf >> nmsl;
    if(!msf.good())
      throw NgException ("Mesh-size file error: No line definitions found\n");

    if(nmsl > 0)
      PrintMessage (4, "Number of mesh-size restriction lines: ", nmsl);

    for (int i = 0; i < nmsl; i++)
      {
        Point3d p1, p2;
        double hi;
        msf >> p1.X() >> p1.Y() >> p1.Z();
        msf >> p2.X() >> p2.Y() >> p2.Z();
        msf >> hi;
        if (!msf.good())
          throw NgException ("Mesh-size file error: Number of line definitions don't match specified list size\n");
        RestrictLocalHLine (p1, p2, hi);
      }

    msf.close();
  }



  void Mesh :: GetBox (Point3d & pmin, Point3d & pmax, int dom) const
  {
    if (points.Size() == 0)
      {
        pmin = pmax = Point3d(0,0,0);
        return;
      }

    if (dom <= 0)
      {
        pmin = Point3d (1e10, 1e10, 1e10);
        pmax = Point3d (-1e10, -1e10, -1e10); 

        // for (PointIndex pi = points.Begin(); pi < points.End(); pi++)
        for (PointIndex pi : points.Range())
          {
            pmin.SetToMin ( (*this) [pi] );
            pmax.SetToMax ( (*this) [pi] );
          }
      }
    else
      {
        int j, nse = GetNSE();
        SurfaceElementIndex sei;

        pmin = Point3d (1e10, 1e10, 1e10);
        pmax = Point3d (-1e10, -1e10, -1e10); 
        for (sei = 0; sei < nse; sei++)
          {
            const Element2d & el = (*this)[sei];
            if (el.IsDeleted() ) continue;

            if (dom == -1 || el.GetIndex() == dom)
              {
                for (j = 0; j < 3; j++)
                  {
                    pmin.SetToMin ( (*this) [el[j]] );
                    pmax.SetToMax ( (*this) [el[j]] );
                  }
              }
          }
      }

    if (pmin.X() > 0.5e10)
      {
        pmin = pmax = Point3d(0,0,0);
      }
  }




  void Mesh :: GetBox (Point3d & pmin, Point3d & pmax, POINTTYPE ptyp) const
  {
    if (points.Size() == 0)
      {
        pmin = pmax = Point3d(0,0,0);
        return;
      }

    pmin = Point3d (1e10, 1e10, 1e10);
    pmax = Point3d (-1e10, -1e10, -1e10); 

    // for (PointIndex pi = points.Begin(); pi < points.End(); pi++)
    for (PointIndex pi : points.Range())
      if (points[pi].Type() <= ptyp)
        {
          pmin.SetToMin ( (*this) [pi] );
          pmax.SetToMax ( (*this) [pi] );
        }
  }




  double Mesh :: ElementError (int eli, const MeshingParameters & mp) const
  {
    const Element & el = volelements[eli-1];
    return CalcTetBadness (points[el[0]], points[el[1]],
                           points[el[2]], points[el[3]], -1, mp);
  }

  void Mesh :: AddLockedPoint (PointIndex pi)
  { 
    lockedpoints.Append (pi); 
  }

  void Mesh :: ClearLockedPoints ()
  { 
    lockedpoints.SetSize (0); 
  }



  void Mesh :: Compress ()
  {
    static Timer t("Mesh::Compress"); RegionTimer reg(t);
    NgLock lock(mutex);
    lock.Lock();
    
    Array<PointIndex,PointIndex> op2np(GetNP());
    Array<bool, PointIndex> pused(GetNP());

    /*
      (*testout) << "volels: " << endl;
      for (i = 1; i <= volelements.Size(); i++)
      {
      for (j = 1; j <= volelements.Get(i).GetNP(); j++)
      (*testout) << volelements.Get(i).PNum(j) << " ";
      (*testout) << endl;
      }
      (*testout) << "np: " << GetNP() << endl;
    */

    for (int i = 0; i < volelements.Size(); i++)
      if (volelements[i][0] <= PointIndex::BASE-1 ||
          volelements[i].IsDeleted())
        {
          volelements.DeleteElement(i);
          i--;
        }


    for (int i = 0; i < surfelements.Size(); i++)
      if (surfelements[i].IsDeleted())
        {
          surfelements.DeleteElement(i);
          i--;
        }

    for (int i = 0; i < segments.Size(); i++)
      if (segments[i][0] <= PointIndex::BASE-1)
        {
          segments.DeleteElement(i);
          i--;
        }

    for(int i=0; i < segments.Size(); i++)
      if(segments[i].edgenr < 0)
          segments.DeleteElement(i--);

    pused = false;
    /*
    for (int i = 0; i < volelements.Size(); i++)
      {
        const Element & el = volelements[i];
        for (int j = 0; j < el.GetNP(); j++)
          pused[el[j]] = true;
      }
    */
    /*
    for (const Element & el : volelements)
      for (PointIndex pi : el.PNums())
        pused[pi] = true;
    */

    ParallelForRange
      (volelements.Range(), [&] (auto myrange)
       {
         for (const Element & el : volelements.Range(myrange))
           for (PointIndex pi : el.PNums())
             pused[pi] = true;
       });

    /*
    for (int i = 0; i < surfelements.Size(); i++)
      {
        const Element2d & el = surfelements[i];
        for (int j = 0; j < el.GetNP(); j++)
          pused[el[j]] = true;
      }
    */
    ParallelForRange
      (surfelements.Range(), [&] (auto myrange)
       {
         for (const Element2d & el : surfelements.Range(myrange))
           for (PointIndex pi : el.PNums())
             pused[pi] = true;
       });
    
    for (int i = 0; i < segments.Size(); i++)
      {
        const Segment & seg = segments[i];
        for (int j = 0; j < seg.GetNP(); j++)
          pused[seg[j]] = true;
      }

    for (int i = 0; i < openelements.Size(); i++)
      {
        const Element2d & el = openelements[i];
        for (int j = 0; j < el.GetNP(); j++)
          pused[el[j]] = true;
      }

    for (int i = 0; i < lockedpoints.Size(); i++)
      pused[lockedpoints[i]] = true;


    /*
    // compress points doesn't work for identified points !
    if (identifiedpoints)
    {
    for (i = 1; i <= identifiedpoints->GetNBags(); i++)
    if (identifiedpoints->GetBagSize(i))
    {
    pused.Set ();
    break;
    }
    }
    */
    //  pused.Set();

    
    {
      Array<MeshPoint> hpoints;
      int npi = PointIndex::BASE;
      for (PointIndex pi : points.Range())
        if (pused[pi])
          {
            op2np[pi] = npi;
            npi++;
            hpoints.Append (points[pi]);
          }
        else
          {
            op2np[pi].Invalidate(); 
          }
      
      points.SetSize(0);
      for (int i = 0; i < hpoints.Size(); i++)
        points.Append (hpoints[i]);
    }
    
    /*
    for (int i = 1; i <= volelements.Size(); i++)
      {
        Element & el = VolumeElement(i);
        for (int j = 0; j < el.GetNP(); j++)
          el[j] = op2np[el[j]];
      }
    */
    ParallelForRange
      (volelements.Range(), [&] (auto myrange)
       {
         for (Element & el : volelements.Range(myrange))
           for (PointIndex & pi : el.PNums())
             pi = op2np[pi];
       });

    /*
    for (int i = 1; i <= surfelements.Size(); i++)
      {
        Element2d & el = SurfaceElement(i);
        for (int j = 0; j < el.GetNP(); j++)
          el[j] = op2np[el[j]];
      }
    */
    ParallelForRange
      (surfelements.Range(), [&] (auto myrange)
       {
         for (Element2d & el : surfelements.Range(myrange))
           for (PointIndex & pi : el.PNums())
             pi = op2np[pi];
       });

    
    for (int i = 0; i < segments.Size(); i++)
      {
        Segment & seg = segments[i];
        for (int j = 0; j < seg.GetNP(); j++)
          seg[j] = op2np[seg[j]];
      }

    for (int i = 1; i <= openelements.Size(); i++)
      {
        Element2d & el = openelements.Elem(i);
        for (int j = 0; j < el.GetNP(); j++)
          el[j] = op2np[el[j]];
      }  


    for (int i = 0; i < lockedpoints.Size(); i++)
      lockedpoints[i] = op2np[lockedpoints[i]];
    /*
    for (int i = 0; i < facedecoding.Size(); i++)
      facedecoding[i].firstelement = -1;
    for (int i = surfelements.Size()-1; i >= 0; i--)
      {
        int ind = surfelements[i].GetIndex();
        surfelements[i].next = facedecoding[ind-1].firstelement;
        facedecoding[ind-1].firstelement = i;
      }
    */
    RebuildSurfaceElementLists ();
    CalcSurfacesOfNode();


    //  FindOpenElements();
    timestamp = NextTimeStamp();
    lock.UnLock();
  }

  void Mesh :: OrderElements()
  {
    for (auto & el : surfelements)
      {
        if (el.GetType() == TRIG)
          while (el[0] > el[1] || el[0] > el[2])
            { // rotate element
              auto hp = el[0];
              el[0] = el[1];
              el[1] = el[2];
              el[2] = hp;
              auto hgi = el.GeomInfoPi(1);
              el.GeomInfoPi(1) = el.GeomInfoPi(2);
              el.GeomInfoPi(2) = el.GeomInfoPi(3);
              el.GeomInfoPi(3) = hgi;
            }
      }

    for (auto & el : volelements)
      if (el.GetType() == TET)
        {
          // lowest index first ...
          int mini = 0;
          for (int i = 1; i < 4; i++)
            if (el[i] < el[mini]) mini = i;
          if (mini != 0)
            { // swap 0 with mini, and the other two ...
              int i3 = -1, i4 = -1;
              for (int i = 1; i < 4; i++)
                if (i != mini)
                  {
                    i4 = i3;
                    i3 = i;
                  }
              swap (el[0], el[mini]);
              swap (el[i3], el[i4]);
            }
          
          while (el[1] > el[2] || el[1] > el[3])
            { // rotate element to move second index to second position
              auto hp = el[1];
              el[1] = el[2];
              el[2] = el[3];
              el[3] = hp;
            }
        }
  }

  int Mesh :: CheckConsistentBoundary () const
  {
    int nf = GetNOpenElements();
    INDEX_2_HASHTABLE<int> edges(nf+2);
    INDEX_2 i2, i2s, edge;
    int err = 0;

    for (int i = 1; i <= nf; i++)
      {
        const Element2d & sel = OpenElement(i);

        for (int j = 1; j <= sel.GetNP(); j++)
          {
            i2.I1() = sel.PNumMod(j);
            i2.I2() = sel.PNumMod(j+1);

            int sign = (i2.I2() > i2.I1()) ? 1 : -1;
            i2.Sort();
            if (!edges.Used (i2))
              edges.Set (i2, 0);
            edges.Set (i2, edges.Get(i2) + sign);
          }
      }

    for (int i = 1; i <= edges.GetNBags(); i++)
      for (int j = 1; j <= edges.GetBagSize(i); j++)
        {
          int cnt = 0;
          edges.GetData (i, j, i2, cnt);
          if (cnt)
            {
              PrintError ("Edge ", i2.I1() , " - ", i2.I2(), " multiple times in surface mesh");

              (*testout) << "Edge " << i2 << " multiple times in surface mesh" << endl;
              i2s = i2;
              i2s.Sort();
              for (int k = 1; k <= nf; k++)
                {
                  const Element2d & sel = OpenElement(k);
                  for (int l = 1; l <= sel.GetNP(); l++)
                    {
                      edge.I1() = sel.PNumMod(l);
                      edge.I2() = sel.PNumMod(l+1);
                      edge.Sort();

                      if (edge == i2s) 
                        (*testout) << "edge of element " << sel << endl;
                    }
                }


              err = 2;
            }
        }

    return err;
  }



  int Mesh :: CheckOverlappingBoundary () 
  {
    static Timer t("Mesh::CheckOverlappingBoundary"); RegionTimer reg(t);
    
    Point3d pmin, pmax;
    GetBox (pmin, pmax);
    BoxTree<3, SurfaceElementIndex> setree(pmin, pmax);
    // NgArray<SurfaceElementIndex> inters;

    bool overlap = 0;
    bool incons_layers = 0;

    for (Element2d & el : SurfaceElements())
      el.badel = false;

    for (SurfaceElementIndex sei : Range(SurfaceElements()))
      {
        const Element2d & tri = SurfaceElement(sei);

        Box<3> box(Box<3>::EMPTY_BOX);
        for (PointIndex pi : tri.PNums())
          box.Add (Point(pi));

        box.Increase(1e-3*box.Diam());
        setree.Insert (box, sei);
      }

    std::mutex m;
    // for (SurfaceElementIndex sei : Range(SurfaceElements()))
    ParallelForRange
      (Range(SurfaceElements()), [&] (auto myrange)
       {
         for (SurfaceElementIndex sei : myrange)
           {
             const Element2d & tri = SurfaceElement(sei);
             
             Box<3> box(Box<3>::EMPTY_BOX);
             for (PointIndex pi : tri.PNums())
               box.Add (Point(pi));
             
             setree.GetFirstIntersecting
               (box.PMin(), box.PMax(),
                [&] (SurfaceElementIndex sej) 
                {
                  const Element2d & tri2 = SurfaceElement(sej);	  
                  
                  if ( (*this)[tri[0]].GetLayer() != (*this)[tri2[0]].GetLayer())
                    return false;
                  
                  if ( (*this)[tri[0]].GetLayer() != (*this)[tri[1]].GetLayer() ||
                       (*this)[tri[0]].GetLayer() != (*this)[tri[2]].GetLayer())
                    {
                      incons_layers = 1;
                      // cout << "inconsistent layers in triangle" << endl;
                    }
                  
                  const netgen::Point<3> *trip1[3], *trip2[3];	  
                  for (int k = 0; k < 3; k++)
                    {
                      trip1[k] = &Point (tri[k]);
                      trip2[k] = &Point (tri2[k]);
                    }
                  
                  if (IntersectTriangleTriangle (&trip1[0], &trip2[0]))
                    {
                      overlap = 1;
                      lock_guard<std::mutex> guard(m);
                      if(!incons_layers)
                        {
                          PrintWarning ("Intersecting elements "
                                        ,int(sei), " and ", int(sej));
                      
                          (*testout) << "Intersecting: " << endl;
                          (*testout) << "openelement " << sei << " with open element " << sej << endl;
                      
                          cout << "el1 = " << tri << endl;
                          cout << "el2 = " << tri2 << endl;
                          cout << "layer1 = " <<  (*this)[tri[0]].GetLayer() << endl;
                          cout << "layer2 = " <<  (*this)[tri2[0]].GetLayer() << endl;
                        }
                      
                      for (int k = 1; k <= 3; k++)
                        (*testout) << tri.PNum(k) << "  ";
                      (*testout) << endl;
                      for (int k = 1; k <= 3; k++)
                        (*testout) << tri2.PNum(k) << "  ";
                      (*testout) << endl;
                      
                      for (int k = 0; k <= 2; k++)
                        (*testout) << *trip1[k] << "   ";
                      (*testout) << endl;
                      for (int k = 0; k <= 2; k++)
                        (*testout) << *trip2[k] << "   ";
                      (*testout) << endl;

                      (*testout) << "Face1 = " << GetFaceDescriptor(tri.GetIndex()) << endl;
                      (*testout) << "Face1 = " << GetFaceDescriptor(tri2.GetIndex()) << endl;
                      
                      SurfaceElement(sei).badel = 1;
                      SurfaceElement(sej).badel = 1;
                    }
                  return false;
                });
           }
       });
    // bug 'fix'
    if (incons_layers) overlap = 0;

    return overlap;
  }


  int Mesh :: CheckVolumeMesh () const
  {
    PrintMessage (3, "Checking volume mesh");

    int ne = GetNE();
    DenseMatrix dtrans(3,3);
    int i, j;

    PrintMessage (5, "elements: ", ne);
    for (i = 1; i <= ne; i++)
      {
        Element & el = (Element&) VolumeElement(i);
        el.flags.badel = 0;
        int nip = el.GetNIP();
        for (j = 1; j <= nip; j++)
          {
            el.GetTransformation (j, Points(), dtrans);
            double det = dtrans.Det();
            if (det > 0)
              {
                PrintError ("Element ", i , " has wrong orientation");
                el.flags.badel = 1;
              }
          }
      }

    return 0;
  }

  // Search for surface trigs with same vertices ( may happen for instance with close surfaces in stl geometies )
  int Mesh :: FindIllegalTrigs ()
  {
    // Temporary table to store the vertex numbers of all triangles
    INDEX_3_CLOSED_HASHTABLE<int> temp_tab(3*GetNSE() + 1);
    size_t cnt = 0;
    for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
      {
        const Element2d & sel = surfelements[sei];
        if (sel.IsDeleted()) continue;

        INDEX_3 i3(sel[0], sel[1], sel[2]);
        i3.Sort();
        if(temp_tab.Used(i3))
          {
            temp_tab.Set (i3, -1);
            cnt++;
          }
        else
          {
            temp_tab.Set (i3, sei);
          }
      }

    illegal_trigs = make_unique<INDEX_3_CLOSED_HASHTABLE<int>> (2*cnt+1);
    for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
      {
        const Element2d & sel = surfelements[sei];
        if (sel.IsDeleted()) continue;

        INDEX_3 i3(sel[0], sel[1], sel[2]);
        i3.Sort();
        if(temp_tab.Get(i3)==-1)
            illegal_trigs -> Set (i3, 1);
      }
    return cnt;
  }

  bool Mesh :: LegalTrig (const Element2d & el) const
  {
      if(illegal_trigs)
      {
          INDEX_3 i3 (el[0], el[1], el[2]);
          i3.Sort();
          if(illegal_trigs->Used(i3))
              return false;
      }

    return 1;
    if ( /* hp */ 1)  // needed for old, simple hp-refinement
      { 
        // trigs with 2 or more segments are illegal
        int i;
        int nseg = 0;

        if (!segmentht)
          {
            cerr << "no segmentht allocated" << endl;
            return 0;
          }

        //      Point3d cp(0.5, 0.5, 0.5);
        for (i = 1; i <= 3; i++)
          {
            INDEX_2 i2(el.PNumMod (i), el.PNumMod (i+1));
            i2.Sort();
            if (segmentht -> Used (i2))
              nseg++;
          }
        if (nseg >= 2) 
          return 0;
      }
    return 1;
  }

  double Mesh :: CalcTotalBad (const MeshingParameters & mp )
  {
    static Timer t("CalcTotalBad"); RegionTimer reg(t);
    static constexpr int n_classes = 20;

    double sum = 0;

    tets_in_qualclass.SetSize(n_classes);
    tets_in_qualclass = 0;

    ParallelForRange( IntRange(volelements.Size()), [&] (auto myrange)
       {
         double local_sum = 0.0;
         double teterrpow = mp.opterrpow;

         std::array<int,n_classes> classes_local{};

         for (auto i : myrange)
           {
             double elbad = pow (max2(CalcBad (points, volelements[i], 0, mp),1e-10), 1/teterrpow);

             int qualclass = int (n_classes / elbad + 1);
             if (qualclass < 1) qualclass = 1;
             if (qualclass > n_classes) qualclass = n_classes;
             classes_local[qualclass-1]++;

             local_sum += elbad;
           }

         AtomicAdd(sum, local_sum);

         for (auto i : Range(n_classes))
             AsAtomic(tets_in_qualclass[i]) += classes_local[i];
    });

    return sum;
  }




  ///
  bool Mesh :: LegalTet2 (Element & el) const
  {
    // static int timer1 = NgProfiler::CreateTimer ("Legaltet2");

    // Test, whether 4 points have a common surface plus
    // at least 4 edges at the boundary

    if(!boundaryedges)
      const_cast<Mesh *>(this)->BuildBoundaryEdges();


    // non-tets are always legal
    if (el.GetType() != TET)
      {
        el.SetLegal (1);
        return 1;
      }

    POINTTYPE pointtype[4];
    for(int i = 0; i < 4; i++)
      pointtype[i] = (*this)[el[i]].Type();



    // element has at least 2 inner points ---> legal
    int cnti = 0;
    for (int j = 0; j < 4; j++)
      if ( pointtype[j] == INNERPOINT)
        {
          cnti++;
          if (cnti >= 2)
            {
              el.SetLegal (1);
              return 1;
            }
        }



    // which faces are boundary faces ?
    int bface[4];
    for (int i = 0; i < 4; i++)
      {
        bface[i] = surfelementht->Used (INDEX_3::Sort(el[gftetfacesa[i][0]],
                                                      el[gftetfacesa[i][1]],
                                                      el[gftetfacesa[i][2]]));
      }

    int bedge[4][4];
    int segedge[4][4];
    static const int pi3map[4][4] = { { -1,  2,  1,  1 },
                                      {  2, -1,  0,  0 },
                                      {  1,  0, -1,  0 },
                                      {  1,  0,  0, -1 } };

    static const int pi4map[4][4] = { { -1,  3,  3,  2 },
                                      {  3, -1,  3,  2 },
                                      {  3,  3, -1,  1 },
                                      {  2,  2,  1, -1 } };


    for (int i = 0; i < 4; i++)
      for (int j = 0; j < i; j++)
        {
          bool sege = false, be = false;

          int pos = boundaryedges -> Position0(INDEX_2::Sort(el[i], el[j]));
          if (pos != -1)
            {
              be = true;
              if (boundaryedges -> GetData0(pos) == 2)
                sege = true;
            }

          segedge[j][i] = segedge[i][j] = sege;
          bedge[j][i] = bedge[i][j] = be;
        }

    // two boundary faces and no edge is illegal
    for (int i = 0; i < 3; i++)
      for (int j = i+1; j < 4; j++)
        {
          if (bface[i] && bface[j])
            if (!segedge[pi3map[i][j]][pi4map[i][j]])
              {
                // 2 boundary faces withoud edge in between
                el.SetLegal (0);
                return 0;
              }
        }

    // three boundary edges meeting in a Surface point
    for (int i = 0; i < 4; i++)
      {
        if ( pointtype[i] == SURFACEPOINT)
          {
            bool alledges = 1;
            for (int j = 0; j < 4; j++)
              if (j != i && !bedge[i][j])
                {
                  alledges = 0;
                  break;
                }
            if (alledges)
              {
                // cout << "tet illegal due to unmarked node" << endl;
                el.SetLegal (0);
                return 0;
              }
          }
      }



    for (int fnr = 0; fnr < 4; fnr++)
      if (!bface[fnr])
        for (int i = 0; i < 4; i++)
          if (i != fnr)
            {
              int pi1 = pi3map[i][fnr];
              int pi2 = pi4map[i][fnr];

              if ( pointtype[i] == SURFACEPOINT)
                {
                  // two connected edges on surface, but no face
                  if (bedge[i][pi1] && bedge[i][pi2])
                    {
                      el.SetLegal (0);
                      return 0;
                    }
                }

              if ( pointtype[i] == EDGEPOINT)
                {
                  // connected surface edge and edge edge, but no face
                  if ( (bedge[i][pi1] && segedge[i][pi2]) ||
                       (bedge[i][pi2] && segedge[i][pi1]) )
                    {
                      el.SetLegal (0);
                      return 0;
                    }
                }

            }


    el.SetLegal (1);
    return 1;

  }



  int Mesh :: GetNDomains() const
  {
    int ndom = 0;

    for (int k = 0; k < facedecoding.Size(); k++)
      {
        if (facedecoding[k].DomainIn() > ndom)
          ndom = facedecoding[k].DomainIn();
        if (facedecoding[k].DomainOut() > ndom)
          ndom = facedecoding[k].DomainOut();
      }

    return ndom;
  }



  void Mesh :: SurfaceMeshOrientation ()
  {
    int i, j;
    int nse = GetNSE();

    NgBitArray used(nse);
    used.Clear();
    INDEX_2_HASHTABLE<int> edges(nse+1);

    bool haschanged = 0;


    const Element2d & tri = SurfaceElement(1);
    for (j = 1; j <= 3; j++)
      {
        INDEX_2 i2(tri.PNumMod(j), tri.PNumMod(j+1));
        edges.Set (i2, 1);
      }
    used.Set(1);

    bool unused;
    do
      {
        bool changed;
        do
          {
            changed = 0;
            for (i = 1; i <= nse; i++)
              if (!used.Test(i))
                {
                  Element2d & el = surfelements[i-1];
                  int found = 0, foundrev = 0;
                  for (j = 1; j <= 3; j++)
                    {
                      INDEX_2 i2(el.PNumMod(j), el.PNumMod(j+1));
                      if (edges.Used(i2))
                        foundrev = 1;
                      swap (i2.I1(), i2.I2());
                      if (edges.Used(i2))
                        found = 1;
                    }

                  if (found || foundrev)
                    {
                      if (foundrev)
                        swap (el.PNum(2), el.PNum(3));

                      changed = 1;
                      for (j = 1; j <= 3; j++)
                        {
                          INDEX_2 i2(el.PNumMod(j), el.PNumMod(j+1));
                          edges.Set (i2, 1);
                        }
                      used.Set (i);
                    }
                }
            if (changed)
              haschanged = 1;
          }
        while (changed);


        unused = 0;
        for (i = 1; i <= nse; i++)
          if (!used.Test(i))
            {
              unused = 1;
              const Element2d & tri = SurfaceElement(i);
              for (j = 1; j <= 3; j++)
                {
                  INDEX_2 i2(tri.PNumMod(j), tri.PNumMod(j+1));
                  edges.Set (i2, 1);
                }
              used.Set(i);
              break;
            }
      }
    while (unused);

    if (haschanged)
      timestamp = NextTimeStamp();
  }


  void Mesh :: Split2Tets()
  {
    PrintMessage (1, "Split To Tets");
    bool has_prisms = 0;

    int oldne = GetNE(); 
    for (int i = 1; i <= oldne; i++)
      {
        Element el = VolumeElement(i);

        if (el.GetType() == PRISM)
          {
            // prism, to 3 tets

            // make minimal node to node 1
            int minpi=0;
            PointIndex minpnum;
            minpnum = GetNP() + 1;

            for (int j = 1; j <= 6; j++)
              {
                if (el.PNum(j) < minpnum)
                  {
                    minpnum = el.PNum(j);
                    minpi = j;
                  }
              }

            if (minpi >= 4)
              {
                for (int j = 1; j <= 3; j++)
                  swap (el.PNum(j), el.PNum(j+3));
                minpi -= 3;
              }

            while (minpi > 1)
              {
                int hi = 0;
                for (int j = 0; j <= 3; j+= 3)
                  {
                    hi = el.PNum(1+j);
                    el.PNum(1+j) = el.PNum(2+j);
                    el.PNum(2+j) = el.PNum(3+j);
                    el.PNum(3+j) = hi;
                  }
                minpi--;
              }

            /*
              version 1: edge from pi2 to pi6,
              version 2: edge from pi3 to pi5,
            */

            static const int ntets[2][12] =
              { { 1, 4, 5, 6, 1, 2, 3, 6, 1, 2, 5, 6 },
                { 1, 4, 5, 6, 1, 2, 3, 5, 3, 1, 5, 6 } };

            const int * min2pi;

            if (min2 (el.PNum(2), el.PNum(6)) <
                min2 (el.PNum(3), el.PNum(5)))
              {
                min2pi = &ntets[0][0];
                // (*testout) << "version 1 ";
              }
            else
              {
                min2pi = &ntets[1][0];
                // (*testout) << "version 2 ";
              }


            int firsttet = 1;
            for (int j = 1; j <= 3; j++)
              {
                Element nel(TET);
                for (int k = 1; k <= 4; k++)
                  nel.PNum(k) = el.PNum(min2pi[4 * j + k - 5]);
                nel.SetIndex (el.GetIndex());

                int legal = 1;
                for (int k = 1; k <= 3; k++)
                  for (int l = k+1; l <= 4; l++)
                    if (nel.PNum(k) == nel.PNum(l))
                      legal = 0;

                // (*testout) << nel << " ";
                if (legal)
                  {
                    if (firsttet)
                      {
                        VolumeElement(i) = nel;
                        firsttet = 0;
                      }
                    else
                      {
                        AddVolumeElement(nel);
                      }
                  }
              }
            if (firsttet) cout << "no legal";
            (*testout) << endl;
          }



        else if (el.GetType() == HEX)
          {
            // hex to A) 2 prisms or B) to 5 tets

            // make minimal node to node 1
            int minpi=0;
            PointIndex minpnum;
            minpnum = GetNP() + 1;

            for (int j = 1; j <= 8; j++)
              {
                if (el.PNum(j) < minpnum)
                  {
                    minpnum = el.PNum(j);
                    minpi = j;
                  }
              }

            if (minpi >= 5)
              {
                for (int j = 1; j <= 4; j++)
                  swap (el.PNum(j), el.PNum(j+4));
                minpi -= 4;
              }

            while (minpi > 1)
              {
                int hi = 0;
                for (int j = 0; j <= 4; j+= 4)
                  {
                    hi = el.PNum(1+j);
                    el.PNum(1+j) = el.PNum(2+j);
                    el.PNum(2+j) = el.PNum(3+j);
                    el.PNum(3+j) = el.PNum(4+j);
                    el.PNum(4+j) = hi;
                  }
                minpi--;
              }



            static const int to_prisms[3][12] =
              { { 0, 1, 2, 4, 5, 6, 0, 2, 3, 4, 6, 7 },
                { 0, 1, 5, 3, 2, 6, 0, 5, 4, 3, 6, 7 },
                { 0, 7, 4, 1, 6, 5, 0, 3, 7, 1, 2, 6 },
              };

            const int * min2pi = 0;
            if (min2 (el[4], el[6]) < min2 (el[5], el[7]))
              min2pi = &to_prisms[0][0];
            else if (min2 (el[3], el[6]) < min2 (el[2], el[7]))
              min2pi = &to_prisms[1][0];
            else if (min2 (el[1], el[6]) < min2 (el[2], el[5]))
              min2pi = &to_prisms[2][0];

            if (min2pi)
              {
                has_prisms = 1;
                for (int j = 0; j < 2; j++)
                  {
                    Element nel(PRISM);
                    for (int k = 0; k < 6; k++)
                      nel[k] = el[min2pi[6*j + k]];
                    nel.SetIndex (el.GetIndex());

                    if (j == 0)
                      VolumeElement(i) = nel;
                    else
                      AddVolumeElement(nel);
                  }
              }
            else
              {
                // split to 5 tets

                static const int to_tets[20] =
                  {
                    1, 2, 0, 5,
                    3, 0, 2, 7,
                    4, 5, 7, 0,
                    6, 7, 5, 2,
                    0, 2, 7, 5
                  };

                for (int j = 0; j < 5; j++)
                  {
                    Element nel(TET);
                    for (int k = 0; k < 4; k++)
                      nel[k] = el[to_tets[4*j + k]];
                    nel.SetIndex (el.GetIndex());

                    if (j == 0)
                      VolumeElement(i) = nel;
                    else
                      AddVolumeElement(nel);
                  }

              }
          }





        else if (el.GetType() == PYRAMID)
          {
            // pyramid, to 2 tets

            // cout << "pyramid: " << el << endl;

            static const int ntets[2][8] =
              { { 1, 2, 3, 5, 1, 3, 4, 5 },
                { 1, 2, 4, 5, 4, 2, 3, 5 }};

            const int * min2pi;

            if (min2 (el[0], el[2]) < min2 (el[1], el[3]))
              min2pi = &ntets[0][0];
            else
              min2pi = &ntets[1][0];

            bool firsttet = 1;
            for (int j = 0; j < 2; j++)
              {
                Element nel(TET);
                for (int k = 0; k < 4; k++)
                  nel[k] = el[min2pi[4*j + k]-1];
                nel.SetIndex (el.GetIndex());

                // cout << "pyramid-tet: " << nel << endl;

                bool legal = 1;
                for (int k = 0; k < 3; k++)
                  for (int l = k+1; l < 4; l++)
                    if (nel[k] == nel[l])
                      legal = 0;

                if (legal)
                  {
                    (*testout) << nel << " ";
                    if (firsttet)
                      VolumeElement(i) = nel;
                    else
                      AddVolumeElement(nel);

                    firsttet = 0;
                  }
              }
            if (firsttet) cout << "no legal";
            (*testout) << endl;
          }
      }


    int oldnse = GetNSE(); 
    for (int i = 1; i <= oldnse; i++)
      {
        Element2d el = SurfaceElement(i);
        if (el.GetNP() == 4)
          {
            (*testout) << "split el: " << el << " to ";

            static const int ntris[2][6] =
              { { 1, 2, 3, 1, 3, 4 },
                { 1, 2, 4, 4, 2, 3 }};

            const int * min2pi;

            if (min2 (el.PNum(1), el.PNum(3)) <
                min2 (el.PNum(2), el.PNum(4)))
              min2pi = &ntris[0][0];
            else
              min2pi = &ntris[1][0];

            for (int j = 0; j <6; j++)
              (*testout) << min2pi[j] << " ";


            int firsttri = 1;
            for (int j = 1; j <= 2; j++)
              {
                Element2d nel(3);
                for (int k = 1; k <= 3; k++)
                  nel.PNum(k) = el.PNum(min2pi[3 * j + k - 4]);
                nel.SetIndex (el.GetIndex());

                int legal = 1;
                for (int k = 1; k <= 2; k++)
                  for (int l = k+1; l <= 3; l++)
                    if (nel.PNum(k) == nel.PNum(l))
                      legal = 0;

                if (legal)
                  {
                    (*testout) << nel << " ";
                    if (firsttri)
                      {
                        SurfaceElement(i) = nel;
                        firsttri = 0;
                      }
                    else
                      {
                        AddSurfaceElement(nel);
                      }
                  }
              }
            (*testout) << endl;

          }
      }


    if (has_prisms)

      Split2Tets();

    else
      {
        for (int i = 1; i <= GetNE(); i++)
          {
            Element & el = VolumeElement(i);
            const Point3d & p1 = Point (el.PNum(1));
            const Point3d & p2 = Point (el.PNum(2));
            const Point3d & p3 = Point (el.PNum(3));
            const Point3d & p4 = Point (el.PNum(4));

            double vol = (Vec3d (p1, p2) * 
                          Cross (Vec3d (p1, p3), Vec3d(p1, p4)));
            if (vol > 0)
              swap (el.PNum(3), el.PNum(4));
          }



        UpdateTopology();
        timestamp = NextTimeStamp();
      }

    RebuildSurfaceElementLists();
  }

  void Mesh :: BuildElementSearchTree ()
  {
    if (elementsearchtreets == GetTimeStamp()) return;

    {
      std::lock_guard<std::mutex> guard(buildsearchtree_mutex);
      if (elementsearchtreets != GetTimeStamp())
        {
          NgLock lock(mutex);
          lock.Lock();
          
          PrintMessage (4, "Rebuild element searchtree");
          
          elementsearchtree = nullptr;
          
          int ne = (dimension == 2) ? GetNSE() : GetNE();
          if (dimension == 3 && !GetNE() && GetNSE())
            ne = GetNSE();

          if (ne) 
            {
              if (dimension == 2 || (dimension == 3 && !GetNE()) )
                {
                  Box<3> box (Box<3>::EMPTY_BOX);
                  for (SurfaceElementIndex sei = 0; sei < ne; sei++)
                    // box.Add (points[surfelements[sei].PNums()]);
                    for (auto pi : surfelements[sei].PNums())
                      box.Add (points[pi]);
                  
                  box.Increase (1.01 * box.Diam());
                  elementsearchtree = make_unique<BoxTree<3>> (box);
                  
                  for (SurfaceElementIndex sei = 0; sei < ne; sei++)
                    {
                      //  box.Set (points[surfelements[sei].PNums()]);

                      Box<3> box (Box<3>::EMPTY_BOX);
                      for (auto pi : surfelements[sei].PNums())
                        box.Add (points[pi]);
                      
                      elementsearchtree -> Insert (box, sei+1);
                    }
                }
              else
                {
                  Box<3> box (Box<3>::EMPTY_BOX);
                  for (ElementIndex ei = 0; ei < ne; ei++)
                    // box.Add (points[volelements[ei].PNums()]);
                    for (auto pi : volelements[ei].PNums())
                      box.Add (points[pi]);
                  
                  box.Increase (1.01 * box.Diam());
                  elementsearchtree = make_unique<BoxTree<3>> (box);
                  
                  for (ElementIndex ei = 0; ei < ne; ei++)
                    {
                      // box.Set (points[volelements[ei].PNums()]);

                      Box<3> box (Box<3>::EMPTY_BOX);
                      for (auto pi : volelements[ei].PNums())
                        box.Add (points[pi]);

                      auto & el = volelements[ei];
                      if(el.IsCurved())
                        box.Increase(1.2*box.Diam());


                      elementsearchtree -> Insert (box, ei+1);
                    }
                }
              
              elementsearchtreets = GetTimeStamp();
            }
        }
    }
  }

  
  int SolveLinearSystemLS (const Vec3d & col1,
                           const Vec3d & col2,
                           const Vec3d & rhs,
                           Vec2d & sol)
  {
    double a11 = col1 * col1;
    double a12 = col1 * col2;
    double a22 = col2 * col2;
    
    double det = a11 * a22 - a12 * a12;
    
    if (det*det <= 1e-24 * a11 * a22)
      {
        sol = Vec2d (0, 0);
        return 1;
      }
    
    Vec2d aTrhs;
    aTrhs.X() = col1*rhs;
    aTrhs.Y() = col2*rhs;

    sol.X() = ( a22 * aTrhs.X() - a12 * aTrhs.Y()) / det;
    sol.Y() = (-a12 * aTrhs.X() + a11 * aTrhs.Y()) / det;
    return 0;
  }

  bool ValidBarCoord(double lami[3], double eps=1e-12)
  {
    return (lami[0]<=1.+eps && lami[0]>=0.-eps && lami[1]<=1.+eps && lami[1]>=0.-eps && lami[2]<=1.+eps && lami[2]>=0.-eps );
  }

  bool Mesh :: PointContainedIn2DElement(const Point3d & p,
                                         double lami[3],
                                         const int element,
                                         bool consider3D) const
  {
    Vec3d col1, col2, col3;
    Vec3d rhs, sol;
    const double eps = 1e-6;

    NgArray<Element2d> loctrigs;

    
    //SZ 
    if(SurfaceElement(element).GetType()==QUAD)
      {
        const Element2d & el = SurfaceElement(element); 

        const Point3d & p1 = Point(el.PNum(1)); 
        const Point3d & p2 = Point(el.PNum(2));
        const Point3d & p3 = Point(el.PNum(3));
        const Point3d & p4 = Point(el.PNum(4));

        // Coefficients of Bilinear Mapping from Ref-Elem to global Elem
        // X = a + b x + c y + d x y 
        Vec3d a = p1; 
        Vec3d b = p2 - a; 
        Vec3d c = p4 - a; 
        Vec3d d = p3 - a - b - c;

        /*cout << "p = " << p << endl;
        cout << "p1 = " << p1 << endl;
        cout << "p2 = " << p2 << endl;
        cout << "p3 = " << p3 << endl;
        cout << "p4 = " << p4 << endl;

        cout << "a = " << a << endl;
        cout << "b = " << b << endl;
        cout << "c = " << c << endl;
        cout << "d = " << d << endl;*/


        Vec3d pa = p-a;
        double dxb = d.X()*b.Y()-d.Y()*b.X();
        double dxc = d.X()*c.Y()-d.Y()*c.X();
        double bxc = b.X()*c.Y()-b.Y()*c.X();
        double bxpa = b.X()*pa.Y()-b.Y()*pa.X();
        double cxpa = c.X()*pa.Y()-c.Y()*pa.X();
        double dxpa = d.X()*pa.Y()-d.Y()*pa.X();

        /*cout << "dxb = " << dxb << endl;
        cout << "dxc = " << dxc << endl;
        cout << "bxc = " << bxc << endl;
        cout << "bxpa = " << bxpa << endl;
        cout << "cxpa = " << cxpa << endl;
        cout << "dxpa = " << dxpa << endl;*/

        /*
          P = a + b x + c y + d x y
          1) P1 = a1 + b1 x + c1 y + d1 x y
          2) P2 = a2 + b2 x + c2 y + d2 x y
          
          -> det(x,d) = det(a,d) + det(b,d) x + det(c,d) y
            -> x = 1/det(b,d) *( det(P-a,d)-det(c,d) y )
            -> y = 1/det(c,d) *( det(P-a,d)-det(b,d) x )
          
          -> x = (P1 - a1 - c1 y)/(b1 + d1 y)
            -> det(c,d) y**2 + [det(d,P-a) + det(c,b)] y + det(b,P-a) = 0
          ( same if we express x = (P2 - a2 - c2 y)/(b2 + d2 y) )

          -> y = (P1 - a1 - b1 x)/(c1 + d1 x)
            -> det(b,d) x**2 + [det(d,P-a) + det(b,c)] x + det(c,P-a) = 0
          ( same if we express y = (P2 - a2 - b2 x)/(c2 + d2 x)
         */

        lami[2]=0.; 
        double eps = 1.E-12;
        double c1,c2,r;

        //First check if point is "exactly" a vertex point
        Vec3d d1 = p-p1;
        Vec3d d2 = p-p2;
        Vec3d d3 = p-p3;
        Vec3d d4 = p-p4;

        //cout << " d1 = " << d1 << ", d2 = " << d2 << ", d3 = " << d3 << ", d4 = " << d4 << endl;
        
        if (d1.Length2() < sqr(eps)*d2.Length2() && d1.Length2() < sqr(eps)*d3.Length2() && d1.Length2() < sqr(eps)*d4.Length2())
          {
            lami[0] = lami[1] = 0.;
            return true;
          }
        else if (d2.Length2() < sqr(eps)*d1.Length2() && d2.Length2() < sqr(eps)*d3.Length2() && d2.Length2() < sqr(eps)*d4.Length2())
          {
            lami[0] = 1.;
            lami[1] = 0.;
            return true;
          }
        else if (d3.Length2() < sqr(eps)*d1.Length2() && d3.Length2() < sqr(eps)*d2.Length2() && d3.Length2() < sqr(eps)*d4.Length2())
          {
            lami[0] = lami[1] = 1.;
            return true;
          }
        else if (d4.Length2() < sqr(eps)*d1.Length2() && d4.Length2() < sqr(eps)*d2.Length2() && d4.Length2() < sqr(eps)*d3.Length2())
          {
            lami[0] = 0.;
            lami[1] = 1.;
            return true;
          }//if d is nearly 0: solve resulting linear system
        else if (d.Length2() < sqr(eps)*b.Length2() && d.Length2() < sqr(eps)*c.Length2())
          {
            Vec2d sol;
            SolveLinearSystemLS (b, c, p-a, sol);
            lami[0] = sol.X();
            lami[1] = sol.Y();
	    return ValidBarCoord(lami, eps);
          }// if dxc is nearly 0: solve resulting linear equation for y and compute x
        else if (fabs(dxc) < sqr(eps))
          {
            lami[1] = -bxpa/(dxpa-bxc);
            lami[0] = (dxpa-dxc*lami[1])/dxb;
            return ValidBarCoord(lami, eps);
          }// if dxb is nearly 0: solve resulting linear equation for x and compute y
        else if (fabs(dxb) < sqr(eps))
          {
            lami[0] = -cxpa/(dxpa+bxc);
            lami[1] = (dxpa-dxb*lami[0])/dxc;
            return ValidBarCoord(lami, eps);
          }//if dxb >= dxc: solve quadratic equation in y and compute x
        else if (fabs(dxb) >= fabs(dxc))
          {
            c1 = (bxc-dxpa)/dxc;
            c2 = -bxpa/dxc;
            r = c1*c1/4.0-c2;

            //quadratic equation has only 1 (unstable) solution
            if (fabs(r) < eps) //not eps^2!
              {
                lami[1] = -c1/2;
                lami[0] = (dxpa-dxc*lami[1])/dxb;
                return ValidBarCoord(lami, eps);
              }
            if (r < 0) return false;

            lami[1] = -c1/2+sqrt(r);
            lami[0] = (dxpa-dxc*lami[1])/dxb;

            if (ValidBarCoord(lami, eps))
                return true;
            else
              {
                lami[1] = -c1/2-sqrt(r);
                lami[0] = (dxpa-dxc*lami[1])/dxb;
                return ValidBarCoord(lami, eps);
              }
          }//if dxc > dxb: solve quadratic equation in x and compute y
        else
          {
            c1 = (-bxc-dxpa)/dxb;
            c2 = -cxpa/dxb;
            r = c1*c1/4.0-c2;

            //quadratic equation has only 1 (unstable) solution
            if (fabs(r) < eps) //not eps^2!
              {
                lami[0] = -c1/2;
                lami[1] = (dxpa-dxb*lami[0])/dxc;
                return ValidBarCoord(lami, eps);
              }
            if (r < 0) return false;

            lami[0] = -c1/2+sqrt(r);
            lami[1] = (dxpa-dxb*lami[0])/dxc;

            if (ValidBarCoord(lami, eps))
                return true;
            else
              {
                lami[0] = -c1/2-sqrt(r);
                lami[1] = (dxpa-dxb*lami[0])/dxc;
                return ValidBarCoord(lami, eps);
              }
          }
        
        /*
        double dxa = d.X()*a.Y()-d.Y()*a.X(); 
        double dxp = d.X()*p.Y()-d.Y()*p.X();
	
	
        double c0,c1,c2; // ,rt; 
        

	Vec3d dp13 = p3-p1;
	Vec3d dp24 = p4-p2;
	double d1 = dp13.Length2();
	double d2 = dp24.Length2();

	// if(fabs(d.X()) <= eps && fabs(d.Y())<= eps)
	//if (d.Length2() < sqr(eps))
        if (d.Length2() < sqr(eps)*d1 && d.Length2() < sqr(eps)*d2)
          {
	    //Solve Linear System
	    Vec2d sol;
            SolveLinearSystemLS (b, c, p-a, sol);
            lami[0] = sol.X();
            lami[1] = sol.Y();

	    if(lami[1]<=1.+eps && lami[1]>=0.-eps && lami[0]<=1.+eps && lami[0]>=0.-eps)
	      return true;
	    
            
              //lami[0]=(c.Y()*(p.X()-a.X())-c.X()*(p.Y()-a.Y()))/
              //(b.X()*c.Y() -b.Y()*c.X()); 
            //lami[1]=(-b.Y()*(p.X()-a.X())+b.X()*(p.Y()-a.Y()))/
             // (b.X()*c.Y() -b.Y()*c.X()); 
            
          } 
        else
          if(fabs(dxb) <= eps*fabs(dxc))
            {
	      lami[1] = (dxp-dxa)/dxc;
              if(fabs(b.X()+d.X()*lami[1])>=fabs(b.Y()+d.Y()*lami[1]))
                lami[0] = (p.X()-a.X() - c.X()*lami[1])/(b.X()+d.X()*lami[1]); 
              else
                lami[0] = (p.Y()-a.Y() - c.Y()*lami[1])/(b.Y()+d.Y()*lami[1]);

	      if(lami[1]<=1.+eps && lami[1]>=0.-eps && lami[0]<=1.+eps && lami[0]>=0.-eps)
		return true;
            }
          else
            if(fabs(dxc) <= eps*fabs(dxb))
              {
		lami[0] = (dxp-dxa)/dxb;
                if(fabs(c.X()+d.X()*lami[0])>=fabs(c.Y()+d.Y()*lami[0]))
                  lami[1] = (p.X()-a.X() - b.X()*lami[0])/(c.X()+d.X()*lami[0]); 
                else
                  lami[1] = (p.Y()-a.Y() - b.Y()*lami[0])/(c.Y()+d.Y()*lami[0]);

		if(lami[1]<=1.+eps && lami[1]>=0.-eps && lami[0]<=1.+eps && lami[0]>=0.-eps)
		  return true;
              }
            else //Solve quadratic equation
              {
		c2 = -d.X()*dxb;
		c1 = b.X()*dxc - c.X()*dxb + d.X()*(dxp-dxa);
		c0 = c.X()*(dxp-dxa) + (a.X()-p.X())*dxc;
		double rt =  c1*c1 - 4*c2*c0;
		
		if (rt < 0.) return false; 
		lami[1] = (-c1 + sqrt(rt))/2/c2;


		if(lami[1]<=1.+eps && lami[1]>=0.-eps)
		  {
		    lami[0] = (dxp - dxa -dxb*lami[1])/dxc;
		    
		    if(lami[0]<=1.+eps && lami[0]>=0.-eps)
		      return true;
		  }
		lami[1] = (-c1 - sqrt(rt))/2/c2;

		lami[0] = (dxp - dxa -dxb*lami[1])/dxc;

		if(lami[1]<=1.+eps && lami[1]>=0.-eps && lami[0]<=1.+eps && lami[0]>=0.-eps)
		  return true;

		c2 = d.Y()*dxb;
		c1 = b.Y()*dxc - c.Y()*dxb + d.Y()*(dxp-dxa);
		c0 = c.Y()*(dxp -dxa) + (a.Y()-p.Y())*dxc;
		rt =  c1*c1 - 4*c2*c0;
		
		if (rt < 0.) return false; 
		lami[1] = (-c1 + sqrt(rt))/2/c2;

		if(lami[1]<=1.+eps && lami[1]>=0.-eps)
		  {
		    lami[0] = (dxp - dxa -dxb*lami[1])/dxc;

		    if(lami[0]<=1.+eps && lami[0]>=0.-eps)
		      return true;
		  }
		lami[1] = (-c1 - sqrt(rt))/2/c2;

		lami[0] = (dxp - dxa -dxb*lami[1])/dxc;

		if(lami[1]<=1.+eps && lami[1]>=0.-eps && lami[0]<=1.+eps && lami[0]>=0.-eps)
		  return true;

		c2 = -d.X()*dxc;
		c1 = -b.X()*dxc + c.X()*dxb + d.X()*(dxp-dxa);
		c0 = b.X()*(dxp -dxa) + (a.X()-p.X())*dxb;
		rt =  c1*c1 - 4*c2*c0;
		
		if (rt < 0.) return false; 
		lami[1] = (-c1 + sqrt(rt))/2/c2;

		if(lami[1]<=1.+eps && lami[1]>=0.-eps)
		  {
		    lami[0] = (dxp - dxa -dxc*lami[1])/dxb;

		    if(lami[0]<=1.+eps && lami[0]>=0.-eps)
		      return true;
		  }
		lami[1] = (-c1 - sqrt(rt))/2/c2;

		lami[0] = (dxp - dxa -dxc*lami[1])/dxb;

		if(lami[1]<=1.+eps && lami[1]>=0.-eps && lami[0]<=1.+eps && lami[0]>=0.-eps)
		  return true;
                  }*/

      
        //cout << "lam0,1 = " << lami[0] << ", " << lami[1] << endl;
        
        /*if( lami[0] <= 1.+eps  && lami[0] >= -eps && lami[1]<=1.+eps && lami[1]>=-eps)
          {
            if(consider3D)
              {
                Vec3d n = Cross(b,c);
                lami[2] = 0;
                for(int i=1; i<=3; i++)
                  lami[2] +=(p.X(i)-a.X(i)-lami[0]*b.X(i)-lami[1]*c.X(i)) * n.X(i);
                if(lami[2] >= -eps && lami[2] <= eps)
                  return true;
              }
            else
              return true;
	      }*/

        return false;

      }
    else
      {
        //	  SurfaceElement(element).GetTets (loctets);
        loctrigs.SetSize(1);
        loctrigs.Elem(1) = SurfaceElement(element);



        for (int j = 1; j <= loctrigs.Size(); j++)
          {
            const Element2d & el = loctrigs.Get(j);


            const Point3d & p1 = Point(el.PNum(1));
            const Point3d & p2 = Point(el.PNum(2));
            const Point3d & p3 = Point(el.PNum(3));
            /*
              Box3d box;
              box.SetPoint (p1);
              box.AddPoint (p2);
              box.AddPoint (p3);
              box.AddPoint (p4);
              if (!box.IsIn (p))
              continue;
            */
            col1 = p2-p1;
            col2 = p3-p1;
            col3 = Cross(col1,col2);
            //col3 = Vec3d(0, 0, 1);
            rhs = p - p1;

            // int retval = 
            SolveLinearSystem (col1, col2, col3, rhs, sol);

            //(*testout) << "retval " << retval << endl;

            //(*testout) << "col1 " << col1 << " col2 " << col2 << " col3 " << col3 << " rhs " << rhs << endl;
            //(*testout) << "sol " << sol << endl;

            if (SurfaceElement(element).GetType() ==TRIG6 || curvedelems->IsSurfaceElementCurved(element-1))
              {
                netgen::Point<2> lam(1./3,1./3);
                Vec<3> rhs;
                Vec<2> deltalam;
                netgen::Point<3> x;
                Mat<3,2> Jac,Jact;
                
                double delta=1;
                
                bool retval;
                
                int i = 0;
                
                const int maxits = 30;
                while(delta > 1e-16 && i<maxits)
                  {
                    curvedelems->CalcSurfaceTransformation(lam,element-1,x,Jac);
                    rhs = p-x;
                    Jac.Solve(rhs,deltalam);
                    
                    lam += deltalam;
                    
                    delta = deltalam.Length2();
                    
                    i++;
                    //(*testout) << "pcie i " << i << " delta " << delta << " p " << p << " x " << x << " lam " << lam << endl;
                    //<< "Jac " << Jac << endl;
                  }
                
                if(i==maxits)
                  return false;
                
                sol.X() = lam(0);
                sol.Y() = lam(1);

                if (SurfaceElement(element).GetType() !=TRIG6 )
                  {
                    sol.Z() = sol.X();
                    sol.X() = sol.Y();
                    sol.Y() = 1.0 - sol.Z() - sol.X();
                  }

              }
            if (sol.X() >= -eps && sol.Y() >= -eps && 
                sol.X() + sol.Y() <= 1+eps)
              {
                if(!consider3D || (sol.Z() >= -eps && sol.Z() <= eps))
                  {
                    lami[0] = sol.X();
                    lami[1] = sol.Y();
                    lami[2] = sol.Z();

                    return true;
                  }
              }
          }
      }

    return false;

  }




  bool Mesh :: PointContainedIn3DElement(const Point3d & p,
                                         double lami[3],
                                         const int element) const
  {
    //bool oldresult = PointContainedIn3DElementOld(p,lami,element);
    //(*testout) << "old result: " << oldresult
    //       << " lam " << lami[0] << " " << lami[1] << " " << lami[2] << endl;

    //if(!curvedelems->IsElementCurved(element-1))
    //  return PointContainedIn3DElementOld(p,lami,element);


    const double eps = 1.e-4;
    const Element & el = VolumeElement(element);

    netgen::Point<3> lam = 0.0;

    if (el.GetType() == TET || el.GetType() == TET10)
      {
        lam = 0.25;
      }
    else if (el.GetType() == PRISM)
      {
        lam(0) = 0.33; lam(1) = 0.33; lam(2) = 0.5;
      }
    else if (el.GetType() == PYRAMID)
      {
        lam(0) = 0.4; lam(1) = 0.4; lam(2) = 0.2;
      }
    else if (el.GetType() == HEX)
      {
        lam = 0.5;
      }


    Vec<3> deltalam,rhs;
    netgen::Point<3> x;
    Mat<3,3> Jac,Jact;

    double delta=1;

    bool retval;

    int i = 0;

    const int maxits = 30;
    while(delta > 1e-16 && i<maxits)
      {
        curvedelems->CalcElementTransformation(lam,element-1,x,Jac);
        rhs = p-x;
        Jac.Solve(rhs,deltalam);

        lam += deltalam;

        delta = deltalam.Length2();

        i++;
        //(*testout) << "pcie i " << i << " delta " << delta << " p " << p << " x " << x << " lam " << lam << endl;
        //<< "Jac " << Jac << endl;
      }

    if(i==maxits)
      return false;

    for(i=0; i<3; i++)
      lami[i] = lam(i);



    if (el.GetType() == TET || el.GetType() == TET10)
      {
        retval = (lam(0) > -eps && 
                  lam(1) > -eps && 
                  lam(2) > -eps && 
                  lam(0) + lam(1) + lam(2) < 1+eps);
      }
    else if (el.GetType() == PRISM || el.GetType() == PRISM15)
      {
        retval = (lam(0) > -eps &&
                  lam(1) > -eps &&
                  lam(2) > -eps &&
                  lam(2) < 1+eps &&
                  lam(0) + lam(1) < 1+eps);
      }
    else if (el.GetType() == PYRAMID || el.GetType() == PYRAMID13)
      {
        retval = (lam(0) > -eps &&
                  lam(1) > -eps &&
                  lam(2) > -eps &&
                  lam(0) + lam(2) < 1+eps &&
                  lam(1) + lam(2) < 1+eps);
      }
    else if (el.GetType() == HEX || el.GetType() == HEX20)
      {
        retval = (lam(0) > -eps && lam(0) < 1+eps &&
                  lam(1) > -eps && lam(1) < 1+eps &&
                  lam(2) > -eps && lam(2) < 1+eps);
      }
    else
      throw NgException("Da haun i wos vagessn");

    return retval;
  }



  bool Mesh :: PointContainedIn3DElementOld(const Point3d & p,
                                            double lami[3],
                                            const int element) const
  {
    Vec3d col1, col2, col3;
    Vec3d rhs, sol;
    const double eps = 1.e-4;

    NgArray<Element> loctets;

    VolumeElement(element).GetTets (loctets);

    for (int j = 1; j <= loctets.Size(); j++)
      {
        const Element & el = loctets.Get(j);

        const Point3d & p1 = Point(el.PNum(1));
        const Point3d & p2 = Point(el.PNum(2));
        const Point3d & p3 = Point(el.PNum(3));
        const Point3d & p4 = Point(el.PNum(4));

        Box3d box;
        box.SetPoint (p1);
        box.AddPoint (p2);
        box.AddPoint (p3);
        box.AddPoint (p4);
        if (!box.IsIn (p))
          continue;

        col1 = p2-p1;
        col2 = p3-p1;
        col3 = p4-p1;
        rhs = p - p1;

        SolveLinearSystem (col1, col2, col3, rhs, sol);

        if (sol.X() >= -eps && sol.Y() >= -eps && sol.Z() >= -eps &&
            sol.X() + sol.Y() + sol.Z() <= 1+eps)
          {
            NgArray<Element> loctetsloc;
            NgArray<netgen::Point<3> > pointsloc;

            VolumeElement(element).GetTetsLocal (loctetsloc);
            VolumeElement(element).GetNodesLocalNew (pointsloc);

            const Element & le = loctetsloc.Get(j);


            Point3d pp = 
              pointsloc.Get(le.PNum(1)) 
              + sol.X() * Vec3d (pointsloc.Get(le.PNum(1)), pointsloc.Get(le.PNum(2))) 
              + sol.Y() * Vec3d (pointsloc.Get(le.PNum(1)), pointsloc.Get(le.PNum(3))) 
              + sol.Z() * Vec3d (pointsloc.Get(le.PNum(1)), pointsloc.Get(le.PNum(4))) ;

            lami[0] = pp.X();
            lami[1] = pp.Y();
            lami[2] = pp.Z();
            return true;
          }
      }
    return false;
  }


  int Mesh :: GetElementOfPoint (const netgen::Point<3> & p,
                                 double lami[3],
                                 bool build_searchtree,
                                 const int index,
                                 const bool allowindex) const
  {
    if(index != -1) 
      {
        NgArray<int> dummy(1);
        dummy[0] = index;
        return GetElementOfPoint(p,lami,&dummy,build_searchtree,allowindex);
      }
    else
      return GetElementOfPoint(p,lami,NULL,build_searchtree,allowindex);
  }




  int Mesh :: GetElementOfPoint (const netgen::Point<3> & p,
                                 double lami[3],
                                 const NgArray<int> * const indices,
                                 bool build_searchtree,
                                 const bool allowindex) const
  {
    // const double pointtol = 1e-12;
    // netgen::Point<3> pmin = p - Vec<3> (pointtol, pointtol, pointtol);
    // netgen::Point<3> pmax = p + Vec<3> (pointtol, pointtol, pointtol);

    if ( (dimension == 2 && !GetNSE()) ||
    	 (dimension == 3 && !GetNE() && !GetNSE()) )
      return -1;

    if (dimension == 2 || (dimension==3 && !GetNE() && GetNSE()))
      {
        int ne;
        int ps_startelement = 0;  // disable global buffering

        if(ps_startelement != 0 && ps_startelement <= GetNSE() && PointContainedIn2DElement(p,lami,ps_startelement))
          return ps_startelement;

        NgArray<int> locels;
        if (elementsearchtree || build_searchtree)
          {
            // update if necessary:
            const_cast<Mesh&>(*this).BuildElementSearchTree ();
            // double tol = elementsearchtree->Tolerance();
            // netgen::Point<3> pmin = p - Vec<3> (tol, tol, tol);
            // netgen::Point<3> pmax = p + Vec<3> (tol, tol, tol);
            elementsearchtree->GetIntersecting (p, p, locels);
            ne = locels.Size();
          }
        else
          ne = GetNSE();

        for (int i = 1; i <= ne; i++)
          {
            int ii;

            if (elementsearchtree)
              ii = locels.Get(i);
            else
              ii = i;

            if(ii == ps_startelement) continue;

            if(indices != NULL && indices->Size() > 0)
              {
                bool contained = indices->Contains(SurfaceElement(ii).GetIndex());
                if((allowindex && !contained) || (!allowindex && contained)) continue;
              }

            if(PointContainedIn2DElement(p,lami,ii)) return ii;

          }
        return 0;
      }
    else

      {
        int ps_startelement = 0;  // disable global buffering
        // int i, j;
        int ne;

        if(ps_startelement != 0 && PointContainedIn3DElement(p,lami,ps_startelement))
          return ps_startelement;

        NgArray<int> locels;
        if (elementsearchtree || build_searchtree)
          {
            // update if necessary:
            const_cast<Mesh&>(*this).BuildElementSearchTree ();
            // double tol = elementsearchtree->Tolerance();            
            // netgen::Point<3> pmin = p - Vec<3> (tol, tol, tol);
            // netgen::Point<3> pmax = p + Vec<3> (tol, tol, tol);
            elementsearchtree->GetIntersecting (p, p, locels);
            ne = locels.Size();
          }
        else
          ne = GetNE();

        for (int i = 1; i <= ne; i++)
          {
            int ii;

            if (elementsearchtree)
              ii = locels.Get(i);
            else
              ii = i;
            if(ii == ps_startelement) continue;

            if(indices != NULL && indices->Size() > 0)
              {
                bool contained = indices->Contains(VolumeElement(ii).GetIndex());
                if((allowindex && !contained) || (!allowindex && contained)) continue;
              }

            if(PointContainedIn3DElement(p,lami,ii)) 
              {
                ps_startelement = ii;
                return ii;
              }
          }

        // Not found, try uncurved variant:
        for (int i = 1; i <= ne; i++)
          {
            int ii;

            if (elementsearchtree)
              ii = locels.Get(i);
            else
              ii = i;

            if(indices != NULL && indices->Size() > 0)
              {
                bool contained = indices->Contains(VolumeElement(ii).GetIndex());
                if((allowindex && !contained) || (!allowindex && contained)) continue;
              }


            if(PointContainedIn3DElementOld(p,lami,ii)) 
              {
                ps_startelement = ii;
                (*testout) << "WARNING: found element of point " << p <<" only for uncurved mesh" << endl;
                return ii;
              }
          }


        return 0;
      }
  }



  int Mesh :: GetSurfaceElementOfPoint (const netgen::Point<3> & p,
                                        double lami[3],
                                        bool build_searchtree,
                                        const int index,
                                        const bool allowindex) const
  {
    if(index != -1) 
      {
        NgArray<int> dummy(1);
        dummy[0] = index;
        return GetSurfaceElementOfPoint(p,lami,&dummy,build_searchtree,allowindex);
      }
    else
      return GetSurfaceElementOfPoint(p,lami,NULL,build_searchtree,allowindex);
  }




  int Mesh :: GetSurfaceElementOfPoint (const netgen::Point<3> & p,
                                        double lami[3],
                                        const NgArray<int> * const indices,
                                        bool build_searchtree,
                                        const bool allowindex) const
  {
    if (dimension == 2)
      {
        double vlam[3];
        int velement = GetElementOfPoint(p, vlam, NULL, build_searchtree, allowindex);
        if(velement == 0)
          return 0;

        vlam[2] = 1.-vlam[0] - vlam[1];
        NgArray<int> edges;
        topology.GetSurfaceElementEdges(velement, edges);
        Array<SegmentIndex> segs(edges.Size());
        for(auto i : Range(edges))
          segs[i] = topology.GetSegmentOfEdge(edges[i]);

        for(auto i : Range(segs))
          {
            if(IsInvalid(segs[i]))
              continue;
            auto& el = SurfaceElement(velement);
            if(el.GetType() == TRIG)
              {
                double seg_lam;
                double lam;
                auto seg = LineSegment(segs[i]);
                    for(auto k : Range(3))
                      {
                        if(seg[0] == el[k])
                          lam = vlam[k];
                        if(seg[1] == el[k])
                          seg_lam = vlam[k];
                      }
                if(1.- seg_lam - lam < 1e-5)
                  {
                    // found point close to segment -> use barycentric coordinates directly
                    lami[0] = lam;
                    return int(segs[i])+1;
                  }
              }
            else
              throw NgException("Quad not implemented yet!");
          }

        return 0;
        throw NgException("GetSurfaceElementOfPoint not yet implemented for 2D meshes");
      }
    else
      {
        double vlam[3];
        int velement = GetElementOfPoint(p,vlam,NULL,build_searchtree,allowindex);

        //(*testout) << "p " << p << endl;
        //(*testout) << "velement " << velement << endl;

        if (!GetNE() && GetNSE() )
          {
            lami[0] = vlam[0];
            lami[1] = vlam[1];
            lami[2] = vlam[2];
            return velement;
          }
        
        NgArray<int> faces;
        topology.GetElementFaces(velement,faces);

        //(*testout) << "faces " << faces << endl;

        for(int i=0; i<faces.Size(); i++)
          faces[i] = topology.GetFace2SurfaceElement(faces[i]);

        //(*testout) << "surfel " << faces << endl;

        for(int i=0; i<faces.Size(); i++)
          {
            if(faces[i] == 0)
              continue;

            auto & el = VolumeElement(velement);

            if (el.GetType() == TET)
            {
              double lam4[4] = { vlam[0], vlam[1], vlam[2], 1.0-vlam[0]-vlam[1]-vlam[2] };
              double face_lam = lam4[i];
              if(face_lam < 1e-5)
              {
                // found volume point very close to a face -> use barycentric coordinates directly
                lami[2] = 0.0;
                auto sel = SurfaceElement(faces[i]);

                for(auto j : Range(1,3))
                  for(auto k : Range(4))
                    if(sel[j] == el[k])
                      lami[j-1] = lam4[k]/(1.0-face_lam);
                return faces[i];
              }
            }

            if(indices && indices->Size() != 0)
              {
                if(indices->Contains(SurfaceElement(faces[i]).GetIndex()) &&
                   PointContainedIn2DElement(p,lami,faces[i],true))
                  return faces[i];
              }
            else
              {
                if(PointContainedIn2DElement(p,lami,faces[i],true))
                  {
                    //(*testout) << "found point " << p << " in sel " << faces[i]
                    //	       << ", lam " << lami[0] << ", " << lami[1] << ", " << lami[2] << endl;
                    return faces[i];
                  }
              }
          }

        NgArray<int> faces2;
        topology.GetElementFaces(velement,faces2);
        /*
        cout << "no matching surf element" << endl
             << "p = " << p << endl
             << "faces-orig = " << faces2 << endl
             << "faces = " << faces << endl
             << ", vol el = " << velement
             << ", vlam = " << vlam[0] << "," << vlam[1] << "," << vlam[2] << endl;
        */
      }

    return 0;
  }


  void Mesh::GetIntersectingVolEls(const Point3d& p1, const Point3d& p2, 
                                   NgArray<int> & locels) const
  {
    elementsearchtree->GetIntersecting (p1, p2, locels);
  }

  void Mesh :: SplitIntoParts()
  {
    int i, j, dom;
    int ne = GetNE();
    int np = GetNP();
    int nse = GetNSE();

    NgBitArray surfused(nse);
    NgBitArray pused (np);

    surfused.Clear();

    dom = 0;

    while (1)
      {
        int cntd = 1;

        dom++;

        pused.Clear();

        int found = 0;
        for (i = 1; i <= nse; i++)
          if (!surfused.Test(i))
            {
              SurfaceElement(i).SetIndex (dom);
              for (j = 1; j <= 3; j++)
                pused.Set (SurfaceElement(i).PNum(j));
              found = 1;
              cntd = 1;
              surfused.Set(i);
              break;
            }

        if (!found)
          break;

        int change;
        do
          {
            change = 0;
            for (i = 1; i <= nse; i++)
              {
                int is = 0, isnot = 0;
                for (j = 1; j <= 3; j++)
                  if (pused.Test(SurfaceElement(i).PNum(j)))
                    is = 1;
                  else
                    isnot = 1;

                if (is && isnot)
                  {
                    change = 1;
                    for (j = 1; j <= 3; j++)
                      pused.Set (SurfaceElement(i).PNum(j));
                  }

                if (is) 
                  {
                    if (!surfused.Test(i))
                      {
                        surfused.Set(i);
                        SurfaceElement(i).SetIndex (dom);
                        cntd++;
                      }
                  }
              }


            for (i = 1; i <= ne; i++)
              {
                int is = 0, isnot = 0;
                for (j = 1; j <= 4; j++)
                  if (pused.Test(VolumeElement(i).PNum(j)))
                    is = 1;
                  else
                    isnot = 1;

                if (is && isnot)
                  {
                    change = 1;
                    for (j = 1; j <= 4; j++)
                      pused.Set (VolumeElement(i).PNum(j));
                  }

                if (is)
                  {
                    VolumeElement(i).SetIndex (dom);
                  }
              }
          }
        while (change);

        PrintMessage (3, "domain ", dom, " has ", cntd, " surfaceelements");
      }

    /*
      facedecoding.SetSize (dom);
      for (i = 1; i <= dom; i++)
      {
      facedecoding.Elem(i).surfnr = 0;
      facedecoding.Elem(i).domin = i;
      facedecoding.Elem(i).domout = 0;
      }
    */
    ClearFaceDescriptors();
    for (i = 1; i <= dom; i++)
      AddFaceDescriptor (FaceDescriptor (0, i, 0, 0));
    CalcSurfacesOfNode();
    timestamp = NextTimeStamp();
  }

  void Mesh :: SplitSeparatedFaces ()
  {
    PrintMessage (3, "SplitSeparateFaces");
    int fdi;
    int np = GetNP();

    NgBitArray usedp(np);
    Array<SurfaceElementIndex> els_of_face;

    fdi = 1;
    while (fdi <= GetNFD())
      {
        GetSurfaceElementsOfFace (fdi, els_of_face);

        if (els_of_face.Size() == 0) continue;

        SurfaceElementIndex firstel = els_of_face[0];

        usedp.Clear();
        for (int j = 1; j <= SurfaceElement(firstel).GetNP(); j++)
          usedp.Set (SurfaceElement(firstel).PNum(j));

        bool changed;
        do
          {
            changed = false;

            for (int i = 0; i < els_of_face.Size(); i++)
              {
                const Element2d & el = SurfaceElement(els_of_face[i]);

                bool has = 0;
                bool hasno = 0;
                for (int j = 0; j < el.GetNP(); j++)
                  {
                    if (usedp.Test(el[j]))
                      has = true;
                    else
                      hasno = true;
                  }

                if (has && hasno)
                  changed = true;

                if (has)
                  for (int j = 0; j < el.GetNP(); j++)
                    usedp.Set (el[j]);
              }
          }
        while (changed);

        int nface = 0;
        for (int i = 0; i < els_of_face.Size(); i++)
          {
            Element2d & el = SurfaceElement(els_of_face[i]);

            int hasno = 0;
            for (int j = 1; j <= el.GetNP(); j++)
              if (!usedp.Test(el.PNum(j)))
                hasno = 1;

            if (hasno)
              {
                if (!nface)
                  {
                    FaceDescriptor nfd = GetFaceDescriptor(fdi);
                    nface = AddFaceDescriptor (nfd);
                  }

                el.SetIndex (nface);
              }
          }

        // reconnect list
        if (nface)
          {
            facedecoding[nface-1].firstelement = -1;
            facedecoding[fdi-1].firstelement = -1;

            for (int i = 0; i < els_of_face.Size(); i++)
              {
                int ind = SurfaceElement(els_of_face[i]).GetIndex();
                SurfaceElement(els_of_face[i]).next = facedecoding[ind-1].firstelement;
                facedecoding[ind-1].firstelement = els_of_face[i];
              }

            // map the segments
            for(auto& seg : segments)
              if(!usedp.Test(seg[0]) || !usedp.Test(seg[1]))
                if(seg.si == fdi)
                  seg.si = nface;
          }

        fdi++;
      }


    /*
      fdi = 1;
      while (fdi <= GetNFD())
      {
      int firstel = 0;
      for (int i = 1; i <= GetNSE(); i++)
      if (SurfaceElement(i).GetIndex() == fdi)
      {
      firstel = i;
      break;
      }
      if (!firstel) continue;

      usedp.Clear();
      for (int j = 1; j <= SurfaceElement(firstel).GetNP(); j++)
      usedp.Set (SurfaceElement(firstel).PNum(j));

      int changed;
      do
      {
      changed = 0;
      for (int i = 1; i <= GetNSE(); i++)
      {
      const Element2d & el = SurfaceElement(i);
      if (el.GetIndex() != fdi)
      continue;

      int has = 0;
      int hasno = 0;
      for (int j = 1; j <= el.GetNP(); j++)
      {
      if (usedp.Test(el.PNum(j)))
      has = 1;
      else
      hasno = 1;
      }
      if (has && hasno)
      changed = 1;

      if (has)
      for (int j = 1; j <= el.GetNP(); j++)
      usedp.Set (el.PNum(j));
      }
      }
      while (changed);

      int nface = 0;
      for (int i = 1; i <= GetNSE(); i++)
      {
      Element2d & el = SurfaceElement(i);
      if (el.GetIndex() != fdi)
      continue;	  

      int hasno = 0;
      for (int j = 1; j <= el.GetNP(); j++)
      {
      if (!usedp.Test(el.PNum(j)))
      hasno = 1;
      }

      if (hasno)
      {
      if (!nface)
      {
      FaceDescriptor nfd = GetFaceDescriptor(fdi);
      nface = AddFaceDescriptor (nfd);
      }

      el.SetIndex (nface);
      }
      }
      fdi++;
      }
    */
  }



  void Mesh :: RebuildSurfaceElementLists ()
  {
    static Timer t("Mesh::LinkSurfaceElements"); RegionTimer reg (t);    
    
    for (int i = 0; i < facedecoding.Size(); i++)
      facedecoding[i].firstelement = -1;
    for (int i = surfelements.Size()-1; i >= 0; i--)
      {
        int ind = surfelements[i].GetIndex();
        surfelements[i].next = facedecoding[ind-1].firstelement;
        facedecoding[ind-1].firstelement = i;
      }
  }

  void Mesh :: GetSurfaceElementsOfFace (int facenr, Array<SurfaceElementIndex> & sei) const
  {
    static int timer = NgProfiler::CreateTimer ("GetSurfaceElementsOfFace");
    NgProfiler::RegionTimer reg (timer);

     /*
     sei.SetSize (0);
     for (SurfaceElementIndex i = 0; i < GetNSE(); i++)
     {
        if ( (*this)[i].GetIndex () == facenr && (*this)[i][0] >= PointIndex::BASE &&
           !(*this)[i].IsDeleted() )
        {
           sei.Append (i);
        }
     }
     */

     /* Philippose - 01/10/2009
     Commented out the following lines, and activated the originally 
     commented out lines above because of a bug which causes corruption 
     of the variable "facedecoding" when a mesh is converted to second order
     */

     //      int size1 = sei.Size();
     sei.SetSize(0);

     SurfaceElementIndex si = facedecoding[facenr-1].firstelement;
     while (si != -1)
     {
        if ( (*this)[si].GetIndex () == facenr && (*this)[si][0] >= PointIndex::BASE &&
             !(*this)[si].IsDeleted() )
        {
           sei.Append (si);
        }

        si = (*this)[si].next;
     }
     
     /*
     // *testout << "with list = " << endl << sei << endl;

     if (size1 != sei.Size()) 
     {
        cout << "size mismatch" << endl;
        exit(1);
     }
     */
  }




  void Mesh :: CalcMinMaxAngle (double badellimit, double * retvalues) 
  {
    int i, j;
    int lpi1, lpi2, lpi3, lpi4;
    double phimax = 0, phimin = 10;
    double facephimax = 0, facephimin = 10;
    int illegaltets = 0, negativetets = 0, badtets = 0;

    for (i = 1; i <= GetNE(); i++)
      {
        int badel = 0;

        Element & el = VolumeElement(i);

        if (el.GetType() != TET)
          {
            VolumeElement(i).flags.badel = 0;
            continue;
          }

        if (el.Volume(Points()) < 0)
          {
            badel = 1;
            negativetets++;
          }


        if (!LegalTet (el)) 
          {
            badel = 1;
            illegaltets++;
            (*testout) << "illegal tet: " << i << " ";
            for (j = 1; j <= el.GetNP(); j++)
              (*testout) << el.PNum(j) << " ";
            (*testout) << endl;
          }


        // angles between faces
        for (lpi1 = 1; lpi1 <= 3; lpi1++)
          for (lpi2 = lpi1+1; lpi2 <= 4; lpi2++)
            {
              lpi3 = 1;
              while (lpi3 == lpi1 || lpi3 == lpi2)
                lpi3++;
              lpi4 = 10 - lpi1 - lpi2 - lpi3;

              const Point3d & p1 = Point (el.PNum(lpi1));
              const Point3d & p2 = Point (el.PNum(lpi2));
              const Point3d & p3 = Point (el.PNum(lpi3));
              const Point3d & p4 = Point (el.PNum(lpi4));

              Vec3d n(p1, p2);
              n /= n.Length();
              Vec3d v1(p1, p3);
              Vec3d v2(p1, p4);

              v1 -= (n * v1) * n;
              v2 -= (n * v2) * n;

              double cosphi = (v1 * v2) / (v1.Length() * v2.Length());
              double phi = acos (cosphi);
              if (phi > phimax) phimax = phi;
              if (phi < phimin) phimin = phi;

              if ((180/M_PI) * phi > badellimit)
                badel = 1;
            }


        // angles in faces
        for (j = 1; j <= 4; j++)
          {
            Element2d face(TRIG);
            el.GetFace (j, face);
            for (lpi1 = 1; lpi1 <= 3; lpi1++)
              {
                lpi2 = lpi1 % 3 + 1;
                lpi3 = lpi2 % 3 + 1;

                const Point3d & p1 = Point (el.PNum(lpi1));
                const Point3d & p2 = Point (el.PNum(lpi2));
                const Point3d & p3 = Point (el.PNum(lpi3));

                Vec3d v1(p1, p2);
                Vec3d v2(p1, p3);
                double cosphi = (v1 * v2) / (v1.Length() * v2.Length());
                double phi = acos (cosphi);
                if (phi > facephimax) facephimax = phi;
                if (phi < facephimin) facephimin = phi;

                if ((180/M_PI) * phi > badellimit)
                  badel = 1;

              }
          }


        VolumeElement(i).flags.badel = badel;
        if (badel) badtets++;
      }

    if (!GetNE())
      {
        phimin = phimax = facephimin = facephimax = 0;
      }

    if (!retvalues)
      {
        PrintMessage (1, "");
        PrintMessage (1, "between planes:  phimin = ", (180/M_PI) * phimin,
                      " phimax = ", (180/M_PI) *phimax);
        PrintMessage (1, "inside planes:   phimin = ", (180/M_PI) * facephimin,
                      " phimax = ", (180/M_PI) * facephimax);
        PrintMessage (1, "");      
      }
    else
      {
        retvalues[0] = (180/M_PI) * facephimin;
        retvalues[1] = (180/M_PI) * facephimax;
        retvalues[2] = (180/M_PI) * phimin;
        retvalues[3] = (180/M_PI) * phimax;
      }
    PrintMessage (3, "negative tets: ", negativetets);
    PrintMessage (3, "illegal tets:  ", illegaltets);
    PrintMessage (3, "bad tets:      ", badtets);
  }


  int Mesh :: MarkIllegalElements ()
  {
    if(!boundaryedges)
      BuildBoundaryEdges();

    atomic<int> cnt = 0;
    ParallelForRange( Range(volelements), [&] (auto myrange)
    {
      int cnt_local = 0;
      for(auto & el : volelements.Range(myrange))
        if (!LegalTet (el))
          cnt_local++;
      cnt += cnt_local;
    });
    return cnt;
  }

  // #ifdef NONE
  //   void Mesh :: AddIdentification (int pi1, int pi2, int identnr)
  //   {
  //     INDEX_2 pair(pi1, pi2);
  //     //  pair.Sort();
  //     identifiedpoints->Set (pair, identnr);
  //     if (identnr > maxidentnr)
  //       maxidentnr = identnr;
  //     timestamp = NextTimeStamp();
  //   }

  //   int Mesh :: GetIdentification (int pi1, int pi2) const
  //   {
  //     INDEX_2 pair(pi1, pi2);
  //     if (identifiedpoints->Used (pair))
  //       return identifiedpoints->Get(pair);
  //     else
  //       return 0;
  //   }

  //   int Mesh :: GetIdentificationSym (int pi1, int pi2) const
  //   {
  //     INDEX_2 pair(pi1, pi2);
  //     if (identifiedpoints->Used (pair))
  //       return identifiedpoints->Get(pair);

  //     pair = INDEX_2 (pi2, pi1);
  //     if (identifiedpoints->Used (pair))
  //       return identifiedpoints->Get(pair);

  //     return 0;
  //   }


  //   void Mesh :: GetIdentificationMap (int identnr, NgArray<int> & identmap) const
  //   {
  //     int i, j;

  //     identmap.SetSize (GetNP());
  //     for (i = 1; i <= identmap.Size(); i++)
  //       identmap.Elem(i) = 0;

  //     for (i = 1; i <= identifiedpoints->GetNBags(); i++)
  //       for (j = 1; j <= identifiedpoints->GetBagSize(i); j++)
  // 	{
  // 	  INDEX_2 i2;
  // 	  int nr;
  // 	  identifiedpoints->GetData (i, j, i2, nr);

  // 	  if (nr == identnr)
  // 	    {
  // 	      identmap.Elem(i2.I1()) = i2.I2();
  // 	    }
  // 	}
  //   }


  //   void Mesh :: GetIdentificationPairs (int identnr, NgArray<INDEX_2> & identpairs) const
  //   {
  //     int i, j;

  //     identpairs.SetSize(0);

  //     for (i = 1; i <= identifiedpoints->GetNBags(); i++)
  //       for (j = 1; j <= identifiedpoints->GetBagSize(i); j++)
  // 	{
  // 	  INDEX_2 i2;
  // 	  int nr;
  // 	  identifiedpoints->GetData (i, j, i2, nr);

  // 	  if (identnr == 0 || nr == identnr)
  // 	    identpairs.Append (i2);
  // 	}
  //   }
  // #endif

  int Mesh::IdentifyPeriodicBoundaries(const string &s1,
                                       const string &s2,
                                       const Transformation<3> &mapping,
                                       double pointTolerance)
  {
    auto nr = ident->GetMaxNr() + 1;
    ident->SetType(nr, Identifications::PERIODIC);
    double lami[4];
    set<int> identified_points;
    if(pointTolerance < 0.)
      {
        Point3d pmin, pmax;
        GetBox(pmin, pmax);
        pointTolerance = 1e-8 * (pmax-pmin).Length();
      }
    for(const auto& se : surfelements)
      {
        if(GetBCName(se.index-1) != s1)
          continue;

        for(const auto& pi : se.PNums())
          {
            if(identified_points.find(pi) != identified_points.end())
              continue;
            auto pt = (*this)[pi];
            auto mapped_pt = mapping(pt);
            auto other_nr = GetElementOfPoint(mapped_pt, lami, true);
            int index = -1;
            if(other_nr != 0)
              {
                auto other_el = VolumeElement(other_nr);
                for(auto i : Range(other_el.PNums().Size()))
                  if((mapped_pt - (*this)[other_el.PNums()[i]]).Length() < pointTolerance)
                    {
                      index = i;
                      break;
                    }
                if(index == -1)
                  {
                    cout << "point coordinates = " << pt << endl;
                    cout << "mapped coordinates = " << mapped_pt << endl;
                    throw Exception("Did not find mapped point with nr " + ToString(pi) + ", are you sure your mesh is periodic?");
                  }
                auto other_pi = other_el.PNums()[index];
                identified_points.insert(pi);
                ident->Add(pi, other_pi, nr);
              }
            else
              {
                cout << "point coordinates = " << pt << endl;
                cout << "mapped coordinates = " << mapped_pt << endl;
                throw Exception("Mapped point with nr " + ToString(pi) + " is outside of mesh, are you sure your mesh is periodic?");
              }
          }
      }
    return nr;
  }

  void Mesh :: InitPointCurve(double red, double green, double blue) const
  {
    pointcurves_startpoint.Append(pointcurves.Size());
    pointcurves_red.Append(red);
    pointcurves_green.Append(green);
    pointcurves_blue.Append(blue);
  }
  void Mesh :: AddPointCurvePoint(const Point3d & pt) const
  {
    pointcurves.Append(pt);
  }
  int Mesh :: GetNumPointCurves(void) const
  {
    return pointcurves_startpoint.Size();
  }
  int Mesh :: GetNumPointsOfPointCurve(int curve) const
  {
    if(curve == pointcurves_startpoint.Size()-1)
      return (pointcurves.Size() - pointcurves_startpoint.Last());
    else
      return (pointcurves_startpoint[curve+1]-pointcurves_startpoint[curve]);
  }

  Point3d & Mesh :: GetPointCurvePoint(int curve, int n) const
  {
    return pointcurves[pointcurves_startpoint[curve]+n];
  }

  void Mesh :: GetPointCurveColor(int curve, double & red, double & green, double & blue) const
  {
    red = pointcurves_red[curve];
    green = pointcurves_green[curve];
    blue = pointcurves_blue[curve];
  }


  void Mesh :: ComputeNVertices ()
  {
    numvertices = 0;

    for (const Element & el : VolumeElements())
      for (PointIndex v : el.Vertices())
        if (v > numvertices) numvertices = v;
        
    for (const Element2d & el : SurfaceElements())
      for (PointIndex v : el.Vertices())
        if (v > numvertices) numvertices = v;

    numvertices += 1-PointIndex::BASE;
  }

  int Mesh :: GetNV () const
  {
    if (numvertices < 0)
      return GetNP();
    else
      return numvertices;
  }

  void Mesh :: SetNP (int np)
  {
    points.SetSize(np);
    //  ptyps.SetSize(np);

    int mlold = mlbetweennodes.Size();
    mlbetweennodes.SetSize(np);
    if (np > mlold)
      for (int i = mlold+PointIndex::BASE; 
           i < np+PointIndex::BASE; i++)
        {
          mlbetweennodes[i].I1() = PointIndex::BASE-1;
          mlbetweennodes[i].I2() = PointIndex::BASE-1;
        }

    GetIdentifications().SetMaxPointNr (np + PointIndex::BASE-1);
  }


  Table<ElementIndex, PointIndex> Mesh :: CreatePoint2ElementTable(std::optional<BitArray> points) const
  {
    if(points)
      {
        const auto & free_points = *points;
        return ngcore::CreateSortedTable<ElementIndex, PointIndex>( volelements.Range(),
               [&](auto & table, ElementIndex ei)
               {
                 const auto & el = (*this)[ei];

                 for (PointIndex pi : el.PNums())
                   if(free_points[pi])
                     table.Add (pi, ei);
               }, GetNP());
      }
    else
        return ngcore::CreateSortedTable<ElementIndex, PointIndex>( volelements.Range(),
               [&](auto & table, ElementIndex ei)
               {
                 const auto & el = (*this)[ei];

                 for (PointIndex pi : el.PNums())
                   table.Add (pi, ei);
               }, GetNP());
  }

  Table<SurfaceElementIndex, PointIndex> Mesh :: CreatePoint2SurfaceElementTable( int faceindex ) const
  {
    static Timer timer("Mesh::CreatePoint2SurfaceElementTable"); RegionTimer rt(timer);

    if(faceindex==0)
      {
        return ngcore::CreateSortedTable<SurfaceElementIndex, PointIndex>( surfelements.Range(),
               [&](auto & table, SurfaceElementIndex ei)
               {
                 for (PointIndex pi : (*this)[ei].PNums())
                   table.Add (pi, ei);
               }, GetNP());
      }

    Array<SurfaceElementIndex> face_els;
    GetSurfaceElementsOfFace(faceindex, face_els);
    return ngcore::CreateSortedTable<SurfaceElementIndex, PointIndex>( face_els.Range(),
           [&](auto & table, size_t i)
           {
             for (PointIndex pi : (*this)[face_els[i]].PNums())
               table.Add (pi, face_els[i]);
           }, GetNP());
  }

  

  /*
    void Mesh :: BuildConnectedNodes ()
    {
    if (PureTetMesh())
    {
    connectedtonode.SetSize(0);
    return;
    }


    int i, j, k;
    int np = GetNP();
    int ne = GetNE();
    TABLE<int> conto(np);
    for (i = 1; i <= ne; i++)
    {
    const Element & el = VolumeElement(i);

    if (el.GetType() == PRISM)
    {
    for (j = 1; j <= 6; j++)
    {
    int n1 = el.PNum (j);
    int n2 = el.PNum ((j+2)%6+1);
    //	    if (n1 != n2)
    {
    int found = 0;
    for (k = 1; k <= conto.EntrySize(n1); k++)
    if (conto.Get(n1, k) == n2)
    {
    found = 1;
    break;
    }
    if (!found)
    conto.Add (n1, n2);
    }
    }
    }
    else if (el.GetType() == PYRAMID)
    {
    for (j = 1; j <= 4; j++)
    {
    int n1, n2;
    switch (j)
    {
    case 1: n1 = 1; n2 = 4; break;
    case 2: n1 = 4; n2 = 1; break;
    case 3: n1 = 2; n2 = 3; break;
    case 4: n1 = 3; n2 = 2; break;
    }

    int found = 0;
    for (k = 1; k <= conto.EntrySize(n1); k++)
    if (conto.Get(n1, k) == n2)
    {
    found = 1;
    break;
    }
    if (!found)
    conto.Add (n1, n2);
    }
    }
    }

    connectedtonode.SetSize(np);
    for (i = 1; i <= np; i++)
    connectedtonode.Elem(i) = 0;

    for (i = 1; i <= np; i++)
    if (connectedtonode.Elem(i) == 0)
    {
    connectedtonode.Elem(i) = i;
    ConnectToNodeRec (i, i, conto);
    }



    }

    void Mesh :: ConnectToNodeRec (int node, int tonode, 
    const TABLE<int> & conto)
    {
    int i, n2;
    //  (*testout) << "connect " << node << " to " << tonode << endl;
    for (i = 1; i <= conto.EntrySize(node); i++)
    {
    n2 = conto.Get(node, i);
    if (!connectedtonode.Get(n2))
    {
    connectedtonode.Elem(n2) = tonode;
    ConnectToNodeRec (n2, tonode, conto);
    }
    }
    }
  */


  bool Mesh :: PureTrigMesh (int faceindex) const
  {
    // if (!faceindex) return !mparam.quad;
    
    if (!faceindex)
      {
	for (int i = 1; i <= GetNSE(); i++)
	  if (SurfaceElement(i).GetNP() != 3)
	    return false;
	return true;
      }

    for (int i = 1; i <= GetNSE(); i++)
      if (SurfaceElement(i).GetIndex() == faceindex &&
          SurfaceElement(i).GetNP() != 3)
        return false;
    return true;
  }

  bool Mesh :: PureTetMesh () const
  {
    for (ElementIndex ei = 0; ei < GetNE(); ei++)
      if (VolumeElement(ei).GetNP() != 4)
        return 0;
    return 1;
  }

  void Mesh :: UpdateTopology (NgTaskManager tm,
                               NgTracer tracer)
  {
    static Timer t("Update Topology"); RegionTimer reg(t);
    topology.Update(tm, tracer);
    (*tracer)("call update clusters", false);
    clusters->Update();
    (*tracer)("call update clusters", true);
#ifdef PARALLEL
    if (paralleltop)
      {
        paralleltop->Reset();
        paralleltop->UpdateCoarseGrid();
      }
#endif
    updateSignal.Emit();
  }

  void Mesh :: BuildCurvedElements  (const Refinement * ref, int aorder, bool arational)
  {
    GetCurvedElements().BuildCurvedElements (ref, aorder, arational);

    for (SegmentIndex seg = 0; seg < GetNSeg(); seg++)
      (*this)[seg].SetCurved (GetCurvedElements().IsSegmentCurved (seg));
    for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
      (*this)[sei].SetCurved (GetCurvedElements().IsSurfaceElementCurved (sei));
    for (ElementIndex ei = 0; ei < GetNE(); ei++)
      (*this)[ei].SetCurved (GetCurvedElements().IsElementCurved (ei));
    
    SetNextMajorTimeStamp();
  }

  void Mesh :: BuildCurvedElements (int aorder)
  {
    if (!GetGeometry())
      throw NgException ("don't have a geometry for mesh curving");
    
    GetCurvedElements().BuildCurvedElements (&GetGeometry()->GetRefinement(), aorder, false);

    for (SegmentIndex seg = 0; seg < GetNSeg(); seg++)
      (*this)[seg].SetCurved (GetCurvedElements().IsSegmentCurved (seg));
    for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
      (*this)[sei].SetCurved (GetCurvedElements().IsSurfaceElementCurved (sei));
    for (ElementIndex ei = 0; ei < GetNE(); ei++)
      (*this)[ei].SetCurved (GetCurvedElements().IsElementCurved (ei));
    
    SetNextMajorTimeStamp();
  }

  void Mesh :: SetMaterial (int domnr, const string & mat)
  {
    if (domnr > materials.Size())
      {
        int olds = materials.Size();
        materials.SetSize (domnr);
        for (int i = olds; i < domnr-1; i++)
          materials[i] = new string("default");
      }
    /*
    materials.Elem(domnr) = new char[strlen(mat)+1];
    strcpy (materials.Elem(domnr), mat);
    */
    materials.Elem(domnr) = new string(mat);
  }

  string Mesh :: defaultmat = "default";
  const string & Mesh :: GetMaterial (int domnr) const
  {
    if (domnr <= materials.Size())
      return *materials.Get(domnr);
    static string emptystring("default");
    return emptystring;
  }

  void Mesh ::SetNBCNames ( int nbcn )
  {
    if ( bcnames.Size() )
      for ( int i = 0; i < bcnames.Size(); i++)
        if ( bcnames[i] ) delete bcnames[i];
    bcnames.SetSize(nbcn);
    bcnames = 0;
  }

  void Mesh ::SetBCName ( int bcnr, const string & abcname )
  {
    if (bcnr >= bcnames.Size())
      {
        int oldsize = bcnames.Size();
        bcnames.SetSize (bcnr+1);  // keeps contents
        for (int i = oldsize; i <= bcnr; i++)
          bcnames[i] = nullptr;
      }

    if ( bcnames[bcnr] ) delete bcnames[bcnr];
    if ( abcname != "default" )
      bcnames[bcnr] = new string ( abcname );
    else
      bcnames[bcnr] = nullptr;

    for (auto & fd : facedecoding)
      if (fd.BCProperty() <= bcnames.Size())
        fd.SetBCName (bcnames[fd.BCProperty()-1]);
  }

  const string & Mesh ::GetBCName ( int bcnr ) const
  {
    static string defaultstring = "default";

    if ( !bcnames.Size() )
      return defaultstring;

    if (bcnr < 0 || bcnr >= bcnames.Size())
      throw RangeException("Illegal bc number ", bcnr, 0, bcnames.Size());

    if ( bcnames[bcnr] )
      return *bcnames[bcnr];
    else
      return defaultstring;
  }

  void Mesh :: SetNCD2Names( int ncd2n )
  {
    if (cd2names.Size())
      for(int i=0; i<cd2names.Size(); i++)
	if(cd2names[i]) delete cd2names[i];
    cd2names.SetSize(ncd2n);
    cd2names = 0;
  }

  void Mesh :: SetCD2Name ( int cd2nr, const string & abcname )
  {
    cd2nr--;
    (*testout) << "setCD2Name on edge " << cd2nr << " to " << abcname << endl;
    if (cd2nr >= cd2names.Size())
      {
	int oldsize = cd2names.Size();
	cd2names.SetSize(cd2nr+1);
	for(int i= oldsize; i<= cd2nr; i++)
	  cd2names[i] = nullptr;
      }
    //if (cd2names[cd2nr]) delete cd2names[cd2nr];
    if (abcname != "default")
      cd2names[cd2nr] = new string(abcname);
    else
      cd2names[cd2nr] = nullptr;
  }

  string Mesh :: cd2_default_name = "default";
  string Mesh :: default_bc = "default";
  const string & Mesh :: GetCD2Name (int cd2nr) const
  {
    static string defaultstring  = "default";
    if (!cd2names.Size())
      return defaultstring;

    if (cd2nr < 0 || cd2nr >= cd2names.Size())
      return defaultstring;

    if (cd2names[cd2nr])
      return *cd2names[cd2nr];
    else
      return defaultstring;
  }

  void Mesh :: SetNCD3Names( int ncd3n )
  {
    if (cd3names.Size())
      for(int i=0; i<cd3names.Size(); i++)
	if(cd3names[i]) delete cd3names[i];
    cd3names.SetSize(ncd3n);
    cd3names = 0;
  }

  void Mesh :: SetCD3Name ( int cd3nr, const string & abcname )
  {
    cd3nr--;
    (*testout) << "setCD3Name on vertex " << cd3nr << " to " << abcname << endl;
    if (cd3nr >= cd3names.Size())
      {
	int oldsize = cd3names.Size();
	cd3names.SetSize(cd3nr+1);
	for(int i= oldsize; i<= cd3nr; i++)
	  cd3names[i] = nullptr;
      }
    if (abcname != "default")
      cd3names[cd3nr] = new string(abcname);
    else
      cd3names[cd3nr] = nullptr;
  }
  
  int Mesh :: AddCD3Name (const string & aname)
  {
    for (int i = 0; i < cd3names.Size(); i++)
      if (*cd3names[i] == aname)
        return i;
    cd3names.Append (new string(aname));
    return cd3names.Size()-1;
  }
  
  string Mesh :: cd3_default_name = "default";
  const string & Mesh :: GetCD3Name (int cd3nr) const
  {
    static string defaultstring  = "default";
    if (!cd3names.Size())
      return defaultstring;

    if (cd3nr < 0 || cd3nr >= cd3names.Size())
      return defaultstring;

    if (cd3names[cd3nr])
      return *cd3names[cd3nr];
    else
      return defaultstring;
  }


  NgArray<string*> & Mesh :: GetRegionNamesCD (int codim)
  {
    switch (codim)
      {
      case 0: return materials;
      case 1: return bcnames;
      case 2: return cd2names;
      case 3: return cd3names;
      default: throw Exception("don't have regions of co-dimension "+ToString(codim));
      }
  }
  
  

  void Mesh :: SetUserData(const char * id, NgArray<int> & data)
  {
    if(userdata_int.Used(id))
      delete userdata_int[id];

    NgArray<int> * newdata = new NgArray<int>(data);

    userdata_int.Set(id,newdata);      
  }
  bool Mesh :: GetUserData(const char * id, NgArray<int> & data, int shift) const
  {
    if(userdata_int.Used(id))
      {
        if(data.Size() < (*userdata_int[id]).Size()+shift)
          data.SetSize((*userdata_int[id]).Size()+shift);
        for(int i=0; i<(*userdata_int[id]).Size(); i++)
          data[i+shift] = (*userdata_int[id])[i];
        return true;
      }
    else
      {
        data.SetSize(0);
        return false;
      }
  }
  void Mesh :: SetUserData(const char * id, NgArray<double> & data)
  {
    if(userdata_double.Used(id))
      delete userdata_double[id];

    NgArray<double> * newdata = new NgArray<double>(data);

    userdata_double.Set(id,newdata);      
  }
  bool Mesh :: GetUserData(const char * id, NgArray<double> & data, int shift) const
  {
    if(userdata_double.Used(id))
      {
        if(data.Size() < (*userdata_double[id]).Size()+shift)
          data.SetSize((*userdata_double[id]).Size()+shift);
        for(int i=0; i<(*userdata_double[id]).Size(); i++)
          data[i+shift] = (*userdata_double[id])[i];
        return true;
      }
    else
      {
        data.SetSize(0);
        return false;
      }
  }



  void Mesh :: PrintMemInfo (ostream & ost) const
  {
    ost << "Mesh Mem:" << endl;

    ost << GetNP() << " Points, of size " 
        << sizeof (Point3d) << " + " << sizeof(POINTTYPE) << " = "
        << GetNP() * (sizeof (Point3d) + sizeof(POINTTYPE)) << endl;

    ost << GetNSE() << " Surface elements, of size " 
        << sizeof (Element2d) << " = " 
        << GetNSE() * sizeof(Element2d) << endl;

    ost << GetNE() << " Volume elements, of size " 
        << sizeof (Element) << " = " 
        << GetNE() * sizeof(Element) << endl;

    // ost << "surfs on node:";
    // surfacesonnode.PrintMemInfo (cout);

    ost << "boundaryedges: ";
    if (boundaryedges)
      boundaryedges->PrintMemInfo (cout);

    ost << "surfelementht: ";
    if (surfelementht)
      surfelementht->PrintMemInfo (cout);
  }

  shared_ptr<Mesh> Mesh :: Mirror ( netgen::Point<3> p_plane, Vec<3> n_plane )
  {
    Mesh & m = *this;
    auto nm_ = make_shared<Mesh>();
    Mesh & nm = *nm_;
    nm = m;

    Point3d pmin, pmax;
    GetBox(pmin, pmax);
    auto v = pmax-pmin;
    double eps = v.Length()*1e-8;

    auto onPlane = [&] (const MeshPoint & p) -> bool
    {
      auto v = p_plane-p;
      auto l = v.Length();
      if(l<eps) return true;

      auto ret = fabs(v*n_plane)/l;
      return fabs(v*n_plane) < eps;
    };

    auto mirror = [&] (PointIndex pi) -> PointIndex
    {
      auto & p = m[pi];

      auto v = p_plane-p;
      auto l = v.Length();
      if(l<eps)
        return pi;

      if(fabs(v*n_plane)/l < eps)
        return pi;

      auto new_point = p + 2*(v*n_plane)*n_plane;
      return nm.AddPoint( new_point, p.GetLayer(), p.Type() );
    };

    Array<PointIndex, PointIndex> point_map;
    point_map.SetSize(GetNP());
    point_map = -1;

    for(auto pi : Range(points))
      point_map[pi] = mirror(pi);

    for(auto & el : VolumeElements())
    {
      auto nel = el;
      for(auto i : Range(el.GetNP()))
        nel[i] = point_map[el[i]];
      nm.AddVolumeElement(nel);
    }

    for (auto ei : Range(SurfaceElements()))
    {
      auto & el = m[ei];
      auto nel = el;
      for(auto i : Range(el.GetNP()))
        nel[i] = point_map[el[i]];

      if(!(nel==el))
        nm.AddSurfaceElement(nel);
    }

    for (auto ei : Range(LineSegments()))
    {
      auto & el = LineSegments()[ei];
      auto nel = el;
      bool is_same = true;

      for(auto i : Range(el.GetNP()))
      {
        auto pi = el[i];
        nel[i] = point_map[pi];
        if(point_map[pi]!=pi)
          is_same = false;
      }

      if(!is_same)
        nm.AddSegment(nel);
    }

    return nm_;
  }

}
