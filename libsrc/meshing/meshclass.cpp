#include <algorithm>
#include <mystdlib.h>
#include <atomic>
#include <map>
#include <regex>
#include <set>
#include "core/array.hpp"
#include "meshing.hpp"
#include "../general/gzstream.h"

#include <core/register_archive.hpp>
#include "basegeom.hpp"

namespace netgen
{
  ElementIndex Find3dElement (const Mesh& mesh,
                              const netgen::Point<3> & p,
                              double * lami,
                              optional<FlatArray<int>> indices,
                              BoxTree<3, ElementIndex> * searchtree,
                              const bool allowindex,
                              double tol=1e-4)
  {
    int ne = 0;
    Array<ElementIndex> locels;
    if (searchtree)
      {
        searchtree->GetIntersecting (p, p, locels);
        ne = locels.Size();
      }
    else
      ne = mesh.GetNE();

    for (auto i : Range(ne))
      {
        ElementIndex ei;

        if (searchtree)
          ei = locels[i];
        else
          ei = ElementIndex::FromNr0(i);

        if(indices && indices->Size() > 0)
          {
            bool contained = indices->Contains(mesh[ei].GetIndex().Nr1());
            if((allowindex && !contained) || (!allowindex && contained)) continue;
          }

        if(mesh.PointContainedIn3DElement(p,lami,ei, tol))
          return ei;
      }

    // Not found, try uncurved variant:
    for (auto i : Range(ne))
      {
        ElementIndex ei;

        if (searchtree)
          ei = locels[i];
        else
          ei = ElementIndex::FromNr0(i);

        if(indices && indices->Size() > 0)
          {
            bool contained = indices->Contains(mesh[ei].GetIndex().Nr1());
            if((allowindex && !contained) || (!allowindex && contained)) continue;
          }


        if(mesh.PointContainedIn3DElementOld(p,lami,ei, tol))
          {
            (*testout) << "WARNING: found element of point " << p <<" only for uncurved mesh" << endl;
            return ei;
          }
      }
    return ElementIndex::INVALID;
  }

  SurfaceElementIndex
  Find2dElement (const Mesh& mesh,
                 const netgen::Point<3> & p,
                 double * lami,
                 std::optional<FlatArray<int>> indices,
                 BoxTree<3, SurfaceElementIndex> * searchtree,
                 bool allowindex)
  {
    double vlam[3];
    ElementIndex velement = ElementIndex::INVALID;

    if(mesh.GetNE())
      {
        if(searchtree)
          const_cast<Mesh&>(mesh).BuildElementSearchTree(3);
        velement = Find3dElement(mesh, p,vlam, nullopt,searchtree ? mesh.GetElementSearchTree() : nullptr,allowindex);
      }

    //(*testout) << "p " << p << endl;
    //(*testout) << "velement " << velement << endl;

    // first try to find a volume element containing p and project to face
    if(velement.IsValid())
    {
      auto & topology = mesh.GetTopology();
      const auto & fnrs = topology.GetFaces(velement);
      auto faces = ArrayMem<SurfaceElementIndex,4>();
      for(auto face : fnrs)
        faces.Append(topology.GetFace2SurfaceElement(face));

      for(int i=0; i<faces.Size(); i++)
        {
          if(!faces[i].IsValid())
            continue;
          auto sel = mesh.SurfaceElement(faces[i]);
          if(indices && indices->Size() > 0 && !indices->Contains(sel.GetIndex().Nr1()))
            continue;

          auto & el = mesh[velement];
          if (el.GetType() == TET)
          {
            double lam4[4] = { vlam[0], vlam[1], vlam[2], 1.0-vlam[0]-vlam[1]-vlam[2] };
            double face_lam = lam4[i];
            if(face_lam < 1e-5)
            {
              // found volume point very close to a face -> use barycentric coordinates directly
              lami[2] = 0.0;
              for(auto j : Range(1,3))
                for(auto k : Range(4))
                  if(sel[j] == el[k])
                    lami[j-1] = lam4[k]/(1.0-face_lam);
              return SurfaceElementIndex(faces[i]);
            }
          }

          if(mesh.PointContainedIn2DElement(p,lami,faces[i],true))
            return faces[i];
        }
    }

    // Did't find any matching face of a volume element, search 2d elements directly
    int ne;

    Array<SurfaceElementIndex> locels;
    if (searchtree)
      {
        searchtree->GetIntersecting (p, p, locels);
        ne = locels.Size();
      }
    else
      ne = mesh.GetNSE();

    for (auto i : Range(ne))
      {
        SurfaceElementIndex ii;

        if (locels.Size())
          ii = locels[i];
        else
          ii = SurfaceElementIndex::FromNr0(i);

        if(indices && indices->Size() > 0)
          {
            bool contained = indices->Contains(mesh[ii].GetIndex().Nr1());
            if((allowindex && !contained) || (!allowindex && contained)) continue;
          }
        if(mesh.PointContainedIn2DElement(p,lami,ii))
          return ii;
      }
    return SurfaceElementIndex::INVALID;
  }

  SegmentIndex Find1dElement (const Mesh& mesh,
                              const netgen::Point<3> & p,
                              double * lami,
                              std::optional<FlatArray<int>> indices,
                              BoxTree<3> * searchtree,
                              const bool allowindex = true)
  {
    double vlam[3];
    if(searchtree)
      const_cast<Mesh&>(mesh).BuildElementSearchTree(2);
    auto velement = Find2dElement(mesh, p, vlam, nullopt, searchtree ? mesh.GetSurfaceElementSearchTree() : nullptr, allowindex);
    if(!velement.IsValid())
      return SegmentIndex::INVALID;

    vlam[2] = 1.-vlam[0] - vlam[1];
    // Array<int> edges;
    auto & topology = mesh.GetTopology();

    /*
    topology.GetSurfaceElementEdges(velement, edges);
    Array<SegmentIndex> segs(edges.Size());
    for(auto i : Range(edges))
      segs[i] = topology.GetSegmentOfEdge(edges[i]);
    */
    auto hedges = topology.GetEdges(velement);
    Array<SegmentIndex> segs(hedges.Size());
    for(auto i : Range(hedges))
      segs[i] = topology.GetSegmentOfEdge(hedges[i]+1);
    
    
    for(auto i : Range(segs))
      {
        if(IsInvalid(segs[i]))
          continue;
        auto& el = mesh.SurfaceElement(velement);
        if(el.GetType() == TRIG)
          {
            double seg_lam=-1;
            double lam=-1;
            auto seg = mesh.LineSegment(segs[i]);
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
                return segs[i];
              }
          }
        else
          throw NgException("Quad not implemented yet!");
      }

    return SegmentIndex::INVALID;
  }

  static mutex buildsearchtree_mutex;

  Mesh :: Mesh ()
    : topology(*this), surfarea(*this)
  {
    lochfunc = {nullptr};
    for(auto i : Range(4))
      elementsearchtreets[i] = NextTimeStamp();
    majortimestamp = timestamp = NextTimeStamp();
    hglob = 1e10;
    hmin = 0;
    numvertices = -1;
    dimension = 3;

    curvedelems = make_unique<CurvedElements> (*this);
    clusters = make_unique<AnisotropicClusters> (*this);
    ident = make_unique<Identifications> (*this);

    ps_startelement = 0;

    geomtype = NO_GEOM;

#ifdef PARALLEL
    paralleltop = make_unique<ParallelMeshTopology> (*this);
#endif
  }


  Mesh :: ~Mesh()
  {
    for(int i = 0; i < userdata_int.Size(); i++)
      delete userdata_int[i];
    for(int i = 0; i < userdata_double.Size(); i++)
      delete userdata_double[i];

    // #ifdef PARALLEL
    // delete paralleltop;
    // #endif
  }

  shared_ptr<NetgenGeometry> Mesh :: GetGeometry() const
  {
    static auto global_geometry = make_shared<NetgenGeometry>();
    return geometry ? geometry : global_geometry;
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
    regions = mesh2.regions;
    dimension = mesh2.dimension;
    hglob = mesh2.hglob;
    hmin = mesh2.hmin;
    maxhdomain = mesh2.maxhdomain;
    pointelements = mesh2.pointelements;

    numvertices = mesh2.numvertices;

    return *this;
  }


  void Mesh :: DeleteMesh()
  {
    std::lock_guard<std::mutex> lock(mutex);
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
    Regions<2>().SetSize(0);
    Regions<1>() = RegionArray<1>();

    ident = make_unique<Identifications> (*this);
    topology = MeshTopology (*this);
    curvedelems = make_unique<CurvedElements> (*this);
    clusters = make_unique<AnisotropicClusters> (*this);

    Regions<1>().SetSize(0);

#ifdef PARALLEL
    paralleltop = make_unique<ParallelMeshTopology> (*this);
#endif

    timestamp = NextTimeStamp();
  }


  void Mesh :: ClearSurfaceElements()
  { 
    surfelements.SetSize(0);
    /*
    for (int i = 0; i < Regions<2>().Size(); i++)
      Regions<2>()[i].firstelement = SurfaceElementIndex::INVALID;
    */
    for (auto & fd : Regions<2>())
      fd.firstelement = SurfaceElementIndex::INVALID;
    
    timestamp = NextTimeStamp();
  }



  PointIndex Mesh :: AddPoint (const netgen::Point<3> & p, int layer)
  { 
    return AddPoint (p, layer, INNERPOINT);
  }

  PointIndex Mesh :: AddPoint (const netgen::Point<3> & p, int layer, POINTTYPE type)
  { 

    // PointIndex pi = points.End();
    PointIndex pi = *points.Range().end();
    if (points.Size() == points.AllocSize())
      {
        std::lock_guard<std::mutex> lock(mutex);
        points.Append ( MeshPoint (p, layer, type) ); 
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
    std::lock_guard<std::mutex> lock(mutex);
    timestamp = NextTimeStamp();

    // int maxn = max2 (s[0], s[1]);
    // maxn += 1-PointIndex::BASE;
    int maxn = max2 (s[0].Nr1(),
                     s[1].Nr1());

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

    SegmentIndex si = IndexBASE<SegmentIndex>() + segments.Size();
    segments.Append (s); 
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

    
    SurfaceElementIndex si = IndexBASE<SurfaceElementIndex>() + surfelements.Size();
    if (surfelements.AllocSize() == surfelements.Size())
      {
        std::lock_guard<std::mutex> lock(mutex);
        surfelements.Append (el);
      }
    else
      {
        surfelements.Append (el);        
      }

    if (!HasFaceDescriptor(el))
      cerr << "has no face descriptor: fd.size = " << Regions<2>().Size() << ", ind = " << el.index << endl;

    surfelements.Last().next = Regions<2>()[el.index].firstelement;
    Regions<2>()[el.index].firstelement = si;

    if (SurfaceArea().Valid())
      SurfaceArea().Add (el);

    return si;
  }

  void Mesh :: SetSurfaceElement (SurfaceElementIndex sei, const Element2d & el)
  {
    /*
    int maxn = el[0];
    for (int i = 1; i < el.GetNP(); i++)
      if (el[i] > maxn) maxn = el[i];

    maxn += 1-PointIndex::BASE;
    */
    PointIndex maxpi = el[0];
    for (int i = 1; i < el.GetNP(); i++)
      if (el[i] > maxpi) maxpi = el[i];
    int maxn = maxpi.Nr1();

    
    if (maxn <= points.Size())
      {
        for (int i = 0; i < el.GetNP(); i++)
          if (points[el[i]].Type() > SURFACEPOINT)
            points[el[i]].SetType(SURFACEPOINT);
      }

    surfelements[sei] = el;
    if (!HasFaceDescriptor(el))
      cerr << "has no face descriptor: fd.size = " << Regions<2>().Size() << ", ind = " << el.index << endl;

    // add lock-free to list ... slow, call RebuildSurfaceElementLists later
    /*
    surfelements[sei].next = Regions<2>()[el.index].firstelement;
    auto & head = reinterpret_cast<atomic<SurfaceElementIndex>&> (Regions<2>()[el.index].firstelement);
    while (!head.compare_exchange_weak (surfelements[sei].next, sei))
      ;
    */

    /*
    if (SurfaceArea().Valid())
      SurfaceArea().Add (el);
    */
  }


  ElementIndex Mesh :: AddVolumeElement (const ElementRef & el)
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

    ElementIndex ve = IndexBASE<ElementIndex>() + volelements.Size();

    if (volelements.Size() == volelements.AllocSize())
      {
        std::lock_guard<std::mutex> lock(mutex);
        volelements.Append (el);
      }
    else
      {
        volelements.Append (el);
      }
    volelements.Last().Touch();
    volelements.Last().Flags().fixed = 0;
    volelements.Last().Flags().deleted = 0;

    // while (volelements.Size() > eltyps.Size())
    // eltyps.Append (FREEELEMENT);

    timestamp = NextTimeStamp();

    return ve;
  }

  void Mesh :: SetVolumeElement (ElementIndex ei, const ElementRef & el)
  {
    /*
    int maxn = el[0];
    for (int i = 1; i < el.GetNP(); i++)
      if (el[i] > maxn) maxn = el[i];

    maxn += 1-PointIndex::BASE;
    */

    if (size_t(el.GetNP()) > volelements.Width())
      volelements.SetWidth (el.GetNP());
    volelements[ei]  = el;
    volelements[ei].Touch();
    volelements[ei].Flags().fixed = 0;
    volelements[ei].Flags().deleted = 0;
  }





  void Mesh :: Save (const filesystem::path & filename) const
  {
    string ext0 = filename.stem().extension().string();
    string ext = filename.extension().string();

    if (ext0 == ".vol" && ext == ".bin")
    {
        BinaryOutArchive in(filename);
        in & const_cast<Mesh&>(*this);
        return;
    }

    ostream * outfile;
    if (ext0 == ".vol" && ext == ".gz")
      outfile = new ogzstream(filename);
    else if (ext == ".vol")
      outfile = new ofstream(filename);
    else
      outfile = new ogzstream(filesystem::path(filename).concat(".vol.gz"));

    Save(*outfile);
    delete outfile;
  }



  void Mesh :: Save (ostream & outfile) const
  {
    static Timer timer("Mesh::Save"); RegionTimer rt(timer);
    /*
    auto seg_fdi = [this](const Segment& s) -> int {
      if (HasEdgeDescriptor(s))
        { auto fdi = Regions<1>()[s.GetIndex()].GetIndex(); if (fdi.IsValid()) return fdi.Nr1(); }
      return -1;
    };
    */

    double scale = 1;  // globflags.GetNumFlag ("scale", 1);
    int inverttets = 0;  // globflags.GetDefineFlag ("inverttets");
    int invertsurf = 0;  // globflags.GetDefineFlag ("invertsurfacemesh");

    outfile << "# Generated by NETGEN " << GetLibraryVersion("netgen") << endl << endl;


    outfile << "mesh3d" << "\n";

    outfile << "dimension\n" << GetDimension() << "\n";

    outfile << "geomtype\n" << int(geomtype) << "\n";

    outfile << "\n";
    outfile << "# surfnr\tdomin\tdomout\ttlosurf\tbcprop\n";
    outfile << "facedescriptors\n";
    outfile << GetNFD() << "\n";
    for(auto & fd : FaceDescriptors())
        outfile << fd.SurfNr() << ' ' << fd.DomainIn() << ' ' << fd.DomainOut() << ' ' << fd.TLOSurface() << ' ' << fd.BCProperty() << '\n';


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

    for (auto & el : SurfaceElements())
      {
        if (el.GetIndex().IsValid())
          {
            outfile << " " << GetFaceDescriptor(el.GetIndex ()).SurfNr()+1;
            outfile << " " << GetFaceDescriptor(el.GetIndex ()).BCProperty();
            outfile << " " << GetFaceDescriptor(el.GetIndex ()).DomainIn();
            outfile << " " << GetFaceDescriptor(el.GetIndex ()).DomainOut();
          }
        else
          outfile << " 0 0 0";

        Element2d sel = el;
        if (invertsurf)
          sel.Invert();

        outfile << " " << sel.GetNP();
        for (int j = 0; j < sel.GetNP(); j++)
          outfile << " " << sel[j];

        switch (geomtype)
          {
          case GEOM_STL:
            for (int j = 1; j <= sel.GetNP(); j++)
              outfile << " " << sel.GeomInfoPi(j).trignum;
            break;
          case GEOM_OCC: case GEOM_ACIS:
            for (int j = 1; j <= sel.GetNP(); j++)
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

    for (auto el2 : VolumeElements())
      {
        outfile << el2.GetIndex();
        outfile << " " << el2.GetNP();

        Element el (el2);
        if (inverttets) el.Invert();

        for (int j = 0; j < el.GetNP(); j++)
          outfile << " " << el[j];
        outfile << "\n";
      }


    outfile << "\n" << "\n";
    //     outfile << "   surf1   surf2      p1      p2" << "\n";
    outfile << "# p1   p2   trignum1   trignum2   dist1   dist2   edsi \n";
    outfile << "edgesegmentsgi3" << "\n";
    outfile << GetNSeg() << "\n";

    for (auto & seg : LineSegments())
      {
        outfile.width(8);
        outfile << seg[0];
        outfile.width(8);
        outfile << seg[1];
        outfile << " ";
        outfile.width(8);
        outfile << seg.GeomInfo(0).trignum;
        outfile << " ";
        outfile.width(8);
        outfile << seg.GeomInfo(1).trignum;
        outfile << " ";
        outfile.width(12);
        outfile.precision(16);
        outfile << seg.EPGeomInfo(0).dist;
        outfile << " ";
        outfile.width(12);
        outfile << seg.EPGeomInfo(1).dist;
        outfile << " ";
        outfile.width(8);
        outfile << seg.GetIndex() - 1;
        outfile << "\n";
      }


    outfile << "\n" << "\n";
    outfile << "#          X             Y             Z" << "\n";
    outfile << "points" << "\n";
    outfile << GetNP() << "\n";
    outfile.precision(16);
    outfile.setf (ios::fixed, ios::floatfield);
    outfile.setf (ios::showpoint);

    /*
    for (pi = PointIndex::BASE; 
         pi < GetNP()+PointIndex::BASE; pi++)
    */
    for (PointIndex pi : (*this).Points().Range())
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

    for (int i = 0; i < pointelements.Size(); i++)
      {
        outfile.width(8);
        outfile << pointelements[i].pnum << "  ";
        outfile.width(8);
        outfile << pointelements[i].index << "\n";
      }

    if (ident -> GetMaxNr() > 0)
      {
        outfile << "identifications\n";
        Array<PointIndices<2>> identpairs;
        int cnt = 0;
        for (int i = 1; i <= ident -> GetMaxNr(); i++)
          {
            ident -> GetPairs (i, identpairs);
            cnt += identpairs.Size();
          }
        outfile << cnt << "\n";
        for (int i = 1; i <= ident -> GetMaxNr(); i++)
          {
            ident -> GetPairs (i, identpairs);
            for (auto pair : identpairs)
              {
                outfile.width (8);
                outfile << pair[0];
                outfile.width (8);
                outfile << pair[1];
                outfile.width (8);
                outfile << i << "\n";
              }
          }

        outfile << "identificationtypes\n";
        outfile << ident -> GetMaxNr() << "\n";
        for (int i = 1; i <= ident -> GetMaxNr(); i++)
          {
            int type = ident -> GetType(i);
            outfile << " " << type;
          }
        outfile << "\n";
        outfile << "identificationnames\n";
        outfile << ident -> GetMaxNr() << "\n";
        for (int i = 1; i <= ident -> GetMaxNr(); i++)
          {
            string name = ident -> GetName(i);
            if(name == "")
              name = "default";
            outfile << name << "\n";
          }
      }

    {
      auto domnames = DomainNames();
      int cntmat = 0;
      for (auto & n : domnames)
        if (n && n->length())
          cntmat++;

      if (cntmat)
        {
          outfile << "materials" << endl;
          outfile << cntmat << endl;
          for (int i = 0; i < domnames.Size(); i++)
            if (domnames[i] && domnames[i]->length())
              outfile << i+1 << " " << *domnames[i] << endl;
        }
    }


    {
      auto names = BCNamesByNumber();
      if (names.Size())
        {
          outfile << "\n\nbcnames" << endl << names.Size() << endl;
          for (int i = 0; i < names.Size(); i++)
            outfile << i+1 << "\t" << names[i] << endl;
          outfile << endl << endl;
        }
    }
    int ncd2 = dimension >= 2 ? GetNRegions(dimension-2) : 0;
    int cntcd2names = 0;
    for (int ii = 0; ii < ncd2; ii++)
      {
        auto n = GetRegionName(dimension-2, ii+1);
        if (n != "default" && !n.empty()) cntcd2names++;
      }

    if(cntcd2names)
      {
        outfile << "\n\ncd2names" << endl << ncd2 << endl;
        for (int i=0; i<ncd2; i++)
          outfile << i+1 << "\t" << GetRegionName(dimension-2, i+1) << endl;
        outfile << endl << endl;
      }

    if (Regions<1>().Size())
      {
        outfile << "\n\nedgedescriptors" << endl << Regions<1>().Size() << endl;
        for (int ii = 0; ii < Regions<1>().Size(); ii++)
          {
            const EdgeRegion & ed = Regions<1>()[EdgeRegionIndex::FromNr0(ii)];
            outfile << ed.EdgeNr() << " "
                    << ed.SurfNr(0) << " " << ed.SurfNr(1) << " "
                    << ed.SingEdgeLeft() << " " << ed.SingEdgeRight() << " "
                    << ed.TLOSurface() << " "
                    << ed.DomainIn() << " " << ed.DomainOut() << " "
                    << ed.GetName() << endl;
          }
        outfile << endl << endl;
      }

    int ncd3 = GetNCD3Names();
    int cntcd3names = 0;
    for (int ii = 0; ii < ncd3; ii++)
      if (Regions<0>()[VertexRegionIndex::FromNr0(ii)].HasName()) cntcd3names++;

    if(cntcd3names)
      {
        outfile << "\n\ncd3names" << endl << ncd3 << endl;
        for (int i=0; i<ncd3; i++)
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
            outfile << pi << "\t" << (*this)[pi].Singularity() << endl;
      }

    auto sing_left = [&](SegmentIndex si)
    { return HasEdgeDescriptor(segments[si]) ? GetEdgeDescriptor(segments[si].GetIndex()).SingEdgeLeft() : 0.0; };
    auto sing_right = [&](SegmentIndex si)
    { return HasEdgeDescriptor(segments[si]) ? GetEdgeDescriptor(segments[si].GetIndex()).SingEdgeRight() : 0.0; };

    cnt_sing = 0;
    for (SegmentIndex si : LineSegments().Range())
      if (sing_left(si)) cnt_sing++;
    if (cnt_sing)
      {
        outfile << "singular_edge_left" << endl << cnt_sing << endl;
        for (SegmentIndex si : LineSegments().Range())
          if (sing_left(si))
            outfile << si << "\t" << sing_left(si) << endl;
      }

    cnt_sing = 0;
    for (SegmentIndex si : LineSegments().Range())
      if (sing_right(si)) cnt_sing++;
    if (cnt_sing)
      {
        outfile << "singular_edge_right" << endl << cnt_sing << endl;
        for (SegmentIndex si : LineSegments().Range())
          if (sing_right(si))
            outfile << si << "\t" << sing_right(si) << endl;
      }


    cnt_sing = 0;
    for (auto & el : SurfaceElements())
      if ( GetFaceDescriptor (el.GetIndex()).domin_singular) 
        cnt_sing++;

    if (cnt_sing)
      {
        outfile << "singular_face_inside" << endl << cnt_sing << endl;
        for (SurfaceElementIndex sei : SurfaceElements().Range())
          if ( GetFaceDescriptor ((*this)[sei].GetIndex()).domin_singular) 
            outfile << sei  << "\t" << 
              GetFaceDescriptor ((*this)[sei].GetIndex()).domin_singular  << endl;
      }

    cnt_sing = 0;
    for (auto & el : SurfaceElements())
      if ( GetFaceDescriptor (el.GetIndex()).domout_singular) cnt_sing++;
    if (cnt_sing)
      {
        outfile << "singular_face_outside" << endl << cnt_sing << endl;
        for (SurfaceElementIndex sei : SurfaceElements().Range())
          if ( GetFaceDescriptor ((*this)[sei].GetIndex()).domout_singular) 
            outfile << sei << "\t" 
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

       for(int i = 1; i <= cnt_facedesc; i++)
       {
          outfile.width(8);
          outfile << GetFaceDescriptor(FaceRegionIndex::FromNr1(i)).SurfNr()+1 << " ";
          outfile.width(12);
          outfile << GetFaceDescriptor(FaceRegionIndex::FromNr1(i)).SurfColour()[0] << " ";
          outfile.width(12);
          outfile << GetFaceDescriptor(FaceRegionIndex::FromNr1(i)).SurfColour()[1] << " ";
          outfile.width(12);
          outfile << GetFaceDescriptor(FaceRegionIndex::FromNr1(i)).SurfColour()[2];
          outfile << endl;
       }

       outfile << "face_transparencies" << endl << cnt_facedesc << endl;
       for(int i = 1; i <= cnt_facedesc; i++)
         {
           outfile.width(8);
           outfile << GetFaceDescriptor(FaceRegionIndex::FromNr1(i)).SurfNr()+1 << " ";
           outfile.width(12);
           outfile << GetFaceDescriptor(FaceRegionIndex::FromNr1(i)).SurfColour()[3] << endl;
         }
    }


    if (curvedelems && curvedelems->IsHighOrder())
      {
        if (level_nv.Size() > 1 && (MeshTopology().HasParentEdges() || MeshTopology().HasParentFaces()))
          cerr << "Waring: cannot store curvedelements on refined meshes with full hierarchy" << endl;
        else
          {
            outfile << "curvedelements" << endl;
            shared_ptr<std::ostream> spoutfile(&outfile,  [](void*) noexcept {});
            TextOutArchive out(std::move(spoutfile));
            out & (*curvedelems);
          }
      }



    
    outfile << endl << endl << "endmesh" << endl << endl;
    if (geometry)
      geometry -> SaveToMeshFile (outfile);
  }



  void Mesh :: Load (const filesystem::path & filename)
  {
    PrintMessage (1, "filename = ", filename);

    string ext0 = filename.stem().extension().string();
    string ext = filename.extension().string();

    if (ext0 == ".vol" && ext == ".bin")
    {
        BinaryInArchive in(filename);
        in & (*this);
        return;
    }

    istream * infile = NULL;

    if (ext0 == ".vol" && ext == ".gz")
      infile = new igzstream (filename);
    else
      infile = new ifstream (filename);

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
        if(iline >> i)
          {
            empty_line = false;
            // skip a single whitespace character after the number, then read the rest
            if(iline.peek() == ' ' || iline.peek() == '\t')
              iline.get();
            std::getline(iline, s);
          }
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

    // int rank = GetCommunicator().Rank();
    int ntasks = GetCommunicator().Size();
    
    char str[100];
    int n;

    double scale = 1;  // globflags.GetNumFlag ("scale", 1);
    int inverttets = 0;  // globflags.GetDefineFlag ("inverttets");
    int invertsurf = 0;  // globflags.GetDefineFlag ("invertsurfacemesh");


    Regions<2>().SetSize(0);

    bool endmesh = false;

    bool has_facedescriptors = false;
    // per-segment data read alongside the segments (edgesegmentsgi2 format)
    Array<std::pair<int,int>, SegmentIndex> seg_surfnrs;
    Array<int, SegmentIndex> seg_edgenrs;
    Array<int, SegmentIndex> seg_sis;
    Array<string> bcnames2d;
    

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

        if (strcmp (str, "facedescriptors") == 0)
          {
            has_facedescriptors = true;
            int nfd;
            infile >> nfd;
            for([[maybe_unused]] auto i : Range(nfd))
            {
                int surfnr, domin, domout, tlosurf, bcprop;
                infile >> surfnr >> domin >> domout >> tlosurf >> bcprop;
                auto faceind = AddFaceDescriptor (FaceRegion(surfnr, domin, domout, tlosurf));
                GetFaceDescriptor(faceind).SetBCProperty(bcprop);
            }
          }


        if (strcmp (str, "surfaceelements") == 0 || strcmp (str, "surfaceelementsgi")==0 || strcmp (str, "surfaceelementsuv") == 0)
          {
            static Timer t1("read surface elements"); RegionTimer rt1(t1);
            infile >> n;
            PrintMessage (3, n, " surface elements");

            bool geominfo = strcmp (str, "surfaceelementsgi") == 0;
            bool uv = strcmp (str, "surfaceelementsuv") == 0;


            for (int i = 0; i < n; i++)
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
                
                for (int j = 1; j <= Regions<2>().Size(); j++)
                  if (GetFaceDescriptor(FaceRegionIndex::FromNr1(j)).SurfNr() == surfnr &&
                      GetFaceDescriptor(FaceRegionIndex::FromNr1(j)).BCProperty() == bcp &&
                      GetFaceDescriptor(FaceRegionIndex::FromNr1(j)).DomainIn() == domin &&
                      GetFaceDescriptor(FaceRegionIndex::FromNr1(j)).DomainOut() == domout)
                    faceind = j;

                // if (Regions<2>().Size()) faceind = 1;   // for timing 

                if (!faceind)
                  {
                    faceind = AddFaceDescriptor (FaceRegion(surfnr, domin, domout, 0)).Nr1();
                    GetFaceDescriptor(FaceRegionIndex::FromNr1(faceind)).SetBCProperty (bcp);
                  }

                infile >> nep;
                if (!nep) nep = 3;

                Element2d tri(nep);
                tri.SetIndex(FaceRegionIndex::FromNr1(faceind));

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
            for (int i = 0; i < n; i++)
              {
                Element el(TET);
                int hi, nep;
                infile >> hi;
                if (hi == 0) hi = 1;
                el.SetIndex(VolumeRegionIndex::FromNr1(hi));
                infile >> nep;
                el.SetNP(nep);
                el.SetCurved (nep != 4);
                for (int j = 0; j < nep; j++)
                  infile >> el[j];

                if (inverttets)
                  el.Invert();

                AddVolumeElement (el);
              }
          }


        if (strcmp (str, "edgesegments") == 0)
          {
            static Timer t1("read edge segments"); RegionTimer rt1(t1);
            infile >> n;
            for (int i = 0; i < n; i++)
              {
                Segment seg;
                int hi;
                int si_tmp;
                infile >> si_tmp >> hi >> seg[0] >> seg[1];
                seg.SetIndex(EdgeRegionIndex::FromNr1(si_tmp));
                AddSegment (seg);
              }
          }



        if (strcmp (str, "edgesegmentsgi") == 0)
          {
            static Timer t1("read edge segmentsgi"); RegionTimer rt1(t1);
            infile >> n;
            for (int i = 0; i < n; i++)
              {
                Segment seg;
                int hi;
                int si_tmp;
                infile >> si_tmp >> hi >> seg[0] >> seg[1]
                       >> seg.GeomInfo(0).trignum
                       >> seg.GeomInfo(1).trignum;
                seg.SetIndex(EdgeRegionIndex::FromNr1(si_tmp));
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

            for (int i = 0; i < n; i++)
              {
                Segment seg;
                int hi;
                int surfnr1_tmp, surfnr2_tmp;
                int edgenr_tmp;
                int si_tmp;
                int epgi_edgenr_tmp;
                infile >> si_tmp >> hi >> seg[0] >> seg[1]
                       >> seg.GeomInfo(0).trignum
                       >> seg.GeomInfo(1).trignum
                       >> surfnr1_tmp >> surfnr2_tmp
                       >> edgenr_tmp
                       >> seg.EPGeomInfo(0).dist
                       >> epgi_edgenr_tmp
                       >> seg.EPGeomInfo(1).dist;

                if (geomtype == GEOM_OCC)
                  seg.SetIndex(EdgeRegionIndex::FromNr1(edgenr_tmp));
                else if (geomtype == GEOM_CSG)
                  seg.SetIndex(EdgeRegionIndex::FromNr1(edgenr_tmp));
                else
                  seg.SetIndex(EdgeRegionIndex::FromNr1(si_tmp));

                surfnr1_tmp--;
                surfnr2_tmp--;

                seg_edgenrs.Append(edgenr_tmp);
                seg_surfnrs.Append({surfnr1_tmp, surfnr2_tmp});
                seg_sis.Append(si_tmp);
                AddSegment (seg);
              }
          }

        if (strcmp (str, "edgesegmentsgi3") == 0)
          {
            static Timer t1("read edge segmentsgi3"); RegionTimer rt1(t1);
            infile >> n;
            PrintMessage (3, n, " curve elements (gi3)");

            for (int i = 0; i < n; i++)
              {
                Segment seg;
                int edsi;
                infile >> seg[0] >> seg[1]
                       >> seg.GeomInfo(0).trignum
                       >> seg.GeomInfo(1).trignum
                       >> seg.EPGeomInfo(0).dist
                       >> seg.EPGeomInfo(1).dist
                       >> edsi;
                seg.SetIndex(EdgeRegionIndex::FromNr0(edsi));
                AddSegment (seg);
              }
          }

        if (strcmp (str, "points") == 0)
          {
            static Timer t1("read points"); RegionTimer rt1(t1);
            infile >> n;
            PrintMessage (3, n, " points");
            for (int i = 0; i < n; i++)
              {
                netgen::Point<3> p;
                infile >> p(0) >> p(1) >> p(2);
                p(0) *= scale;
                p(1) *= scale;
                p(2) *= scale;
                AddPoint (p);
              }
            PrintMessage (3, n, " points done");
          }

        if (strcmp (str, "pointelements") == 0)
          {
            static Timer t1("read point elements"); RegionTimer rt1(t1);
            infile >> n;
            PrintMessage (3, n, " pointelements");
            for (int i = 0; i < n; i++)
              {
                Element0d el;
                int index;
                infile >> el.pnum >> index;
                el.SetIndex(VertexRegionIndex::FromNr1(index));
                pointelements.Append (el);
              }
            PrintMessage (3, n, " pointelements done");
          }

        if (strcmp (str, "identifications") == 0)
          {
            infile >> n;
            PrintMessage (3, n, " identifications");
            for (int i = 0; i < n; i++)
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
            for (int i = 1; i <= n; i++)
              {
                int type;
                infile >> type;
                ident -> SetType(i,Identifications::ID_TYPE(type));
              }
          }
        if (strcmp (str, "identificationnames") == 0)
          {
            infile >> n;
            PrintMessage (3, n, " identificationnames");
            for (int i = 1; i <= n; i++)
              {
                string name;
                infile >> name;
                ident -> SetName(i,name);
              }
          }

        if (strcmp (str, "materials") == 0)
          {
            infile >> n;
            for ([[maybe_unused]] auto i : Range(n) )
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
            Array<string> names(n);
            names = "default";
            for ( [[maybe_unused]] auto i : Range(n) )
              {
                int nr;
                string nextbcname;
                ReadNumberAndName( infile, nr, nextbcname );
                if (nr >= 1 && nr <= n) names[nr-1] = nextbcname;
              }

            if ( GetDimension() == 3 )
              {
                // the file keys names by bc number
                for (auto & el : SurfaceElements())
                  if (el.GetIndex().IsValid())
                    {
                      int bcp = GetFaceDescriptor(el.GetIndex ()).BCProperty();
                      GetFaceDescriptor(el.GetIndex ()).SetBCName((bcp >= 1 && bcp <= n) ? names[bcp-1] : "default");
                    }
              }
            else if ( GetDimension() == 2 )
              bcnames2d = std::move(names);
            else
              for (auto i : Range(n))
                SetBCName(i, names[i]);
          }

        if ( strcmp (str, "cd2names" ) == 0)
          {
            infile >> n;
            Array<int> cd2nrs(n);
            for ( auto i : Range(n) )
              {
                string nextcd2name;
                ReadNumberAndName( infile, cd2nrs[i], nextcd2name );
                SetCD2NameCompat(cd2nrs[i], nextcd2name);
              }
            if (GetDimension() < 2)
              {
                throw NgException("co dim 2 elements not implemented for dimension < 2");
              }
          }

        if ( strcmp (str, "edgedescriptors" ) == 0)
          {
            infile >> n;
            Regions<1>().SetSize(n);
            for (int ii = 0; ii < n; ii++)
              {
                EdgeRegion & ed = Regions<1>()[EdgeRegionIndex::FromNr0(ii)];
                int ednr, s0, s1, tlo;
                double sl, sr;
                infile >> ednr >> s0 >> s1 >> sl >> sr >> tlo;
                ed.SetEdgeNr(ednr);
                ed.SetSurfNr(0, s0);
                ed.SetSurfNr(1, s1);
                ed.SetSingEdgeLeft(sl);
                ed.SetSingEdgeRight(sr);
                ed.SetTLOSurface(tlo);
                // try to read domin/domout (new format) or name (old format)
                int di;
                if (infile >> di)
                  {
                    ed.SetDomainIn(di);
                    int dout;
                    infile >> dout;
                    ed.SetDomainOut(dout);
                    string nm;
                    infile >> nm;
                    ed.SetName(nm);
                    // consume rest of line (may contain legacy fdindex - discard)
                    { string rest; getline(infile, rest); }
                  }
                else
                  {
                    // old format: next token is the name (not an int)
                    infile.clear();
                    string nm;
                    infile >> nm;
                    ed.SetName(nm);
                  }
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
                Regions<0>()[VertexRegionIndex::FromNr1(cd3nrs[i])].SetName(nextcd3name);
              }
            if (GetDimension() < 3)
              {
                throw NgException("co dim 3 elements not implemented for dimension < 3");
              }
          }

        if (strcmp (str, "singular_points") == 0)
          {
            infile >> n;
            for (int i = 0; i < n; i++)
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
            for (int i = 0; i < n; i++)
              {
                SegmentIndex si;
                double s; 
                infile >> si;
                infile >> s; 
                auto & seg = (*this)[si];
                if (HasEdgeDescriptor(seg))
                  GetEdgeDescriptor(seg).SetSingEdgeLeft(s);
              }
          }
        if (strcmp (str, "singular_edge_right") == 0)
          {
            infile >> n;
            for (int i = 0; i < n; i++)
              {
                SegmentIndex si;
                double s; 
                infile >> si;
                infile >> s; 
                auto & seg = (*this)[si];
                if (HasEdgeDescriptor(seg))
                  GetEdgeDescriptor(seg).SetSingEdgeRight(s);
              }
          }

        if (strcmp (str, "singular_face_inside") == 0)
          {
            infile >> n;
            for (int i = 0; i < n; i++)
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
            for (int i = 0; i < n; i++)
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
              for(int i = 1; i <= n; i++)
              {
                 int surfnr = 0;
                 Vec<4> surfcolour(0.0,1.0,0.0,1.0);

                 infile >> surfnr 
                        >> surfcolour[0]
                        >> surfcolour[1]
                        >> surfcolour[2];

                 surfnr--;

                 if(has_facedescriptors)
                 {
                    GetFaceDescriptor(FaceRegionIndex::FromNr1(i)).SetSurfColour(surfcolour);
                 }
                 else if(surfnr > 0)
                 {
                    for(int facedesc = 1; facedesc <= cnt_facedesc; facedesc++)
                    {
                       if(surfnr == GetFaceDescriptor(FaceRegionIndex::FromNr1(facedesc)).SurfNr())
                       {
                          GetFaceDescriptor(FaceRegionIndex::FromNr1(facedesc)).SetSurfColour(surfcolour);
                       }
                    }
                 }
              }
           }
        }

        if (strcmp (str, "face_transparencies") == 0)
          {
            int cnt_facedesc = GetNFD();
            infile >> n;
            // int index = 1;
            if(n == cnt_facedesc)
              {
                for(int index = 1; index <= n; index++)
                  {
                    int surfnr;
                    double transp;
                    infile >> surfnr >> transp;
                    surfnr--;
                    if(has_facedescriptors)
                    {
                       auto& fd = GetFaceDescriptor(FaceRegionIndex::FromNr1(index));
                       auto scol = fd.SurfColour();
                       scol[3] = transp;
                       fd.SetSurfColour(scol);
                    }
                    else if(surfnr > 0)
                      {
                        for(int facedesc = 1; facedesc <= cnt_facedesc; facedesc++)
                          {
                            if(surfnr == GetFaceDescriptor(FaceRegionIndex::FromNr1(facedesc)).SurfNr())
                              {
                                auto& fd = GetFaceDescriptor(FaceRegionIndex::FromNr1(facedesc));
                                auto scol = fd.SurfColour();
                                scol[3] = transp;
                                fd.SetSurfColour(scol);
                              }
                          }
                      }
                  }
              }
          }
        
        if (strcmp (str, "curvedelements") == 0)              
          {
            topology.Update();
            shared_ptr<std::istream> spinfile(&infile,  [](void*) noexcept {});
            TextInArchive in(std::move(spinfile));
            in & (*curvedelems);


            for (SurfaceElementIndex sei : SurfaceElements().Range())
              (*this)[sei].SetCurved (GetCurvedElements().IsCurved (sei));
            for (ElementIndex ei : VolumeElements().Range())
              (*this)[ei].SetCurved (GetCurvedElements().IsCurved (ei));
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

    // reconstruct edge descriptors from segment data if not loaded from file
    if (Regions<1>().Size() == 0)
      ReconstructEdgeDescriptors(&seg_surfnrs, &seg_edgenrs);
    else if (seg_edgenrs.Size() > 0)
      {
        // edgesegmentsgi2 (legacy): match each segment to its ED using
        // the redundant per-segment data that was saved alongside.
        // With per-refedge EDs, multiple EDs can share the same edgenr,
        // so we match by (edgenr, surfnr1, surfnr2) using the temp surfnr data.
        for (auto segi : segments.Range())
          {
            auto & seg = segments[segi];
            int seg_edgenr = seg_edgenrs.Range().Contains(segi) ? seg_edgenrs[segi] : -1;
            int snr1 = -1, snr2 = -1;
            if (seg_surfnrs.Range().Contains(segi))
              {
                snr1 = seg_surfnrs[segi].first;
                snr2 = seg_surfnrs[segi].second;
              }
            // Find matching ED: prefer exact (edgenr, surfnr1, surfnr2, fdindex) match
            int best = -1;
            int best_no_fdi = -1;
            for (int j = 0; j < Regions<1>().Size(); j++)
              {
                const auto & ed = Regions<1>()[EdgeRegionIndex::FromNr0(j)];
                if (ed.EdgeNr() == seg_edgenr)
                  {
                    if (ed.SurfNr(0) == snr1 && ed.SurfNr(1) == snr2)
                      {
                        int seg_si_val = seg_sis.Range().Contains(segi) ? seg_sis[segi] : -1;
                        if (ed.GetIndex().IsValid() && ed.GetIndex().Nr1() == seg_si_val)
                          { best = j; break; }  // exact match including index
                        if (best_no_fdi < 0)
                          best_no_fdi = j;  // surfnr match without fdindex
                      }
                    if (best < 0 && best_no_fdi < 0)
                      best_no_fdi = j;  // first edgenr match as fallback
                  }
              }
            if (best < 0) best = best_no_fdi;
            if (best >= 0)
              seg.SetIndex(EdgeRegionIndex::FromNr0(best));
          }
      }
    // else: edgesegmentsgi3 - segments already have correct indices

    for (auto i : Range(bcnames2d))
      SetBCName(i, bcnames2d[i]);

    RebuildFDIndices();

    SetNextMajorTimeStamp();
    //  PrintMemInfo (cout);
  }


  static const string names_in_descriptors_version = "v6.2.2607-105";

  static std::array<Array<optional<string>>, 4> ReadRegionNamesCompat (Archive & archive)
  {
    std::array<Array<string*>, 4> tmp;
    for (auto & t : tmp)
      archive & t;
    std::array<Array<optional<string>>, 4> names;
    std::set<string*> owned;
    for (int k = 0; k < 4; k++)
      {
        names[k].SetSize(tmp[k].Size());
        for (int i = 0; i < tmp[k].Size(); i++)
          {
            if (tmp[k][i]) { names[k][i] = *tmp[k][i]; owned.insert(tmp[k][i]); }
            else names[k][i] = nullopt;
          }
      }
    for (auto p : owned) delete p;
    return names;
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

        // auto rank = comm.Rank();
        
        auto & partop = GetParallelTopology();
        
        // global enumration of points:
        // not used now, but will be needed for refined meshes
        // GridFunciton pickling is not compatible, now
        // should go to paralleltopology
        
        
        
        // merge points
        Array<PointIndex, PointIndex> globnum(points.Size());
        PointIndex maxglob = PointIndex::INVALID;
        for (auto pi : Range(points))
          {
            globnum[pi] = PointIndex::FromNr1(partop.GetGlobalPNum(pi));
            // globnum[pi] = global_pnums[pi];
            maxglob = max(globnum[pi], maxglob);
          }
        
        maxglob = comm.AllReduce (maxglob, NG_MPI_MAX);
        int numglob = maxglob+1-IndexBASE<PointIndex>();
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
        for (auto el : copy_el3d)
          for (auto & pi : el.PNums())
            pi = globnum[pi];

        // strided slots are trivially copyable: send width, size and raw bytes
        if (comm.Rank() > 0)
          {
            Array<size_t> shape { copy_el3d.Width(), copy_el3d.Size() };
            comm.Send(FlatArray<size_t>(shape), 0, 200);
            comm.Send(FlatArray<char>(copy_el3d.Size()*copy_el3d.Stride(), copy_el3d.Data()), 0, 200);
          }
        else
          {
            for (int j = 1; j < comm.Size(); j++)
              {
                Array<size_t> shape(2);
                comm.Recv(FlatArray<size_t>(shape), j, 200);
                T_VOLELEMENTS el3di(shape[1], shape[0]);
                comm.Recv(FlatArray<char>(el3di.Size()*el3di.Stride(), el3di.Data()), j, 200);
                for (auto el : el3di)
                  copy_el3d.Append (el);
              }
            archive & copy_el3d;
          }


        // sending 1D elements
        auto copy_el1d  (segments);
        for (auto & el : copy_el1d)
          for (auto & pi : el.PNums())
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


        // sending 0D elements
        auto copy_el0d  (pointelements);
        for (auto & el : copy_el0d)
          {
            auto & pi = el.pnum;
            if (pi != PointIndex(PointIndex::INVALID))
              pi = globnum[pi];
          }
        
        if (comm.Rank() > 0)
          comm.Send(copy_el0d, 0, 200);
        else
          {
            Array<Element0d> el0di;
            for (int j = 1; j < comm.Size(); j++)
              {
                comm.Recv(el0di, j, 200);
                for (auto & el : el0di)
                  copy_el0d += el;
              }
            archive & copy_el0d;
          }



        
        if (comm.Rank() == 0)
          {
            archive & Regions<2>();
            archive.NeedsVersion("netgen", names_in_descriptors_version);
            ArchiveRegionNames<3>(archive);
            ArchiveRegionNames<0>(archive);
            auto mynv = numglob;
            archive & mynv;   // numvertices;
            archive & *ident;

            if(archive.GetVersion("netgen") >= "v6.2.2103-1")
              {
                archive.NeedsVersion("netgen", "v6.2.2103-1");
                archive & vol_partition & surf_partition & seg_partition;
              }
            
            archive.Shallow(geometry);
            archive & *curvedelems;

            if(archive.GetVersion("netgen") >= "v6.2.2603-26")
              {
                archive.NeedsVersion("netgen", "v6.2.2603-26");
                archive & Regions<1>();
              }
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
    archive & pointelements;
    archive & Regions<2>();
    Array<optional<string>> bcnames2d_compat;
    if (archive.GetVersion("netgen") >= names_in_descriptors_version)
      {
        archive.NeedsVersion("netgen", names_in_descriptors_version);
        ArchiveRegionNames<3>(archive);
        ArchiveRegionNames<0>(archive);
      }
    else
      {
        PrintWarning("Mesh archive written by netgen ", archive.GetVersion("netgen"),
                     " uses the old layout of region names. It is converted now, but loading it might not be",
                     " supported by future versions. Save it again to update the file.");
        auto [mats, bcnames, cd2names, cd3names] = ReadRegionNamesCompat(archive);
        SetDomainNames(std::move(mats));
        if (dimension == 1)
          SetRegionNames<0>(bcnames);   // bc names of 1D meshes are the vertex names
        if (dimension == 2)
          bcnames2d_compat = std::move(bcnames);
        if (dimension == 3)
          SetRegionNames<0>(cd3names);
        for (int i = 0; i < cd2names.Size(); i++)
          if (cd2names[i])
            SetCD2NameCompat(i+1, *cd2names[i]);
      }
    archive & numvertices;

    archive & *ident;

    // cout << "archive, ngsversion = " << archive.GetVersion("netgen") << endl;
    if(archive.GetVersion("netgen") >= "v6.2.2103-1")
      {
        // cout << "do the partition" << endl;
        archive.NeedsVersion("netgen", "v6.2.2103-1");
        archive & vol_partition & surf_partition & seg_partition;
      }
    // else
    // cout << "no partition" << endl;
    
    archive.Shallow(geometry);
    archive & *curvedelems;

    if(archive.GetVersion("netgen") >= "v6.2.2603-26")
      {
        archive.NeedsVersion("netgen", "v6.2.2603-26");
        archive & Regions<1>();
      }

    if (archive.Input())
      {
        // int rank = GetCommunicator().Rank();
        int ntasks = GetCommunicator().Size();
        
        RebuildSurfaceElementLists();
        if (Regions<1>().Size() == 0)
          ReconstructEdgeDescriptors(nullptr, nullptr);
        for (int i = 0; i < bcnames2d_compat.Size(); i++)
          if (bcnames2d_compat[i])
            SetBCName(i, *bcnames2d_compat[i]);
        RebuildFDIndices();
        
        CalcSurfacesOfNode ();
        if (ntasks == 1) // sequential run only
          {
            topology.Update();
            clusters -> Update();
          }
        SetNextMajorTimeStamp();
      }
  }


  void Mesh :: Merge (const filesystem::path & filename, const int surfindex_offset)
  {
    ifstream infile(filename);
    if (!infile.good())
      throw NgException ("mesh file not found");

    Merge(infile,surfindex_offset);

  }



  void Mesh :: Merge (istream & infile, const int surfindex_offset)
  {
    char str[100];
    int n;

    // per-segment data read alongside the merged segments, aligned with segments
    SegmentIndex first_new_seg = segments.Range().Next();
    Array<std::pair<int,int>, SegmentIndex> merge_seg_surfnrs(segments.Size());
    Array<int, SegmentIndex> merge_seg_edgenrs(segments.Size());
    Array<int, SegmentIndex> merge_seg_sis(segments.Size());
    merge_seg_surfnrs = std::pair<int,int>{-1,-1};
    merge_seg_edgenrs = -1;
    merge_seg_sis = -1;

    int inverttets = 0;  // globflags.GetDefineFlag ("inverttets");

    int oldnp = GetNP();
    int oldne = GetNSeg();
    int oldnd = GetNDomains();
    int oldned = Regions<1>().Size();
    bool merge_has_gi2 = false;

    for (auto & el : SurfaceElements())
      for(int j=1; j<=el.GetNP(); j++) el.GeomInfoPi(j).trignum = -1;

    int max_surfnr = 0;
    for (auto i : FaceDescriptors().Range())
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
            for (int i = 0; i < n; i++)
              {
                int surfnr, bcp, domin, domout, nep, faceind = 0;
                infile >> surfnr >> bcp >> domin >> domout;

                surfnr--;

                if(domin > 0) domin += oldnd;
                if(domout > 0) domout += oldnd;
                surfnr += max_surfnr;


                for (int j = 1; j <= Regions<2>().Size(); j++)
                  if (GetFaceDescriptor(FaceRegionIndex::FromNr1(j)).SurfNr() == surfnr &&
                      GetFaceDescriptor(FaceRegionIndex::FromNr1(j)).BCProperty() == bcp &&
                      GetFaceDescriptor(FaceRegionIndex::FromNr1(j)).DomainIn() == domin &&
                      GetFaceDescriptor(FaceRegionIndex::FromNr1(j)).DomainOut() == domout)
                    faceind = j;

                if (!faceind)
                  {
                    faceind = AddFaceDescriptor (FaceRegion(surfnr, domin, domout, 0)).Nr1();
                    if(GetDimension() == 2) bcp++;
                    GetFaceDescriptor(FaceRegionIndex::FromNr1(faceind)).SetBCProperty (bcp);
                  }

                infile >> nep;
                if (!nep) nep = 3;

                Element2d tri(nep);
                tri.SetIndex(FaceRegionIndex::FromNr1(faceind));

                for (int j = 1; j <= nep; j++)
                  {
                    infile >> tri.PNum(j);
                    tri.PNum(j) = tri.PNum(j) + oldnp;
                  }


                if (strcmp (str, "surfaceelementsgi") == 0)
                  for (int j = 1; j <= nep; j++)
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
            for (int i = 0; i < n; i++)
              {
                Segment seg;
                int hi;
                int si_tmp;
                infile >> si_tmp >> hi >> seg[0] >> seg[1];
                seg.SetIndex(EdgeRegionIndex::FromNr1(si_tmp));
                seg[0] = seg[0] + oldnp;
                seg[1] = seg[1] + oldnp;
                AddSegment (seg);
              }
          }



        if (strcmp (str, "edgesegmentsgi") == 0)
          {
            infile >> n;
            for (int i = 0; i < n; i++)
              {
                Segment seg;
                int hi;
                int si_tmp;
                infile >> si_tmp >> hi >> seg[0] >> seg[1]
                       >> seg.GeomInfo(0).trignum
                       >> seg.GeomInfo(1).trignum;
                seg.SetIndex(EdgeRegionIndex::FromNr1(si_tmp));
                seg[0] = seg[0] + oldnp;
                seg[1] = seg[1] + oldnp;
                AddSegment (seg);
              }
          }
        if (strcmp (str, "edgesegmentsgi2") == 0)
          {
            infile >> n;
            PrintMessage (3, n, " curve elements");

            for (int i = 0; i < n; i++)
              {
                Segment seg;
                int hi;
                int surfnr1_tmp, surfnr2_tmp;
                int edgenr_tmp;
                int si_tmp;
                int epgi_edgenr_tmp;
                infile >> si_tmp >> hi >> seg[0] >> seg[1]
                       >> seg.GeomInfo(0).trignum
                       >> seg.GeomInfo(1).trignum
                       >> surfnr1_tmp >> surfnr2_tmp
                       >> edgenr_tmp
                       >> seg.EPGeomInfo(0).dist
                       >> epgi_edgenr_tmp
                       >> seg.EPGeomInfo(1).dist;
                seg.SetIndex(EdgeRegionIndex::FromNr1(si_tmp));

                surfnr1_tmp--;
                surfnr2_tmp--;

                if(surfnr1_tmp >= 0)  surfnr1_tmp = surfnr1_tmp + max_surfnr;
                if(surfnr2_tmp >= 0)  surfnr2_tmp = surfnr2_tmp + max_surfnr;
                seg[0] = seg[0] +oldnp;
                seg[1] = seg[1] +oldnp;
                *testout << "old edgenr: " << edgenr_tmp << endl;
                edgenr_tmp = edgenr_tmp + oldne;
                *testout << "new edgenr: " << edgenr_tmp << endl;

                merge_seg_edgenrs.Append(edgenr_tmp);
                merge_seg_surfnrs.Append({surfnr1_tmp, surfnr2_tmp});
                merge_seg_sis.Append(si_tmp);
                AddSegment (seg);
                merge_has_gi2 = true;
              }
          }

        if (strcmp (str, "edgesegmentsgi3") == 0)
          {
            infile >> n;
            PrintMessage (3, n, " curve elements (gi3)");
            for (int i = 0; i < n; i++)
              {
                Segment seg;
                int edsi;
                infile >> seg[0] >> seg[1]
                       >> seg.GeomInfo(0).trignum
                       >> seg.GeomInfo(1).trignum
                       >> seg.EPGeomInfo(0).dist
                       >> seg.EPGeomInfo(1).dist
                       >> edsi;
                // index refers to the merged file's edge descriptors, appended below
                seg.SetIndex(EdgeRegionIndex::FromNr1(edsi + 1 + oldned));
                seg[0] = seg[0] + oldnp;
                seg[1] = seg[1] + oldnp;
                AddSegment (seg);
              }
          }

        if (strcmp (str, "edgedescriptors") == 0)
          {
            infile >> n;
            for (int ii = 0; ii < n; ii++)
              {
                EdgeRegion ed;
                int ednr, s0, s1, tlo;
                double sl, sr;
                infile >> ednr >> s0 >> s1 >> sl >> sr >> tlo;
                if (ednr >= 0) ednr += oldne;
                if (s0 >= 0) s0 += max_surfnr;
                if (s1 >= 0) s1 += max_surfnr;
                ed.SetEdgeNr(ednr);
                ed.SetSurfNr(0, s0);
                ed.SetSurfNr(1, s1);
                ed.SetSingEdgeLeft(sl);
                ed.SetSingEdgeRight(sr);
                ed.SetTLOSurface(tlo);
                int di;
                if (infile >> di)
                  {
                    int dout;
                    string nm;
                    infile >> dout >> nm;
                    ed.SetDomainIn(di > 0 ? di + oldnd : di);
                    ed.SetDomainOut(dout > 0 ? dout + oldnd : dout);
                    ed.SetName(nm);
                    { string rest; getline(infile, rest); }
                  }
                else
                  {
                    infile.clear();
                    string nm;
                    infile >> nm;
                    ed.SetName(nm);
                  }
                Regions<1>().Append(ed);
              }
          }

        if (strcmp (str, "volumeelements") == 0)
          {
            infile >> n;
            PrintMessage (3, n, " volume elements");
            for (int i = 0; i < n; i++)
              {
                Element el(TET);
                int hi, nep;
                infile >> hi;
                if (hi == 0) hi = 1;
                el.SetIndex(VolumeRegionIndex::FromNr1(hi+oldnd));
                infile >> nep;
                el.SetNP(nep);

                for (int j = 0; j < nep; j++)
                  {
                    infile >> el[j];
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
            for (int i = 0; i < n; i++)
              {
                netgen::Point<3> p;
                infile >> p(0) >> p(1) >> p(2);
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
            for (int i = 0; i < n; i++)
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

    if (Regions<1>().Size() == 0)
      ReconstructEdgeDescriptors(&merge_seg_surfnrs, &merge_seg_edgenrs);
    else if (merge_has_gi2)
      {
        // edgedescriptors were loaded from file; match each segment to its ED
        // With per-refedge EDs, multiple EDs can share the same edgenr,
        // so we match by (edgenr, surfnr1, surfnr2) using the temp surfnr data.
        for (auto segi : Range(first_new_seg, segments.Range().Next()))
          {
            auto & seg = segments[segi];
            int seg_edgenr = merge_seg_edgenrs.Range().Contains(segi) ? merge_seg_edgenrs[segi] : -1;
            int snr1 = -1, snr2 = -1;
            if (merge_seg_surfnrs.Range().Contains(segi))
              {
                snr1 = merge_seg_surfnrs[segi].first;
                snr2 = merge_seg_surfnrs[segi].second;
              }
            int best = -1;
            int best_no_fdi = -1;
            for (int j = 0; j < Regions<1>().Size(); j++)
              {
                const auto & ed = Regions<1>()[EdgeRegionIndex::FromNr0(j)];
                if (ed.EdgeNr() == seg_edgenr)
                  {
                    if (ed.SurfNr(0) == snr1 && ed.SurfNr(1) == snr2)
                      {
                        int seg_si_val = merge_seg_sis.Range().Contains(segi) ? merge_seg_sis[segi] : -1;
                        if (ed.GetIndex().IsValid() && ed.GetIndex().Nr1() == seg_si_val)
                          { best = j; break; }
                        if (best_no_fdi < 0)
                          best_no_fdi = j;
                      }
                    if (best < 0 && best_no_fdi < 0)
                      best_no_fdi = j;
                  }
              }
            if (best < 0) best = best_no_fdi;
            if (best >= 0)
              seg.SetIndex(EdgeRegionIndex::FromNr0(best));
          }
      }

    RebuildFDIndices();

    SetNextMajorTimeStamp();
  }










  bool Mesh :: TestOk () const
  {
    for (ElementIndex ei : volelements.Range())
      {
        for (int j = 0; j < 4; j++)
          if ( !(*this)[ei][j].IsValid())
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

    boundaryedges = make_unique<ClosedHashTable<SortedPointIndices<2>, int>>
      (3 * (GetNSE() + GetNOpenElements()) + GetNSeg() + 1);


    for (const Element2d & sel : SurfaceElements())
      {
        if (sel.IsDeleted()) continue;

        // int si = sel.GetIndex().Nr1();

        if (sel.GetNP() <= 4)
          for (int j = 0; j < sel.GetNP(); j++)
            {
              boundaryedges->Set ({ sel.PNumMod(j+1), sel.PNumMod(j+2) }, 1);
            }
        else if (sel.GetType()==TRIG6)
          {
            for (int j = 0; j < 3; j++)
              {
                boundaryedges->Set ({ sel[j], sel[(j+1)%3] }, 1);
              }
          }
        else 
          cerr << "illegal element for buildboundaryedges" << endl;
      }

    /*
    for (int i = 0; i < openelements.Size(); i++)
      {
        const Element2d & sel = openelements[i];
        for (int j = 0; j < sel.GetNP(); j++)
          {
            IVec<2> i2;
            i2[0] = sel.PNumMod(j+1);
            i2[1] = sel.PNumMod(j+2);
            i2.Sort();
            boundaryedges->Set (i2, 1);

            points[sel[j]].SetType(FIXEDPOINT);
          }
      }
    */
    for (const Element2d & sel : openelements)
      for (int j = 0; j < sel.GetNP(); j++)
        {
          boundaryedges->Set ({ sel.PNumMod(j+1), sel.PNumMod(j+2) }, 1);

          points[sel[j]].SetType(FIXEDPOINT);
        }

    /*
    for (int i = 0; i < GetNSeg(); i++)
      {
        const Segment & seg = segments[i];
        IVec<2> i2(seg[0], seg[1]);
        i2.Sort();

        boundaryedges -> Set (i2, 2);
        //segmentht -> Set (i2, i);
      }
    */
    for (const Segment & seg : segments)
      {
        boundaryedges -> Set ({ seg[0], seg[1] }, 2);
        //segmentht -> Set (i2, i);
      }

  }

  void Mesh :: ReconstructEdgeDescriptors (const Array<std::pair<int,int>, SegmentIndex> * seg_surfnrs,
                                           const Array<int, SegmentIndex> * seg_edgenrs)
  {
    Array<string> oldnames;   // names set before the reconstruction (readers name first)
    for (const auto & ed : Regions<1>()) oldnames.Append(ed.GetName());
    Regions<1>().SetSize(0);

    // find the max index value across all segments
    int maxindex = 0;
    for (auto & seg : segments)
      if (seg.GetIndex().Nr1() > maxindex)
        maxindex = seg.GetIndex().Nr1();

    if (maxindex < 1) return;

    // create edge descriptors indexed by seg.GetIndex() (1-based)
    Regions<1>().SetSize(maxindex);

    // mark which indices are used
    Array<bool> used(maxindex);
    used = false;

    for (auto segi : segments.Range())
    {
      auto & seg = segments[segi];
      int idx = seg.GetIndex().Nr1();
      if (idx < 1 || idx > maxindex) continue;

      seg.SetIndex(EdgeRegionIndex::FromNr1(idx));

      if (!used[idx-1])
      {
        used[idx-1] = true;
        int snr1 = -1, snr2 = -1;
        if (seg_surfnrs && seg_surfnrs->Range().Contains(segi))
          {
            snr1 = (*seg_surfnrs)[segi].first;
            snr2 = (*seg_surfnrs)[segi].second;
          }
        auto & ed = Regions<1>()[EdgeRegionIndex::FromNr1(idx)];
        int ednr = -1;
        if (seg_edgenrs && seg_edgenrs->Range().Contains(segi))
          ednr = (*seg_edgenrs)[segi];
        ed.SetEdgeNr(ednr);
        ed.SetSurfNr(0, snr1);
        ed.SetSurfNr(1, snr2);
        ed.SetSingEdgeLeft(0);
        ed.SetSingEdgeRight(0);
        ed.SetTLOSurface(-1);
        ed.SetDomainIn(snr1);
        ed.SetDomainOut(snr2);
      }
    }

    for (int i = 0; i < min(oldnames.Size(), Regions<1>().Size()); i++)
      if (oldnames[i] != "default")
        Regions<1>()[EdgeRegionIndex::FromNr0(i)].SetName(oldnames[i]);

    RebuildFDIndices();
  }

  void Mesh :: RebuildFDIndices ()
  {
    // Recompute EdgeRegion::index_ from surfnr + domin/domout vs face descriptors.
    for (int edi = 0; edi < Regions<1>().Size(); edi++)
      {
        auto & ed = Regions<1>()[EdgeRegionIndex::FromNr0(edi)];
        ed.SetIndex(FaceRegionIndex::INVALID);
        for (auto k : FaceDescriptors().Range())
          {
            const auto & fd = GetFaceDescriptor(k);
            if ((fd.SurfNr() == ed.SurfNr(0) || fd.SurfNr() == ed.SurfNr(1)) &&
                fd.DomainIn() == ed.DomainIn()+1 &&
                fd.DomainOut() == ed.DomainOut()+1)
              {
                ed.SetIndex(k);
                break;
              }
          }
        // fallback: match surfnr only (OCC, STL - domin/domout may be unset)
        if (!ed.GetIndex().IsValid())
          {
            for (auto k : FaceDescriptors().Range())
              {
                const auto & fd = GetFaceDescriptor(k);
                if (fd.SurfNr() == ed.SurfNr(0) || fd.SurfNr() == ed.SurfNr(1))
                  {
                    ed.SetIndex(k);
                    break;
                  }
              }
          }
      }
  }

  void Mesh :: CalcSurfacesOfNode ()
  {
    static Timer t("Mesh::CalcSurfacesOfNode"); RegionTimer reg (t);
    static Timer tn2se("Mesh::CalcSurfacesOfNode - surf on node");     
    static Timer tht("Mesh::CalcSurfacesOfNode - surfelementht"); 
    // surfacesonnode.SetSize (GetNP());
    DynamicTable<int,PointIndex> surfacesonnode(GetNP());

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
      surfelementht = make_unique<ClosedHashTable<SortedPointIndices<3>, SurfaceElementIndex>> (3*GetNSE() + 1);
    segmentht = make_unique<ClosedHashTable<SortedPointIndices<2>, SegmentIndex>> (3*GetNSeg() + 1);

    tn2se.Start();
    if (dimension == 3)
      /*
    for (SurfaceElementIndex sei : SurfaceElements().Range())
      {
        const Element2d & sel = surfelements[sei];
      */
      for (const Element2d & sel : surfelements)
        {
        if (sel.IsDeleted()) continue;

        int si = sel.GetIndex().Nr1();

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

      IVec<3> i3;
      i3[0] = sel[0];
      i3[1] = sel[1];
      i3[2] = sel[2];
      i3.Sort();
      surfelementht -> PrepareSet (i3);
      }

      surfelementht -> AllocateElements();
    */
    tn2se.Stop();
    
    tht.Start();
    if (dimension==3)
    for (SurfaceElementIndex sei : SurfaceElements().Range())
      {
        const Element2d & sel = surfelements[sei];
        if (sel.IsDeleted()) continue;

        surfelementht -> Set ({ sel[0], sel[1], sel[2] }, sei);   // war das wichtig ???    sel.GetIndex());
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
            for (SurfaceElementIndex sei : SurfaceElements().Range())
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
            for (const Element2d & sel : SurfaceElements())
              {
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

    for(const auto& pointel : pointelements)
      points[pointel.pnum].SetType(FIXEDPOINT);

    /*
      for (i = 0; i < openelements.Size(); i++)
      {
      const Element2d & sel = openelements[i];
      for (j = 0; j < sel.GetNP(); j++)
      {
      IVec<2> i2;
      i2[0] = sel.PNumMod(j+1);
      i2[1] = sel.PNumMod(j+2);
      i2.Sort();
      boundaryedges->Set (i2, 1);

      points[sel[j]].SetType(FIXEDPOINT);
      }
      }
    */

    // eltyps.SetSize (GetNE());
    // eltyps = FREEELEMENT;

    for (SegmentIndex i : segments.Range())
      {
        const Segment & seg = segments[i];
        //boundaryedges -> Set ({ seg[0], seg[1] }, 2);
        segmentht -> Set ({ seg[0], seg[1] }, i);
      }
  }

  // BitArray base is PointIndex::BASE ... 
  void Mesh :: FixPoints (const TBitArray<PointIndex> & fixpoints)
  {
    if (fixpoints.Size() != GetNP())
      {
        cerr << "Mesh::FixPoints: sizes don't fit" << endl;
        return;
      }
    /*
    int np = GetNP();
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
    // int ne = GetNE();
    int nse = GetNSE();
    
    t_table.Start();

    auto elsonpoint = ngcore::CreateSortedTable<ElementIndex, PointIndex>( volelements.Range(),
           [&](auto & table, ElementIndex ei)
           {
             auto el = (*this)[ei];
             if(el.IsDeleted()) return;
             if (dom == 0 || dom == el.GetIndex().Nr1())
               {
                 if (el.GetNP() == 4)
                   {
                     PointIndices<4> i4(el[0], el[1], el[2], el[3]);
                     i4.Sort();
                     table.Add (i4[0], ei);
                     table.Add (i4[1], ei);
                   }
                 else
                   {
                     for (PointIndex pi : el.PNums())
                       table.Add(pi, ei);
                   }
               }
           }, GetNP());


    Array<int, PointIndex> numonpoint(np);
    /*
    numonpoint = 0;
    for (ElementIndex ei = 0; ei < ne; ei++)
      {
        const Element & el = (*this)[ei];
        if (dom == 0 || dom == el.GetIndex())
          {
            if (el.GetNP() == 4)
              {
                IVec<4> i4(el[0], el[1], el[2], el[3]);
                i4.Sort();
                numonpoint[i4[0]]++;
                numonpoint[i4[1]]++;
              }
            else
              for (int j = 0; j < el.GetNP(); j++)
                numonpoint[el[j]]++;
          }
      }

    DynamicTable<ElementIndex, PointIndex> elsonpoint(np);
    for (ElementIndex ei = 0; ei < ne; ei++)
      {
        const Element & el = (*this)[ei];
        if (dom == 0 || dom == el.GetIndex())
          {
            if (el.GetNP() == 4)
              {
                IVec<4> i4(el[0], el[1], el[2], el[3]);
                i4.Sort();
                elsonpoint.Add (i4[0], ei);
                elsonpoint.Add (i4[1], ei);
              }
            else
              for (int j = 0; j < el.GetNP(); j++)
                elsonpoint.Add (el[j], ei);
          }
      }
    */
    t_table.Stop();


    Array<bool> hasface(GetNFD());

    for (int i = 1; i <= GetNFD(); i++)
      {
        int domin = GetFaceDescriptor(FaceRegionIndex::FromNr1(i)).DomainIn();
        int domout = GetFaceDescriptor(FaceRegionIndex::FromNr1(i)).DomainOut();
        hasface[i-1] = 
          ( dom == 0 && (domin != 0 || domout != 0) ) ||
          ( dom != 0 && (domin == dom || domout == dom) );
      }

    numonpoint = 0;
    for (SurfaceElementIndex sii : T_Range<SurfaceElementIndex>(nse))
      {
        auto ind = surfelements[sii].GetIndex();
        /*
          if (
          GetFaceDescriptor(ind).DomainIn() && 
          (dom == 0 || dom == GetFaceDescriptor(ind).DomainIn())
          ||
          GetFaceDescriptor(ind).DomainOut() && 
          (dom == 0 || dom == GetFaceDescriptor(ind).DomainOut())
          )
        */
        if (hasface[ind.Nr0()])
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

    DynamicTable<SurfaceElementIndex, PointIndex> selsonpoint(np);
    for (SurfaceElementIndex sii : T_Range<SurfaceElementIndex>(nse))
      {
        auto ind = surfelements[sii].GetIndex();

        /*
          if (
          GetFaceDescriptor(ind).DomainIn() && 
          (dom == 0 || dom == GetFaceDescriptor(ind).DomainIn())
          ||
          GetFaceDescriptor(ind).DomainOut() && 
          (dom == 0 || dom == GetFaceDescriptor(ind).DomainOut())
          )
        */
        if (hasface[ind.Nr0()])
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
              auto ind = hel.GetIndex();   

              if (GetFaceDescriptor(ind).DomainIn() && 
                  (dom == 0 || dom == GetFaceDescriptor(ind).DomainIn()) )
                {
                  hel.NormalizeNumbering();
                  if (hel[0] == pi)
                    {
                      IVec<3> i3(hel[0], hel[1], hel[2]);
                      tval i2;
                      i2.index = GetFaceDescriptor(ind).DomainIn();
                      i2.p4 = (hel.GetNP() == 3)
                            ? PointIndex (PointIndex::INVALID)
                      : hel[3];
                      faceht.Set (i3, i2);
                    }
                }
              if (GetFaceDescriptor(ind).DomainOut() &&
                  (dom == 0 || dom == GetFaceDescriptor(ind).DomainOut()) )
                {
                  hel.Invert();
                  hel.NormalizeNumbering();
                  if (hel[0] == pi)
                    {
                      IVec<3> i3(hel[0], hel[1], hel[2]);
                      tval i2;
                      i2.index = GetFaceDescriptor(ind).DomainOut();
                      i2.p4 = (hel.GetNP() == 3)
                        ? PointIndex (PointIndex::INVALID)
                        : hel[3];
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
                          IVec<3> i3(hel[0], hel[1], hel[2]);

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
                                        (*testout) << "p = " << (*this)[PointIndex(i3[jj-1])] << endl;
                                    }
                                }
                            }
                          else
                            {
                              hel.Invert();
                              hel.NormalizeNumbering();
                              IVec<3> i3(hel[0], hel[1], hel[2]);
                              
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
                IVec<3> i3;
                //IVec<2> i2;
                tval i2;
                faceht.GetData (i, i3, i2);
                if (i2.index != PointIndex::BASE-1)
                  {
                    Element2d tri ( (i2.p4 == PointIndex::BASE-1) ? TRIG : QUAD);
                    for (int l = 0; l < 3; l++)
                      tri[l] = i3[l];
                    tri[3] = i2.p4;
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
        // keyed on NormalizeNumbering()ed (rotated, not sorted) triples.
        // The slot walk below builds openelements, so the hash decides their order
        // and thus which mesh comes out; sized to avoid rehashing during the fill.
        ClosedHashTable<PointIndices<3>, tval> faceht(128);
        for (PointIndex pi : myrange)
          if (selsonpoint[pi].Size()+elsonpoint[pi].Size())
            {
              faceht.SetSize (4 * selsonpoint[pi].Size() + 8 * elsonpoint[pi].Size() + 16);

              for (SurfaceElementIndex sei : selsonpoint[pi])
                {
                  Element2d hel = SurfaceElement(sei);
                  if (hel.GetType() == TRIG6) hel.SetType(TRIG);
                  auto ind = hel.GetIndex();       

                  if (GetFaceDescriptor(ind).DomainIn() && 
                      (dom == 0 || dom == GetFaceDescriptor(ind).DomainIn()) )
                    {
                      hel.NormalizeNumbering();
                      if (hel[0] == pi)
                        {
                          PointIndices<3> i3(hel[0], hel[1], hel[2]);
                          tval i2;
                          i2.index = GetFaceDescriptor(ind).DomainIn();
                          i2.p4 = (hel.GetNP() == 3)
                            ? PointIndex (PointIndex::INVALID)
                            : hel[3];
                          faceht.Set (i3, i2);
                        }
                    }
                  if (GetFaceDescriptor(ind).DomainOut() &&
                      (dom == 0 || dom == GetFaceDescriptor(ind).DomainOut()) )
                    {
                      hel.Invert();
                      hel.NormalizeNumbering();
                      if (hel[0] == pi)
                        {
                          PointIndices<3> i3(hel[0], hel[1], hel[2]);
                          tval i2;
                          i2.index = GetFaceDescriptor(ind).DomainOut();
                          i2.p4 = (hel.GetNP() == 3)
                            ? PointIndex (PointIndex::INVALID)
                            : hel[3];
                          faceht.Set (i3, i2);
                        }
                    }
                }
              
              for (ElementIndex ei : elsonpoint[pi])
                {
                  auto el = VolumeElement(ei);
                  if(el.IsDeleted()) continue;
                  
                  if (dom == 0 || el.GetIndex().Nr1() == dom)
                    {
                      for (int j = 1; j <= el.GetNFaces(); j++)
                        {
                          Element2d hel(TRIG);
                          el.GetFace (j, hel);
                          hel.Invert();
                          hel.NormalizeNumbering();
                          
                          if (hel[0] == pi)
                            {
                              PointIndices<3> i3(hel[0], hel[1], hel[2]);
                              
                              if (faceht.Used (i3))
                                {
                                  tval i2 = faceht.Get(i3);
                                  if (i2.index == el.GetIndex().Nr1())
                                    {
                                      i2.index = long(PointIndex::BASE)-1;
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
                                          for (int jj = 0; jj < 3; jj++)
                                            (*testout) << "p = " << (*this)[i3[jj]] << endl;
                                        }
                                    }
                                }
                              else
                                {
                                  hel.Invert();
                                  hel.NormalizeNumbering();
                                  PointIndices<3> i3(hel[0], hel[1], hel[2]);
                                  
                                  tval i2;
                                  i2.index = el.GetIndex().Nr1();
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
                    PointIndices<3> i3;
                    tval i2;
                    faceht.GetData (i, i3, i2);
                    if (i2.index != PointIndex::BASE-1)
                      {
                        Element2d tri ( (!i2.p4.IsValid()) ? TRIG : QUAD);
                        for (int l = 0; l < 3; l++)
                          tri[l] = i3[l];
                        tri[3] = i2.p4;
                        tri.SetIndex (FaceRegionIndex::FromNr1(i2.index));
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


    string treequad;
    if (cnt4)
      treequad = " (" + ToString(cnt3) + " + " + ToString(cnt4) + ")";

    PrintMessage (5, openelements.Size(), treequad, " open elements");

    BuildBoundaryEdges();


    for (int i = 0; i < openelements.Size(); i++)
      {
        const Element2d & sel = openelements[i];

        if (boundaryedges)
          for (int j = 1; j <= sel.GetNP(); j++)
            {
              SortedPointIndices<2> i2 (sel.PNumMod(j), sel.PNumMod(j+1));
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
      IVec<2> i2(seg[0], seg[1]);
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
    auto seg_fdi = [this](const Segment& s) -> int {
      const Mesh & self = *this;
      if (self.HasEdgeDescriptor(s))
        { auto fdi = self.Regions<1>()[s.GetIndex()].GetIndex(); if (fdi.IsValid()) return fdi.Nr1(); }
      return -1;
    };
    // int i, j, k;

    // new version, general elements
    // hash index: pnum1-2, surfnr
    // hash data : surfel-nr (pos) or segment nr(neg)
    // key: the (oriented) segment point pair plus its face index
    ClosedHashTable<std::tuple<PointIndices<2>, int>, int> faceht(2*(4 * GetNSE()+GetNSeg())+8);

    PrintMessage (5, "Test Opensegments");
    for (SegmentIndex i : LineSegments().Range())
      {
        const Segment & seg = (*this)[i];

        if (surfnr == 0 || seg_fdi(seg) == surfnr)
          {
            std::tuple<PointIndices<2>, int> key { { seg[0], seg[1] }, seg_fdi(seg) };
            int data = -i.Nr1();

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
            IVec<2> key(seg[1], seg[0]);
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


    for (SurfaceElementIndex i : SurfaceElements().Range())
      {
        const Element2d & el = (*this)[i];
        if (el.IsDeleted()) continue;

        if (surfnr == 0 || el.GetIndex().Nr1() == surfnr)
          {
            for (int j = 1; j <= el.GetNP(); j++)
              {
                auto [pi1, pi2] = PointIndices<2>(el.PNumMod(j), el.PNumMod(j+1));
                std::tuple<PointIndices<2>, int> seg { { pi1, pi2 }, el.GetIndex().Nr1() };

                if (!pi1.IsValid() || !pi2.IsValid())
                  cerr << "seg = " << pi1 << "-" << pi2 << endl;

                if (faceht.Used(seg))
                  {
                    faceht.Set (seg, 0);
                    /*
                    data = faceht.Get(seg);
                    
                    if (data[0] == el.GetIndex())
                      {
                        data[0] = 0;
                        faceht.Set (seg, data);
                      }
                    else
                      {
                        // buggy = true;
                        PrintWarning ("hash table si not fitting for segment: ",
                                       seg[0], "-", seg[1], " other = ",
                                      data[1], ", surfnr = ", surfnr);
                      }
                    */
                  }
                else
                  {
                    std::get<0>(seg) = PointIndices<2>(pi2, pi1);
                    faceht.Set (seg, i.Nr1());
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
                IVec<2> i2, data;
                faceht.GetData (j, k, i2, data);
                bout << "key = " << i2 << ", data = " << data << endl;
              }
          }
        exit(1);
      }
    */

    (*testout) << "open segments: " << endl;
    opensegments.SetSize(0);
    opensegment_faces.SetSize(0);
    for (auto [key, data] : faceht)
        {
          if (data)  // surfnr
            {
              auto [i2, face] = key;
              Segment seg;
              seg[0] = i2[0];
              seg[1] = i2[1];

              // find geomdata:
              if (data > 0)
                {
                  // segment due to triangle
                  const Element2d & el = (*this)[SurfaceElementIndex::FromNr1(data)];
                  for (int k = 1; k <= el.GetNP(); k++)
                    {
                      if (seg[0] == el.PNum(k))
                        seg.GeomInfo(0) = el.GeomInfoPi(k);
                      if (seg[1] == el.PNum(k))
                        seg.GeomInfo(1) = el.GeomInfoPi(k);
                    }

                  (*testout) << "trig seg: ";
                }
              else
                {
                  // segment due to line
                  const Segment & lseg = (*this)[SegmentIndex::FromNr1(-data)];
                  seg.GeomInfo(0) = lseg.GeomInfo(0);
                  seg.GeomInfo(1) = lseg.GeomInfo(1);

                  (*testout) << "line seg: ";
                }

              (*testout) << seg[0] << " - " << seg[1] 
                         << " len = " << Dist (Point(seg[0]), Point(seg[1]))
                         << endl;

              opensegments.Append (seg);
              opensegment_faces.Append (face);
              if (seg.GeomInfo(0).trignum <= 0 || seg.GeomInfo(1).trignum <= 0)
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
    
    for (auto & seg : LineSegments())
      {
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
    IVec<2> i2;
    i2[0] = sel.PNumMod(j);
    i2[1] = sel.PNumMod(j+1);
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
    TBitArray<PointIndex> frontpoints(np);  // for 0- and 1-based
    frontpoints.Clear();
    
    for (int i = 1; i <= GetNOpenSegments(); i++)
      {
        const Segment & seg = GetOpenSegment(i);
        frontpoints.SetBit (seg[0]);
        frontpoints.SetBit (seg[1]);
      }

    for (Element2d & sel : surfelements)
      {
        bool remove = false;
        for (int j = 0; j < sel.GetNP(); j++)
          if (frontpoints.Test(sel[j]))
            remove = true;
        if (remove)
          sel[0].Invalidate();
      }

    for (int i = surfelements.Size(); i >= 1; i--)
      {
        SurfaceElementIndex sei = SurfaceElementIndex::FromNr1(i);
        if (!surfelements[sei][0].IsValid())
          {
            surfelements[sei] = surfelements.Last();
            surfelements.DeleteLast();
          }
      }

    RebuildSurfaceElementLists ();
    /*
    for (int i = 0; i < Regions<2>().Size(); i++)
      Regions<2>()[i].firstelement = SurfaceElementIndex::INVALID;
    for (int i = surfelements.Size()-1; i >= 0; i--)
      {
        int ind = surfelements[i].GetIndex();
        surfelements[i].next = Regions<2>()[FaceRegionIndex::FromNr1(ind)].firstelement;
        Regions<2>()[FaceRegionIndex::FromNr1(ind)].firstelement = i;
      }
    */

    timestamp = NextTimeStamp();
    //  Compress();
  }





  void Mesh :: FreeOpenElementsEnvironment (int layers)
  {
    static Timer timer("FreeOpenElementsEnvironment"); RegionTimer rt(timer);
    const int large = 9999;
    Array<int,PointIndex> dist(GetNP());

    dist = large;

    for (int i = 1; i <= GetNOpenElements(); i++)
      {
        const Element2d & face = OpenElement(i);
        for (int j = 0; j < face.GetNP(); j++)
          dist[face[j]] = 1;
      }

    for (int k = 1; k <= layers; k++)
      /*
      for (i = 1; i <= GetNE(); i++)
        {
          const Element & el = VolumeElement(i);
      */
      for (auto el : VolumeElements())
        {
          if (!el[0].IsValid() || el.IsDeleted()) continue;

          int elmin = large;
          for (int j = 0; j < el.GetNP(); j++)
            if (dist[el[j]] < elmin)
              elmin = dist[el[j]];
          
          if (elmin < large)
            {
              for (int j = 0; j < el.GetNP(); j++)
                if (dist[el[j]] > elmin+1)
                  dist[el[j]] = elmin+1;
            }
        }

    int cntfree = 0;
    /*
    for (int i = 1; i <= GetNE(); i++)
      {
        Element & el = VolumeElement(i);
    */
    for (auto el : VolumeElements())
      {
        if (!el[0].IsValid() || el.IsDeleted()) continue;

        int elmin = large;
        for (int j = 0; j < el.GetNP(); j++)
          if (dist[el[j]] < elmin)
            elmin = dist[el[j]];

        el.Flags().fixed = elmin > layers;
        // eltyps.Elem(i) = (elmin <= layers) ? 
        // FREEELEMENT : FIXEDELEMENT;
        if (elmin <= layers)
          cntfree++;
      }

    PrintMessage (5, "free: ", cntfree, ", fixed: ", GetNE()-cntfree);
    (*testout) << "free: " << cntfree << ", fixed: " << GetNE()-cntfree << endl;

    for (PointIndex pi = IndexBASE<PointIndex>(); 
         pi < GetNP()+IndexBASE<PointIndex>(); pi++)
      {
        if (dist[pi] > layers+1)
          points[pi].SetType(FIXEDPOINT);
      }
  }



  void Mesh :: SetLocalH (netgen::Point<3> pmin, netgen::Point<3> pmax, double grading, int layer)
  {
    using netgen::Point;
    Point<3> c = Center (pmin, pmax);
    double d = max3 (pmax(0)-pmin(0),
                     pmax(1)-pmin(1),
                     pmax(2)-pmin(2));
    d /= 2;
    Point<3> pmin2 = c - Vec<3> (d, d, d);
    Point<3> pmax2 = c + Vec<3> (d, d, d);

    SetLocalH(make_unique<LocalH> (pmin2, pmax2, grading, dimension), layer);
  }

  void Mesh :: RestrictLocalH (const netgen::Point<3> & p, double hloc, int layer)
  {
    if(hloc < hmin)
      hloc = hmin;

    //cout << "restrict h in " << p << " to " << hloc << endl;
    if (!lochfunc[layer-1])
      {
        PrintWarning("RestrictLocalH called, creating mesh-size tree");

        netgen::Point<3> boxmin, boxmax;
        GetBox (boxmin, boxmax);
        SetLocalH (boxmin, boxmax, 0.8, layer);
      }

    lochfunc[layer-1] -> SetH (p, hloc);
  }

  void Mesh :: RestrictLocalHLine (const netgen::Point<3> & p1, 
                                   const netgen::Point<3> & p2,
                                   double hloc, int layer)
  {
    if(hloc < hmin)
      hloc = hmin;

    // cout << "restrict h along " << p1 << " - " << p2 << " to " << hloc << endl;
    int steps = int (Dist (p1, p2) / hloc) + 2;
    Vec<3> v = p2 - p1;

    for (int i = 0; i <= steps; i++)
      {
        netgen::Point<3> p = p1 + (double(i)/double(steps) * v);
        RestrictLocalH (p, hloc, layer);
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
    if (dom >= 0 && dom < maxhdomain.Size())
      return maxhdomain[dom-1];
    else
      return 1e10;
  }

  void Mesh :: SetMaxHDomain (const Array<double> & mhd)
  {
    maxhdomain.SetSize(mhd.Size());
    for (int i = 0; i < mhd.Size(); i++)
      maxhdomain[i] = mhd[i];
  }


  double Mesh :: GetH (const netgen::Point<3> & p, int layer) const
  {
    const auto& lh = GetLocalH(layer);
    double hmin = hglob;
    if (lh)
      {
        double hl = lh->GetH (p);
        if (hl < hglob)
          hmin = hl;
      }
    return hmin;
  }

  double Mesh :: GetMinH (const netgen::Point<3> & pmin, const netgen::Point<3> & pmax, int layer)
  {
    const auto& lh = GetLocalH(layer);
    double hmin = hglob;
    if (lh)
      {
        double hl = lh->GetMinH (pmin, pmax);
        if (hl < hmin)
          hmin = hl;
      }
    return hmin;
  }





  double Mesh :: AverageH (int surfnr) const
  {
    int n;
    double hi, hsum;
    double maxh = 0, minh = 1e10;

    hsum = 0;
    n = 0;
    for (auto & el : SurfaceElements())
      {
        if (surfnr == 0 || el.GetIndex().Nr1() == surfnr)
          {
            for (int j = 1; j <= 3; j++)
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



  void Mesh :: CalcLocalH (double grading, int layer)
  {
    static Timer t("Mesh::CalcLocalH"); RegionTimer reg(t);
    
    if (!lochfunc[layer-1])
      {
        netgen::Point<3> pmin, pmax;
        GetBox (pmin, pmax);
        // SetLocalH (pmin, pmax, mparam.grading);
        SetLocalH (pmin, pmax, grading, layer);
      }

    PrintMessage (3,
                  "CalcLocalH: ", 
                  GetNP(), " Points ", 
                  GetNE(), " Elements ", 
                  GetNSE(), " Surface Elements");


    for (const Element2d & el : surfelements)
      {

        if (el.GetNP() == 3)
          {
            double hel = -1;
            for (int j = 1; j <= 3; j++)
              {
                const auto & p1 = points[el.PNumMod(j)];
                const auto & p2 = points[el.PNumMod(j+1)];

                /*
                  IVec<2> i21(el.PNumMod(j), el.PNumMod(j+1));
                  IVec<2> i22(el.PNumMod(j+1), el.PNumMod(j));
                  if (! identifiedpoints->Used (i21) &&
                  ! identifiedpoints->Used (i22) )
                */
                if (!ident -> UsedSymmetric (el.PNumMod(j),
                                             el.PNumMod(j+1)))
                  {
                    double hedge = Dist (p1, p2);
                    if (hedge > hel)
                      hel = hedge;
                    //            lochfunc->SetH (Center (p1, p2), 2 * Dist (p1, p2));
                    //            (*testout) << "trigseth, p1,2 = " << el.PNumMod(j) << ", " << el.PNumMod(j+1) 
                    //                       << " h = " << (2 * Dist(p1, p2)) << endl;
                  }
              }

            if (hel > 0)
              {
                const auto & p1 = points[el[0]];
                const auto & p2 = points[el[1]];
                const auto & p3 = points[el[2]];
                lochfunc[layer-1]->SetH (Center (p1, p2, p3), hel);
              }
          }
        else
          {
            {
              const auto & p1 = points[el[0]];
              const auto & p2 = points[el[1]];
              lochfunc[layer-1]->SetH (Center (p1, p2), 2 * Dist (p1, p2));
            }
            {
              const auto & p1 = points[el[2]];
              const auto & p2 = points[el[3]];
              lochfunc[layer-1]->SetH (Center (p1, p2), 2 * Dist (p1, p2));
            }
          }
      }

    for (const Segment & seg : segments)
      {
        const auto & p1 = points[seg[0]];
        const auto & p2 = points[seg[1]];
        /*
          IVec<2> i21(seg[0], seg[1]);
          IVec<2> i22(seg[1], seg[0]);
          if (identifiedpoints)
          if (!identifiedpoints->Used (i21) && !identifiedpoints->Used (i22))
        */
        if (!ident -> UsedSymmetric (seg[0], seg[1]))
          {
            lochfunc[layer-1]->SetH (Center (p1, p2), Dist (p1, p2));
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
      const auto & p1 = Point (el.PNum(j));
      const auto & p2 = Point (el.PNum(k));
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
      Point<3> pi;
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


  void Mesh :: CalcLocalHFromPointDistances(double grading, int layer)
  {
    PrintMessage (3, "Calculating local h from point distances");

    if (!lochfunc[layer-1])
      {
        netgen::Point<3> pmin, pmax;
        GetBox (pmin, pmax);

        // SetLocalH (pmin, pmax, mparam.grading);
        SetLocalH (pmin, pmax, grading, layer);
      }

    // double hl;

    for (PointIndex i = IndexBASE<PointIndex>(); 
         i < GetNP()+IndexBASE<PointIndex>(); i++)
      {
        for(PointIndex j=i+1; j<GetNP()+IndexBASE<PointIndex>(); j++)
          {
            const auto & p1 = points[i];
            const auto & p2 = points[j];
            double hl = Dist(p1,p2);
            RestrictLocalH(p1,hl);
            RestrictLocalH(p2,hl);
            //cout << "restricted h at " << p1 << " and " << p2 << " to " << hl << endl;
          }
      }


  }


  void Mesh :: CalcLocalHFromSurfaceCurvature (double grading, double elperr, int layer) 
  {
    PrintMessage (3, "Calculating local h from surface curvature");

    if (!lochfunc[layer-1])
      {
        netgen::Point<3> pmin, pmax;
        GetBox (pmin, pmax);

        // SetLocalH (pmin, pmax, mparam.grading);
        SetLocalH (pmin, pmax, grading, layer);
      }


    ClosedHashTable<SortedPointIndices<2>, int> edges(4 * GetNP() + 2);
    ClosedHashTable<SortedPointIndices<2>, int> bedges(2 * GetNSeg() + 2);

    for (auto & seg : LineSegments())
      {
        bedges.Set ({ seg[0], seg[1] }, 1);
      }
    for (SurfaceElementIndex i : SurfaceElements().Range())
      {
        const Element2d & sel = (*this)[i];
        if (!sel[0].IsValid())
          continue;
        for (int j = 1; j <= 3; j++)
          {
            SortedPointIndices<2> i2(sel.PNumMod(j), sel.PNumMod(j+1));
            if (bedges.Used(i2)) continue;

            if (edges.Used(i2))
              {
                int other = edges.Get(i2);

                const Element2d & elother = (*this)[SurfaceElementIndex::FromNr1(other)];

                int pi3_ = 1;
                while ( (sel.PNum(pi3_) == i2[0]) || 
                        (sel.PNum(pi3_) == i2[1]))
                  pi3_++;
                PointIndex pi3 = sel.PNum(pi3_);

                int pi4_ = 1;
                while ( (elother.PNum(pi4_) == i2[0]) || 
                        (elother.PNum(pi4_) == i2[1]))
                  pi4_++;
                PointIndex pi4 = elother.PNum(pi4_);

                double rad = ComputeCylinderRadius (Point (i2[0]),
                                                    Point (i2[1]),
                                                    Point (pi3),
                                                    Point (pi4));

                RestrictLocalHLine (Point(PointIndex(i2[0])), Point(PointIndex(i2[1])), rad/elperr);


                /*            
                  (*testout) << "pi1,2, 3, 4 = " << i2[0] << ", " << i2[1] << ", " << pi3 << ", " << pi4
                  << " p1 = " << Point(i2[0]) 
                  << ", p2 = " << Point(i2[1]) 
                  //                     << ", p3 = " << Point(pi3) 
                  //                     << ", p4 = " << Point(pi4) 
                  << ", rad = " << rad << endl;
                */
              }
            else
              edges.Set (i2, i.Nr1());
          }
      }


    // Restrict h due to line segments

    for (auto & seg : LineSegments())
      {
        const auto & p1 = Point(seg[0]);
        const auto & p2 = Point(seg[1]);
        RestrictLocalH (Center (p1, p2),  Dist (p1, p2));
      }



    /*


    int i, j;
    int np = GetNP();
    int nseg = GetNSeg();
    int nse = GetNSE();

    Array<Vec<3>> normals(np);
    BitArray linepoint(np);

    linepoint.Clear();
    for (i = 1; i <= nseg; i++)
    {
    linepoint.Set (LineSegment(i)[0]);
    linepoint.Set (LineSegment(i)[1]);
    }

    for (i = 1; i <= np; i++)
    normals.Elem(i) = Vec<3>(0,0,0);

    for (i = 1; i <= nse; i++)
    {
    Element2d & el = SurfaceElement(i);
    Vec<3> nf = Cross (Vec<3> (Point (el[0]), Point(el[1])),
    Vec<3> (Point (el[0]), Point(el[2])));
    for (j = 1; j <= 3; j++)
    normals.Elem(el.PNum(j)) += nf;
    }

    for (i = 1; i <= np; i++)
    normals.Elem(i) /= (1e-12 + normals.Elem(i).Length());

    for (i = 1; i <= nse; i++)
    {
    Element2d & el = SurfaceElement(i);
    Vec<3> nf = Cross (Vec<3> (Point (el[0]), Point(el[1])),
    Vec<3> (Point (el[0]), Point(el[2])));
    nf /= nf.Length();
    Point<3> c = Center (Point(el[0]),
    Point(el[1]),
    Point(el[2]));

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
    switch (rht)
      {
      case RESTRICTH_FACE:
        {
          for (const Element2d & sel : SurfaceElements())
            if (sel.GetIndex().Nr1() == nr)
              RestrictLocalH (sel, loch);
          break;
        }
      case RESTRICTH_EDGE:
        {
          for (const Segment & seg : LineSegments())
            if (GetEdgeDescriptor(seg.GetIndex()).EdgeNr() == nr)
              RestrictLocalH (seg, loch);
          break;
        }
      case RESTRICTH_POINT:
        {
          RestrictLocalH (Point (PointIndex::FromNr1(nr)), loch);
          break;
        }

      case RESTRICTH_SURFACEELEMENT:
        {
          RestrictLocalH ((*this)[SurfaceElementIndex::FromNr1(nr)], loch);
          break;
        }
      case RESTRICTH_SEGMENT:
        {
          RestrictLocalH ((*this)[SegmentIndex::FromNr1(nr)], loch);
          break;
        }
      }
  }

  void Mesh :: RestrictLocalH (const Element2d & sel, double loch)
  {
    RestrictLocalH (Center (Point(sel[0]), Point(sel[1]), Point(sel[2])), loch);
  }

  void Mesh :: RestrictLocalH (const Segment & seg, double loch)
  {
    RestrictLocalHLine (Point (seg[0]), Point(seg[1]), loch);
  }


  void Mesh :: LoadLocalMeshSize (const filesystem::path &  meshsizefilename)
  {
    // Philippose - 10/03/2009
    // Improve error checking when loading and reading
    // the local mesh size file

    if (meshsizefilename.empty()) return;

    ifstream msf(meshsizefilename);

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
        netgen::Point<3> pi;
        double hi;
        msf >> pi(0) >> pi(1) >> pi(2);
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
        netgen::Point<3> p1, p2;
        double hi;
        msf >> p1(0) >> p1(1) >> p1(2);
        msf >> p2(0) >> p2(1) >> p2(2);
        msf >> hi;
        if (!msf.good())
          throw NgException ("Mesh-size file error: Number of line definitions don't match specified list size\n");
        RestrictLocalHLine (p1, p2, hi);
      }

    msf.close();
  }



  void Mesh :: SetLocalH(shared_ptr<LocalH> loch, int layer)
  {
      if(layer>lochfunc.Size())
      {
          auto pre_size = lochfunc.Size();
          lochfunc.SetSize(layer);
          for(auto & func : lochfunc.Range(pre_size, layer-1))
              func = lochfunc[0];
      }
      lochfunc[layer-1] = loch;
  }

  void Mesh :: GetBox (netgen::Point<3> & pmin, netgen::Point<3> & pmax, int dom) const
  {
    if (points.Size() == 0)
      {
        pmin = pmax = netgen::Point<3>(0,0,0);
        return;
      }

    pmin = netgen::Point<3> (1e10, 1e10, 1e10);
    pmax = netgen::Point<3> (-1e10, -1e10, -1e10);

    auto grow = [&] (const netgen::Point<3> & p)
      {
        for (int j = 0; j < 3; j++)
          {
            pmin(j) = min2 (pmin(j), p(j));
            pmax(j) = max2 (pmax(j), p(j));
          }
      };

    if (dom <= 0)
      {
        for (PointIndex pi : points.Range())
          grow ((*this)[pi]);
      }
    else
      {
        for (auto & sel : SurfaceElements())
          {
            const Element2d & el = sel;
            if (el.IsDeleted() ) continue;

            if (dom == -1 || el.GetIndex().Nr1() == dom)
              for (int j = 0; j < 3; j++)
                grow ((*this)[el[j]]);
          }
      }

    if (pmin(0) > 0.5e10)
      pmin = pmax = netgen::Point<3>(0,0,0);
  }

  void Mesh :: GetBox (netgen::Point<3> & pmin, netgen::Point<3> & pmax, POINTTYPE ptyp) const
  {
    if (points.Size() == 0)
      {
        pmin = pmax = netgen::Point<3>(0,0,0);
        return;
      }

    pmin = netgen::Point<3> (1e10, 1e10, 1e10);
    pmax = netgen::Point<3> (-1e10, -1e10, -1e10);

    for (PointIndex pi : points.Range())
      if (points[pi].Type() <= ptyp)
        for (int j = 0; j < 3; j++)
          {
            pmin(j) = min2 (pmin(j), (*this)[pi](j));
            pmax(j) = max2 (pmax(j), (*this)[pi](j));
          }
  }




  double Mesh :: ElementError (int eli, const MeshingParameters & mp) const
  {
    auto el = volelements[ElementIndex::FromNr1(eli)];
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
    std::lock_guard<std::mutex> lock(mutex);
    
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

    // DeleteElement moves the last element into the hole, so re-check the slot
    for (auto ei = volelements.Range().First(); ei < volelements.Range().Next(); )
      if (!volelements[ei][0].IsValid() || volelements[ei].IsDeleted())
        volelements.DeleteElement(ei);
      else
        ei++;

    for (auto sei = surfelements.Range().First(); sei < surfelements.Range().Next(); )
      if (surfelements[sei].IsDeleted())
        surfelements.DeleteElement(sei);
      else
        sei++;

    for (auto si = segments.Range().First(); si < segments.Range().Next(); )
      if (!segments[si][0].IsValid() || !segments[si].GetIndex().IsValid())
        segments.DeleteElement(si);
      else
        si++;

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
         for (auto el : volelements.Range(myrange))
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
    
    for (const Segment & seg : segments)
      {
        for (int j = 0; j < seg.GetNP(); j++)
          pused[seg[j]] = true;
      }

    for(auto& pe : pointelements)
      pused[pe.pnum] = true;

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
      PointIndex npi = IndexBASE<PointIndex>();
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
         for (auto el : volelements.Range(myrange))
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

    
    for (Segment & seg : segments)
      {
        for (int j = 0; j < seg.GetNP(); j++)
          seg[j] = op2np[seg[j]];
      }

    for(auto& pe : pointelements)
      pe.pnum = op2np[pe.pnum];

    for (int i = 0; i < openelements.Size(); i++)
      {
        Element2d & el = openelements[i];
        for (int j = 0; j < el.GetNP(); j++)
          el[j] = op2np[el[j]];
      }  


    for (int i = 0; i < lockedpoints.Size(); i++)
      lockedpoints[i] = op2np[lockedpoints[i]];

    GetIdentifications().MapPoints(op2np);
    /*
    for (int i = 0; i < Regions<2>().Size(); i++)
      Regions<2>()[i].firstelement = SurfaceElementIndex::INVALID;
    for (int i = surfelements.Size()-1; i >= 0; i--)
      {
        int ind = surfelements[i].GetIndex();
        surfelements[i].next = Regions<2>()[FaceRegionIndex::FromNr1(ind)].firstelement;
        Regions<2>()[FaceRegionIndex::FromNr1(ind)].firstelement = i;
      }
    */
    RebuildSurfaceElementLists ();
    CalcSurfacesOfNode();

    topology.ClearEdges();
    topology.ClearFaces();

    //  FindOpenElements();
    timestamp = NextTimeStamp();
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

    for (auto el : volelements)
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
    ClosedHashTable<SortedPointIndices<2>, int> edges(4*nf+2);
    int err = 0;

    for (int i = 1; i <= nf; i++)
      {
        const Element2d & sel = OpenElement(i);

        for (int j = 1; j <= sel.GetNP(); j++)
          {
            PointIndices<2> e { sel.PNumMod(j), sel.PNumMod(j+1) };
            int sign = (e[1] > e[0]) ? 1 : -1;
            SortedPointIndices<2> i2 = e;
            if (!edges.Used (i2))
              edges.Set (i2, 0);
            edges.Set (i2, edges.Get(i2) + sign);
          }
      }

    for (auto [i2, cnt] : edges)
        {
          if (cnt)
            {
              PrintError ("Edge ", i2[0].Nr1() , " - ", i2[1].Nr1(), " multiple times in surface mesh");

              (*testout) << "Edge " << i2 << " multiple times in surface mesh" << endl;
              for (int k = 1; k <= nf; k++)
                {
                  const Element2d & sel = OpenElement(k);
                  for (int l = 1; l <= sel.GetNP(); l++)
                    {
                      SortedPointIndices<2> edge (sel.PNumMod(l), sel.PNumMod(l+1));

                      if (edge == i2) 
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
    
    netgen::Point<3> pmin, pmax;
    GetBox (pmin, pmax);
    BoxTree<3, SurfaceElementIndex> setree(pmin, pmax);
    // Array<SurfaceElementIndex> inters;

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
                                        ,sei.Nr0(), " and ", sej.Nr0());
                      
                          (*testout) << "Intersecting: " << endl;
                          (*testout) << "openelement " << sei << " with open element " << sej << endl;
                      
                          cout << "el1 = " << tri << endl;
                          cout << "el2 = " << tri2 << endl;
                          cout << "layer1 = " <<  (*this)[tri[0]].GetLayer() << endl;
                          cout << "layer2 = " <<  (*this)[tri2[0]].GetLayer() << endl;
                        }
                      
                      for (int k = 0; k < 3; k++)
                        (*testout) << tri[k] << "  ";
                      (*testout) << endl;
                      for (int k = 0; k < 3; k++)
                        (*testout) << tri2[k] << "  ";
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

    PrintMessage (5, "elements: ", ne);
    for (ElementIndex i : T_Range<ElementIndex>(ne))
      {
        auto el = const_cast<Mesh&>(*this)[i];
        el.Flags().badel = 0;
        int nip = el.GetNIP();
        for (int j = 1; j <= nip; j++)
          {
            el.GetTransformation (j, Points(), dtrans);
            double det = dtrans.Det();
            if (det > 0)
              {
                PrintError ("Element ", i.Nr1() , " has wrong orientation");
                el.Flags().badel = 1;
              }
          }
      }

    return 0;
  }

  // Search for surface trigs with same vertices ( may happen for instance with close surfaces in stl geometies )
  int Mesh :: FindIllegalTrigs ()
  {
    // Temporary table to store the vertex numbers of all triangles
    ClosedHashTable<SortedPointIndices<3>, SurfaceElementIndex> temp_tab(3*GetNSE() + 1);
    size_t cnt = 0;
    for (SurfaceElementIndex sei : SurfaceElements().Range())
      {
        const Element2d & sel = surfelements[sei];
        if (sel.IsDeleted()) continue;

        SortedPointIndices<3> i3(sel[0], sel[1], sel[2]);
        if(temp_tab.Used(i3))
          {
            temp_tab.Set (i3, SurfaceElementIndex::INVALID);
            cnt++;
          }
        else
          {
            temp_tab.Set (i3, sei);
          }
      }

    illegal_trigs = make_unique<ClosedHashTable<SortedPointIndices<3>, int>> (2*cnt+1);
    for (const Element2d & sel : SurfaceElements())
      {
        if (sel.IsDeleted()) continue;

        SortedPointIndices<3> i3(sel[0], sel[1], sel[2]);
        if(!temp_tab.Get(i3).IsValid())
            illegal_trigs -> Set (i3, 1);
      }
    return cnt;
  }

  bool Mesh :: LegalTrig (const Element2d & el) const
  {
      if(illegal_trigs)
      {
          if(illegal_trigs->Used({ el[0], el[1], el[2] }))
              return false;
      }

    return 1;
    // if ( /* hp */ 1)  // needed for old, simple hp-refinement
    //   { 
    //     // trigs with 2 or more segments are illegal
    //     int i;
    //     int nseg = 0;

    //     if (!segmentht)
    //       {
    //         cerr << "no segmentht allocated" << endl;
    //         return 0;
    //       }

    //     //      Point<3> cp(0.5, 0.5, 0.5);
    //     for (i = 1; i <= 3; i++)
    //       {
    //         IVec<2> i2(el.PNumMod (i), el.PNumMod (i+1));
    //         i2.Sort();
    //         if (segmentht -> Used (i2))
    //           nseg++;
    //       }
    //     if (nseg >= 2) 
    //       return 0;
    //   }
    // return 1;
  }

  double Mesh :: CalcTotalBad (const MeshingParameters & mp )
  {
    static Timer t("CalcTotalBad"); RegionTimer reg(t);
    static constexpr int n_classes = 20;

    double sum = 0;

    tets_in_qualclass.SetSize(n_classes);
    tets_in_qualclass = 0;

    ParallelForRange( volelements.Range(), [&] (auto myrange)
       {
         double local_sum = 0.0;
         double teterrpow = mp.opterrpow;

         // std::array<int,n_classes> classes_local{};
         size_t n_classes = tets_in_qualclass.Size();
         Array<int> classes_local(n_classes);
         for (int i = 0; i < n_classes; i++)
           classes_local[i] = 0;

         for (auto i : myrange)
           {
             double bad = max2(CalcBad (points, volelements[i], 0, mp),1e-10);
             double elbad;
             if (teterrpow == 1) elbad = bad;
             else if (teterrpow == 2) elbad = sqrt(bad);
             else elbad = pow(bad, 1/teterrpow);

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
  bool Mesh :: LegalTet2 (ElementRef el) const
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
        bface[i] = surfelementht->Used ({ el[gftetfacesa[i][0]],
                                         el[gftetfacesa[i][1]],
                                         el[gftetfacesa[i][2]] });
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

          size_t pos = boundaryedges -> Position(SortedPointIndices<2>(el[i], el[j]));
          if (pos != size_t(-1))
            {
              be = true;
              if (boundaryedges -> GetData(pos) == 2)
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
                // 2 boundary faces without edge in between
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

    for (int k = 0; k < Regions<2>().Size(); k++)
      {
        if (Regions<2>()[FaceRegionIndex::FromNr0(k)].DomainIn() > ndom)
          ndom = Regions<2>()[FaceRegionIndex::FromNr0(k)].DomainIn();
        if (Regions<2>()[FaceRegionIndex::FromNr0(k)].DomainOut() > ndom)
          ndom = Regions<2>()[FaceRegionIndex::FromNr0(k)].DomainOut();
      }

    return ndom;
  }

  void Mesh :: SetDimension (int dim)
  {
    // domain names of 2D/1D meshes live in the face/edge descriptors, vertex names stay
    if (dim != 3)
      Regions<3>().SetSize(0);
    dimension = dim;
  }

  void Mesh :: SurfaceMeshOrientation ()
  {
    // int i, j;
    int nse = GetNSE();

    BitArray used(nse+1);
    used.Clear();
    ClosedHashTable<PointIndices<2>, int> edges(4*nse+1);

    bool haschanged = 0;


    const Element2d & tri = (*this)[SurfaceElementIndex::FromNr1(1)];
    for (int j = 1; j <= 3; j++)
      {
        PointIndices<2> i2(tri.PNumMod(j), tri.PNumMod(j+1));
        edges.Set (i2, 1);
      }
    used.SetBit(1);

    bool unused;
    do
      {
        bool changed;
        do
          {
            changed = 0;
            for (SurfaceElementIndex i : T_Range<SurfaceElementIndex>(nse))
              if (!used.Test(i.Nr1()))
                {
                  Element2d & el = surfelements[i];
                  int found = 0, foundrev = 0;
                  for (int j = 1; j <= 3; j++)
                    {
                      PointIndices<2> i2(el.PNumMod(j), el.PNumMod(j+1));
                      if (edges.Used(i2))
                        foundrev = 1;
                      swap (i2[0], i2[1]);
                      if (edges.Used(i2))
                        found = 1;
                    }

                  if (found || foundrev)
                    {
                      if (foundrev)
                        swap (el[1], el[2]);

                      changed = 1;
                      for (int j = 1; j <= 3; j++)
                        {
                          PointIndices<2> i2(el.PNumMod(j), el.PNumMod(j+1));
                          edges.Set (i2, 1);
                        }
                      used.SetBit (i.Nr1());
                    }
                }
            if (changed)
              haschanged = 1;
          }
        while (changed);


        unused = 0;
        for (SurfaceElementIndex i : T_Range<SurfaceElementIndex>(nse))
          if (!used.Test(i.Nr1()))
            {
              unused = 1;
              const Element2d & tri = (*this)[i];
              for (int j = 1; j <= 3; j++)
                {
                  PointIndices<2> i2(tri.PNumMod(j), tri.PNumMod(j+1));
                  edges.Set (i2, 1);
                }
              used.SetBit(i.Nr1());
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
    for (ElementIndex i : T_Range<ElementIndex>(oldne))
      {
        Element el ((*this)[i]);

        if (el.GetType() == PRISM)
          {
            // prism, to 3 tets

            // make minimal node to node 1
            int minpi=0;
            PointIndex minpnum = PointIndex::FromNr0(GetNP());

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
                  swap (el.PNum(j), el[j+2]);
                minpi -= 3;
              }

            while (minpi > 1)
              {
                for (int j = 0; j <= 3; j+= 3)
                  {
                    PointIndex hi = el.PNum(1+j);
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

            if (min2 (el[1], el[5]) <
                min2 (el[2], el[4]))
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
                        (*this)[i] = nel;
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
            PointIndex minpnum = GetNP() + IndexBASE<PointIndex>();

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
                  swap (el.PNum(j), el[j+3]);
                minpi -= 4;
              }

            while (minpi > 1)
              {
                for (int j = 0; j <= 4; j+= 4)
                  {
                    PointIndex hi = el.PNum(1+j);
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
                      (*this)[i] = nel;
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
                      (*this)[i] = nel;
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
                      (*this)[i] = nel;
                    else
                      AddVolumeElement(nel);

                    firsttet = 0;
                  }
              }
            if (firsttet) cout << "no legal";
            (*testout) << endl;
          }
      }


    for (SurfaceElementIndex i : SurfaceElements().Range())
      {
        Element2d el = (*this)[i];
        if (el.GetNP() == 4)
          {
            (*testout) << "split el: " << el << " to ";

            static const int ntris[2][6] =
              { { 1, 2, 3, 1, 3, 4 },
                { 1, 2, 4, 4, 2, 3 }};

            const int * min2pi;

            if (min2 (el[0], el[2]) <
                min2 (el[1], el[3]))
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
                        (*this)[i] = nel;
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
        for (auto el : VolumeElements())
          {
            const auto & p1 = Point (el[0]);
            const auto & p2 = Point (el[1]);
            const auto & p3 = Point (el[2]);
            const auto & p4 = Point (el[3]);

            double vol = (Vec<3> (p1, p2) * 
                          Cross (Vec<3> (p1, p3), Vec<3>(p1, p4)));
            if (vol > 0)
              swap (el[2], el[3]);
          }



        UpdateTopology();
        timestamp = NextTimeStamp();
      }

    RebuildSurfaceElementLists();
  }

  void Mesh :: BuildElementSearchTree (int dim)
  {
    if(dim < 2)
      return;
    if (elementsearchtreets[dim] == GetTimeStamp())
      return;

    {
      std::lock_guard<std::mutex> guard(buildsearchtree_mutex);
      // check again to see if some other thread built while waiting for lock
      if (elementsearchtreets[dim] == GetTimeStamp()) return;

      PrintMessage (4, "Rebuild element searchtree dim " + ToString(dim));
          

      netgen::Point<3> pmin, pmax;
      GetBox(pmin, pmax);
      Box<3> box(pmin, pmax);
      box.Scale(1.2);
      if (dim == 3)
        elementsearchtree_vol = make_unique<BoxTree<3, ElementIndex>>(box);
      else
        elementsearchtree_surf = make_unique<BoxTree<3, SurfaceElementIndex>>(box);

      if (dim == 3)
        {
          for(auto ei : volelements.Range())
            {
              const auto& el = volelements[ei];
              Box<3> box (Box<3>::EMPTY_BOX);
              for (auto pi : el.PNums())
                box.Add (points[pi]);

              if(el.IsCurved() && curvedelems->IsCurved(ei))
                {
                  // add edge/face midpoints to box
                  auto eltype = el.GetType();
                  const auto verts = topology.GetVertices(eltype);

                  const auto edges = FlatArray<const ELEMENT_EDGE>(topology.GetNEdges(eltype), topology.GetEdges0(eltype));
                  for (const auto & edge: edges) {
                    netgen::Point<3> lam = netgen::Point<3>(0.5* (Vec<3>(verts[edge[0]]) + Vec<3>(verts[edge[1]])));
                    auto p = netgen::Point<3>(0.0);
                    curvedelems->CalcElementTransformation(lam,ei,p);
                    box.Add(p);
                  }

                  const auto faces = FlatArray<const ELEMENT_FACE>(topology.GetNFaces(eltype), topology.GetFaces0(eltype));
                  for (const auto & face: faces) {
                    netgen::Vec<3> lam = Vec<3>(verts[face[0]]) + Vec<3>(verts[face[1]]) + Vec<3>(verts[face[2]]);
                    if(face[3] != -1) {
                      lam += netgen::Vec<3>(verts[face[3]]);
                      lam *= 0.25;
                    }
                    else
                      lam *= 1.0/3;
                    auto p = netgen::Point<3>(0.0);
                    curvedelems->CalcElementTransformation(netgen::Point<3>(lam),ei,p);
                    box.Add(p);
                  }
                }
              box.Scale(1.2);
              elementsearchtree_vol -> Insert (box, ei);
            }
        }
      else if (dim == 2)
        {
          for (auto ei : Range(surfelements))
            {
              const auto& el = surfelements[ei];
              Box<3> box (Box<3>::EMPTY_BOX);
              for (auto pi : el.PNums())
                box.Add (points[pi]);

              if(el.IsCurved() && curvedelems->IsCurved(ei))
                {
                  netgen::Point<2>  lami [4] = {netgen::Point<2>(0.5,0), netgen::Point<2>(0,0.5), netgen::Point<2>(0.5,0.5), netgen::Point<2>(1./3,1./3)};
                  for (auto lam : lami)
                    {
                      netgen::Point<3> x;
                      Mat<3,2> Jac;

                      curvedelems->CalcSurfaceTransformation(lam,ei,x,Jac);
                      box.Add (x);
                    }
                  box.Scale(1.2);
                }
              elementsearchtree_surf -> Insert (box, ei);
            }
        }
      elementsearchtreets[dim] = GetTimeStamp();
    }
  }

  
  int SolveLinearSystemLS (const Vec<3> & col1,
                           const Vec<3> & col2,
                           const Vec<3> & rhs,
                           Vec<2> & sol)
  {
    double a11 = col1 * col1;
    double a12 = col1 * col2;
    double a22 = col2 * col2;
    
    double det = a11 * a22 - a12 * a12;
    
    if (det*det <= 1e-24 * a11 * a22)
      {
        sol = Vec<2> (0, 0);
        return 1;
      }
    
    Vec<2> aTrhs;
    aTrhs(0) = col1*rhs;
    aTrhs(1) = col2*rhs;

    sol(0) = ( a22 * aTrhs(0) - a12 * aTrhs(1)) / det;
    sol(1) = (-a12 * aTrhs(0) + a11 * aTrhs(1)) / det;
    return 0;
  }

  bool ValidBarCoord(double lami[3], double eps=1e-12)
  {
    return (lami[0]<=1.+eps && lami[0]>=0.-eps && lami[1]<=1.+eps && lami[1]>=0.-eps && lami[2]<=1.+eps && lami[2]>=0.-eps );
  }

  bool Mesh :: PointContainedIn2DElement(const netgen::Point<3> & p,
                                         double lami[3],
                                         SurfaceElementIndex ei,
                                         bool consider3D) const
  {
    Vec<3> col1, col2, col3;
    Vec<3> rhs, sol;
    const double eps = 1e-6;

    Array<Element2d> loctrigs;

    
    //SZ 
    if(surfelements[ei].GetType()==QUAD)
      {
        const Element2d & el = surfelements[ei];

        const auto & p1 = Point(el[0]); 
        const auto & p2 = Point(el[1]);
        const auto & p3 = Point(el[2]);
        const auto & p4 = Point(el[3]);

        if (el.GetOrder() > 1 || el.GetHpElnr() != -1) {
          netgen::Point<2> lam(0.5,0.5);
          Vec<3> rhs;
          Vec<2> deltalam;

          netgen::Point<3> x;
          Mat<3,2> Jac;
          double delta = 1.;
          const int maxits = 30;
          int i = 0;
          while(delta > 1e-16 && i < maxits)
            {
              curvedelems->CalcSurfaceTransformation(lam,ei,x,Jac);
              rhs = p - x;
              Jac.Solve(rhs,deltalam);
              lam += deltalam;
              delta = deltalam.Length2();
              i++;
            }
          if(i == maxits)
            return false;
          lami[0] = lam[0];
          lami[1] = lam[1];
          if(lami[0] < -eps || lami[0] > 1+eps || lami[1] < -eps || lami[1] > 1+eps)
            return false;
          return true;
        }

        // Coefficients of Bilinear Mapping from Ref-Elem to global Elem
        // X = a + b x + c y + d x y 
        Vec<3> a (p1);
        Vec<3> b = Vec<3>(p2) - a;
        Vec<3> c = Vec<3>(p4) - a;
        Vec<3> d = Vec<3>(p3) - a - b - c;

        /*cout << "p = " << p << endl;
        cout << "p1 = " << p1 << endl;
        cout << "p2 = " << p2 << endl;
        cout << "p3 = " << p3 << endl;
        cout << "p4 = " << p4 << endl;

        cout << "a = " << a << endl;
        cout << "b = " << b << endl;
        cout << "c = " << c << endl;
        cout << "d = " << d << endl;*/


        Vec<3> pa = Vec<3>(p) - a;
        double dxb = d(0)*b(1)-d(1)*b(0);
        double dxc = d(0)*c(1)-d(1)*c(0);
        double bxc = b(0)*c(1)-b(1)*c(0);
        double bxpa = b(0)*pa(1)-b(1)*pa(0);
        double cxpa = c(0)*pa(1)-c(1)*pa(0);
        double dxpa = d(0)*pa(1)-d(1)*pa(0);

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
        // double eps = 1.E-12;
        double c1,c2,r;

        //First check if point is "exactly" a vertex point
        Vec<3> d1 = p-p1;
        Vec<3> d2 = p-p2;
        Vec<3> d3 = p-p3;
        Vec<3> d4 = p-p4;

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
            Vec<2> sol;
            SolveLinearSystemLS (b, c, Vec<3>(p)-a, sol);
            lami[0] = sol(0);
            lami[1] = sol(1);
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
        

        Vec<3> dp13 = p3-p1;
        Vec<3> dp24 = p4-p2;
        double d1 = dp13.Length2();
        double d2 = dp24.Length2();

        // if(fabs(d.X()) <= eps && fabs(d.Y())<= eps)
        //if (d.Length2() < sqr(eps))
        if (d.Length2() < sqr(eps)*d1 && d.Length2() < sqr(eps)*d2)
          {
            //Solve Linear System
            Vec<2> sol;
            SolveLinearSystemLS (b, c, Vec<3>(p)-a, sol);
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
                Vec<3> n = Cross(b,c);
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
        //        SurfaceElement(element).GetTets (loctets);
        loctrigs.SetSize(1);
        loctrigs[0] = surfelements[ei];



        for (int j = 0; j < loctrigs.Size(); j++)
          {
            const Element2d & el = loctrigs[j];


            const auto & p1 = Point(el[0]);
            const auto & p2 = Point(el[1]);
            const auto & p3 = Point(el[2]);
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
            //col3 = Vec<3>(0, 0, 1);
            rhs = p - p1;

            // int retval = 
            SolveLinearSystem (col1, col2, col3, rhs, sol);

            //(*testout) << "retval " << retval << endl;

            //(*testout) << "col1 " << col1 << " col2 " << col2 << " col3 " << col3 << " rhs " << rhs << endl;
            //(*testout) << "sol " << sol << endl;

            if (surfelements[ei].GetType() ==TRIG6 || curvedelems->IsCurved(ei))
              {
                // netgen::Point<2> lam(1./3,1./3);
                netgen::Point<2> lam(sol(0), sol(1));
                if(surfelements[ei].GetType() != TRIG6)
                  {
                    lam[0] = 1-sol(0)-sol(1);
                    lam[1] = sol(0);
                  }
                Vec<3> rhs;
                Vec<2> deltalam;
                netgen::Point<3> x;
                Mat<3,2> Jac,Jact;
                
                double delta=1;
                
                // bool retval;
                
                int i = 0;
                
                const int maxits = 30;
                while(delta > 1e-16 && i<maxits)
                  {
                    curvedelems->CalcSurfaceTransformation(lam,ei,x,Jac);
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
                
                sol(0) = lam(0);
                sol(1) = lam(1);

                if (surfelements[ei].GetType() !=TRIG6 )
                  {
                    sol(2) = sol(0);
                    sol(0) = sol(1);
                    sol(1) = 1.0 - sol(2) - sol(0);
                  }

              }
            if (sol(0) >= -eps && sol(1) >= -eps && 
                sol(0) + sol(1) <= 1+eps)
              {
                if(!consider3D || (sol(2) >= -eps && sol(2) <= eps))
                  {
                    lami[0] = sol(0);
                    lami[1] = sol(1);
                    lami[2] = sol(2);

                    return true;
                  }
              }
          }
      }

    return false;

  }




  bool Mesh :: PointContainedIn3DElement(const netgen::Point<3> & p,
                                         double lami[3],
                                         ElementIndex ei,
                                         double eps) const
  {
    //bool oldresult = PointContainedIn3DElementOld(p,lami,element);
    //(*testout) << "old result: " << oldresult
    //       << " lam " << lami[0] << " " << lami[1] << " " << lami[2] << endl;

    //if(!curvedelems->IsCurved(ei))
    //  return PointContainedIn3DElementOld(p,lami,ei);
    auto el = volelements[ei];

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
        curvedelems->CalcElementTransformation(lam,ei,x,Jac);
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

    for (int i = 0; i < 3; i++)
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



  bool Mesh :: PointContainedIn3DElementOld(const netgen::Point<3> & p,
                                            double lami[3],
                                            ElementIndex element,
                                            double eps) const
  {
    Vec<3> col1, col2, col3;
    Vec<3> rhs, sol;

    Array<Element> loctets;

    VolumeElement(element).GetTets (loctets);

    for (int j = 1; j <= loctets.Size(); j++)
      {
        const Element & el = loctets[j-1];

        const auto & p1 = Point(el[0]);
        const auto & p2 = Point(el[1]);
        const auto & p3 = Point(el[2]);
        const auto & p4 = Point(el[3]);

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

        if (sol(0) >= -eps && sol(1) >= -eps && sol(2) >= -eps &&
            sol(0) + sol(1) + sol(2) <= 1+eps)
          {
            Array<ElementTet> loctetsloc;
            Array<netgen::Point<3> > pointsloc;

            VolumeElement(element).GetTetsLocal (loctetsloc);
            VolumeElement(element).GetNodesLocalNew (pointsloc);

            const ElementTet & le = loctetsloc[j-1];


            auto locp = [&](int j) -> const netgen::Point<3> &
              { return pointsloc[le.PNum(j).Nr1()-1]; };
            const auto & lp1 = locp(1);
            netgen::Point<3> pp =
              lp1
              + sol(0) * (locp(2) - lp1)
              + sol(1) * (locp(3) - lp1)
              + sol(2) * (locp(4) - lp1);

            lami[0] = pp(0);
            lami[1] = pp(1);
            lami[2] = pp(2);
            return true;
          }
      }
    return false;
  }


  ElementIndex Mesh :: GetElementOfPoint (const netgen::Point<3> & p,
                                          double* lami,
                                          bool build_searchtree,
                                          int index,
                                          bool allowindex,
                                          double tol) const
  {
    if(index != -1) 
      {
        Array<int> dummy(1);
        dummy[0] = index;
        return GetElementOfPoint(p,lami,dummy,build_searchtree,allowindex, tol);
      }
    else
      return GetElementOfPoint(p,lami,nullopt,build_searchtree,allowindex, tol);
  }




  ElementIndex Mesh :: GetElementOfPoint (const netgen::Point<3> & p,
                                          double* lami,
                                          std::optional<FlatArray<int>> indices,
                                          bool build_searchtree,
                                          bool allowindex,
                                          double tol) const
  {
    if (build_searchtree)
      const_cast<Mesh&>(*this).BuildElementSearchTree (3);
    return Find3dElement(*this, p, lami, indices, elementsearchtree_vol.get(), allowindex, tol);
  }



  SurfaceElementIndex Mesh ::
  GetSurfaceElementOfPoint (const netgen::Point<3> & p,
                            double* lami,
                            bool build_searchtree,
                            int index,
                            bool allowindex) const
  {
    if(index != -1)
      {
        Array<int> dummy(1);
        dummy[0] = index;
        return GetSurfaceElementOfPoint(p,lami,dummy,build_searchtree,allowindex);
      }
    else
      return GetSurfaceElementOfPoint(p,lami,nullopt,build_searchtree,allowindex);
  }

  SurfaceElementIndex Mesh ::
  GetSurfaceElementOfPoint (const netgen::Point<3> & p,
                            double* lami,
                            std::optional<FlatArray<int>> indices,
                            bool build_searchtree,
                            bool allowindex) const
  {
    if (build_searchtree)
      const_cast<Mesh&>(*this).BuildElementSearchTree(2);
    return Find2dElement(*this, p, lami, indices, elementsearchtree_surf.get(), allowindex);
  }


  void Mesh::GetIntersectingVolEls(const netgen::Point<3>& p1, const netgen::Point<3>& p2, 
                                   Array<ElementIndex> & locels) const
  {
    elementsearchtree_vol->GetIntersecting (p1, p2, locels);
  }

  void Mesh :: SplitIntoParts()
  {
    // int i, j, dom;
    int ne = GetNE();
    int np = GetNP();
    int nse = GetNSE();

    BitArray surfused(nse+1);
    TBitArray<PointIndex> pused (np);

    surfused.Clear();

    int dom = 0;

    while (1)
      {
        int cntd = 1;

        dom++;

        pused.Clear();

        int found = 0;
        for (SurfaceElementIndex i : T_Range<SurfaceElementIndex>(nse))
          if (!surfused.Test(i.Nr1()))
            {
              (*this)[i].SetIndex (FaceRegionIndex::FromNr1(dom));
              for (int j = 0; j < 3; j++)
                pused.SetBit ((*this)[i][j]);
              found = 1;
              cntd = 1;
              surfused.SetBit(i.Nr1());
              break;
            }

        if (!found)
          break;

        int change;
        do
          {
            change = 0;
            for (SurfaceElementIndex i : T_Range<SurfaceElementIndex>(nse))
              {
                int is = 0, isnot = 0;
                for (int j = 0; j < 3; j++)
                  if (pused.Test((*this)[i][j]))
                    is = 1;
                  else
                    isnot = 1;

                if (is && isnot)
                  {
                    change = 1;
                    for (int j = 0; j < 3; j++)
                      pused.SetBit ((*this)[i][j]);
                  }

                if (is) 
                  {
                    if (!surfused.Test(i.Nr1()))
                      {
                        surfused.SetBit(i.Nr1());
                        (*this)[i].SetIndex (FaceRegionIndex::FromNr1(dom));
                        cntd++;
                      }
                  }
              }


            for (ElementIndex i : T_Range<ElementIndex>(ne))
              {
                int is = 0, isnot = 0;
                for (int j = 0; j < 4; j++)
                  if (pused.Test((*this)[i][j]))
                    is = 1;
                  else
                    isnot = 1;

                if (is && isnot)
                  {
                    change = 1;
                    for (int j = 0; j < 4; j++)
                      pused.SetBit ((*this)[i][j]);
                  }

                if (is)
                  {
                    (*this)[i].SetIndex (VolumeRegionIndex::FromNr1(dom));
                  }
              }
          }
        while (change);

        PrintMessage (3, "domain ", dom, " has ", cntd, " surfaceelements");
      }

    /*
      Regions<2>().SetSize (dom);
      for (i = 1; i <= dom; i++)
      {
      Regions<2>().Elem(i).surfnr = 0;
      Regions<2>().Elem(i).domin = i;
      Regions<2>().Elem(i).domout = 0;
      }
    */
    ClearFaceDescriptors();
    for (int i = 1; i <= dom; i++)
      AddFaceDescriptor (FaceRegion (0, i, 0, 0));
    CalcSurfacesOfNode();
    timestamp = NextTimeStamp();
  }

  void Mesh :: SplitSeparatedFaces ()
  {
    auto seg_fdi = [this](const Segment& s) -> int {
      const Mesh & self = *this;
      if (self.HasEdgeDescriptor(s))
        { auto fdi = self.Regions<1>()[s.GetIndex()].GetIndex(); if (fdi.IsValid()) return fdi.Nr1(); }
      return -1;
    };
    PrintMessage (3, "SplitSeparateFaces");
    int fdi;
    int np = GetNP();

    TBitArray<PointIndex> usedp(np);
    Array<SurfaceElementIndex> els_of_face;

    fdi = 1;
    while (fdi <= GetNFD())
      {
        GetSurfaceElementsOfFace (fdi, els_of_face);

        if (els_of_face.Size() == 0)
        {
            fdi++;
            continue;
        }

        SurfaceElementIndex firstel = els_of_face[0];

        usedp.Clear();
        for (int j = 0; j < SurfaceElement(firstel).GetNP(); j++)
          usedp.SetBit (SurfaceElement(firstel)[j]);

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
                    usedp.SetBit (el[j]);
              }
          }
        while (changed);

        int nface = 0;
        for (int i = 0; i < els_of_face.Size(); i++)
          {
            Element2d & el = SurfaceElement(els_of_face[i]);

            int hasno = 0;
            for (int j = 0; j < el.GetNP(); j++)
              if (!usedp.Test(el[j]))
                hasno = 1;

            if (hasno)
              {
                if (!nface)
                  {
                    FaceRegion nfd = GetFaceDescriptor(FaceRegionIndex::FromNr1(fdi));
                    nface = AddFaceDescriptor (nfd).Nr1();
                  }

                el.SetIndex (FaceRegionIndex::FromNr1(nface));
              }
          }

        // reconnect list
        if (nface)
          {
            Regions<2>()[FaceRegionIndex::FromNr1(nface)].firstelement = SurfaceElementIndex::INVALID;
            Regions<2>()[FaceRegionIndex::FromNr1(fdi)].firstelement = SurfaceElementIndex::INVALID;

            for (int i = 0; i < els_of_face.Size(); i++)
              {
                auto ind = SurfaceElement(els_of_face[i]).GetIndex();
                SurfaceElement(els_of_face[i]).next = Regions<2>()[ind].firstelement;
                Regions<2>()[ind].firstelement = els_of_face[i];
              }

            // map the segments - also create per-face EDs so edsi stays in sync
            map<pair<int,int>, int> split_ed_cache;
            for(auto& seg : segments)
              if(!usedp.Test(seg[0]) || !usedp.Test(seg[1]))
                {
                  if(seg_fdi(seg) == fdi)
                    {
                      if (HasEdgeDescriptor(seg))
                        {
                          auto key = make_pair(seg.GetIndex().Nr1(), nface);
                          auto it = split_ed_cache.find(key);
                          if (it != split_ed_cache.end())
                            {
                              seg.SetIndex(EdgeRegionIndex::FromNr1(it->second));
                            }
                          else
                            {
                              EdgeRegion new_ed = Regions<1>()[seg.GetIndex()];
                              new_ed.SetIndex(FaceRegionIndex::FromNr1(nface));
                              auto new_edsi = AddEdgeDescriptor(new_ed);
                              split_ed_cache[key] = new_edsi.Nr1();
                              seg.SetIndex(new_edsi);
                            }
                        }
                    }
                }
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
      FaceRegion nfd = GetFaceDescriptor(fdi);
      nface = AddFaceDescriptor (nfd).Nr1();
      }

      el.SetIndex (nface);
      }
      }
      fdi++;
      }
    */
  }

  void Mesh :: ZRefine(const string& name, const Array<double>& slices)
  {
    auto nr = GetIdentifications().GetNr(name);
    auto& identpts = GetIdentifications().GetIdentifiedPoints();

    UpdateTopology();

    std::map<std::pair<PointIndex, PointIndex>,
             Array<PointIndex>> inserted_points;
    TBitArray<PointIndex> mapped_points(GetNV());
    mapped_points = false;

    // Add new points
    for(auto [hash, dummy] : identpts)
      {
        auto [hash_pts, hash_nr] = hash;
        if(hash_nr != nr)
          continue;
        // auto& ipts = inserted_points[{p1p2[0], p1p2[1]}];
        auto& ipts = inserted_points[ { hash_pts[0], hash_pts[1] }];
        auto p1 = Point(hash_pts[0]);
        auto p2 = Point(hash_pts[1]);
        ipts.Append(hash_pts[0]);
        mapped_points.SetBit(hash_pts[0]);
        for(auto slice : slices)
          {
            auto np = p1 + slice * (p2-p1);
            auto npi = AddPoint(np);
            ipts.Append(npi);
          }
        ipts.Append(hash_pts[1]);
      }

    // Store offset-point identifications for curving
    {
      auto & ident = GetIdentifications();
      int offset_nr = ident.GetNr("offset_points");
      ident.SetType(offset_nr, Identifications::OFFSET_POINT);
      for (const auto& [pair, chain] : inserted_points)
        {
          PointIndex base_pi = pair.first;
          for (auto i : Range(size_t(1), chain.Size()-1))  // skip endpoints
            ident.Add(chain[i], base_pi, offset_nr);  // inverse: offset -> base
        }
    }

    // Split segments
    for(auto si : Range(segments))
      {
        auto& seg = segments[si];
        // Copy segment, as reference above might get invalidated in AddSegment()
        auto reference_seg = seg;
        auto p1 = seg[0];
        auto p2 = seg[1];

        auto c1 = inserted_points.count({p1, p2});
        auto c2 = inserted_points.count({p2, p1});

        if(c1 == 0 && c2 == 0)
          continue;

        if(c2)
          Swap(p1,p2);

        const auto& ipts = inserted_points[{p1,p2}];
        if(c2)
          seg[1] = ipts[ipts.Size()-2];
        else
          seg[1] = ipts[1];
        for(auto i : Range(size_t(1), ipts.Size()-1))
          {
            Segment snew = reference_seg;
            if(c2)
              {
                snew[0] = ipts[ipts.Size()-1-i];
                snew[1] = ipts[ipts.Size()-2-i];
              }
            else
              {
                snew[0] = ipts[i];
                snew[1] = ipts[i+1];
              }
            AddSegment(snew);
          }
      }

    TBitArray<SurfaceElementIndex> sel_done(surfelements.Size());
    sel_done = false;

    // Split surface elements
    auto p2sel = CreatePoint2SurfaceElementTable();
    for(const auto& [pair, inserted] : inserted_points)
      {
        for(auto si : p2sel[pair.first])
          {
            if(sel_done[si])
              continue;
            sel_done.SetBit(si);
            auto sel = surfelements[si];
            map<PointIndex, Array<PointIndex>> mapped_points;
            int nmapped = 0;
            for(auto i : Range(sel.GetNP()))
              {
                auto p1 = sel[i];
                auto p2 = sel[(i+1)%sel.GetNP()];
                auto c1 = inserted_points.count({p1, p2});
                auto c2 = inserted_points.count({p2, p1});
                if(c1 == 0 && c2 == 0)
                  continue;
                if(c2)
                  Swap(p1, p2);
                auto& ipts = inserted_points[{p1, p2}];
                auto& a1 = mapped_points[p1];
                auto& a2 = mapped_points[p2];
                a1 = ipts.Range(0, ipts.Size()-1);
                a2 = ipts.Range(1, ipts.Size());
                nmapped = ipts.Size()-1;
              }
            for(auto i : Range(nmapped))
              {
                Element2d nsel = sel;
                for(auto& pi : nsel.PNums())
                  if(mapped_points.count(pi))
                    pi = mapped_points[pi][i];
                AddSurfaceElement(nsel);
              }
            if(nmapped)
              surfelements[si].Delete();
          }
      }

    // Split volume elements
    TBitArray<ElementIndex> vol_done(volelements.Size());
    vol_done = false;
    auto p2el = CreatePoint2ElementTable(); // mapped_points);
    for(const auto& [pair, inserted] : inserted_points)
      {
        for(auto ei : p2el[pair.first])
          {
            if(vol_done[ei])
              continue;
            vol_done.SetBit(ei);
            auto el = volelements[ei];
            map<PointIndex, Array<PointIndex>> mapped_points;
            int nmapped = 0;
            // Array<int> eledges;
            // topology.GetElementEdges(ei+1, eledges);
            // for(auto edgei : eledges)
            for(auto edgei : topology.GetEdges(ElementIndex(ei)))
              {
                // int p1, p2;
                // topology.GetEdgeVertices(edgei+1, p1, p2);
                auto [p1, p2] = topology.GetEdgeVertices(edgei);
                auto c1 = inserted_points.count({p1, p2});
                auto c2 = inserted_points.count({p2, p1});
                if(c1 == 0 && c2 == 0)
                  continue;
                if(c2)
                  Swap(p1, p2);
                auto& ipts = inserted_points[{p1, p2}];
                auto& a1 = mapped_points[p1];
                auto& a2 = mapped_points[p2];
                a1 = ipts.Range(0, ipts.Size()-1);
                a2 = ipts.Range(1, ipts.Size());
                nmapped = ipts.Size()-1;
              }

            for(auto i : Range(nmapped))
              {
                Element nel (el);
                for(auto& pi : nel.PNums())
                  if(mapped_points.count(pi))
                    pi = mapped_points[pi][i];
                AddVolumeElement(nel);
              }
            if(nmapped)
              volelements[ei].Delete();
          }
      }

    Compress();
    SetNextMajorTimeStamp();
  }

  void Mesh :: RebuildSurfaceElementLists ()
  {
    static Timer t("Mesh::LinkSurfaceElements"); RegionTimer reg (t);    
    
    for (int i = 0; i < Regions<2>().Size(); i++)
      Regions<2>()[FaceRegionIndex::FromNr0(i)].firstelement = SurfaceElementIndex::INVALID;
    for (int i = surfelements.Size()-1; i >= 0; i--)
      {
        SurfaceElementIndex sei = SurfaceElementIndex::FromNr0(i);
        auto ind = surfelements[sei].GetIndex();
        surfelements[sei].next = Regions<2>()[ind].firstelement;
        Regions<2>()[ind].firstelement = sei;
      }
  }

  void Mesh :: GetSurfaceElementsOfFace (FaceRegionIndex fi, Array<SurfaceElementIndex> & sei) const
  {
    static Timer timer("GetSurfaceElementsOfFace");
    RegionTimer reg (timer);

    if(!fi.IsValid())
    {
        sei.SetSize(GetNSE());
        ParallelForRange( IntRange(GetNSE()), [&sei] (auto myrange)
            {
                for(auto i : myrange)
                    sei[i] = SurfaceElementIndex::FromNr0(i);
            });
        return;
    }

     sei.SetSize(0);

     SurfaceElementIndex si = Regions<2>()[fi].firstelement;
     while (si.IsValid())
     {
       if ( (*this)[si].GetIndex() == fi && (*this)[si][0].IsValid() &&
            !(*this)[si].IsDeleted() )
        {
           sei.Append (si);
        }

        si = (*this)[si].next;
     }
  }




  void Mesh :: CalcMinMaxAngle (double badellimit, double * retvalues) 
  {
    double phimax = 0, phimin = 10;
    double facephimax = 0, facephimin = 10;
    int illegaltets = 0, negativetets = 0, badtets = 0;

    // for (int i = 1; i <= GetNE(); i++)
    for (ElementIndex ei : Range(VolumeElements()))
      {
        int badel = 0;

        auto el = VolumeElement(ei);

        if (el.GetType() != TET)
          {
            VolumeElement(ei).Flags().badel = 0;
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
            (*testout) << "illegal tet: " << ei << " ";
            for (int j = 0; j < el.GetNP(); j++)
              (*testout) << el[j] << " ";
            (*testout) << endl;
          }


        // angles between faces
        for (int lpi1 = 1; lpi1 <= 3; lpi1++)
          for (int lpi2 = lpi1+1; lpi2 <= 4; lpi2++)
            {
              int lpi3 = 1;
              while (lpi3 == lpi1 || lpi3 == lpi2)
                lpi3++;
              int lpi4 = 10 - lpi1 - lpi2 - lpi3;

              const auto & p1 = Point (el.PNum(lpi1));
              const auto & p2 = Point (el.PNum(lpi2));
              const auto & p3 = Point (el.PNum(lpi3));
              const auto & p4 = Point (el.PNum(lpi4));

              Vec<3> n(p1, p2);
              n /= n.Length();
              Vec<3> v1(p1, p3);
              Vec<3> v2(p1, p4);

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
        for (int j = 1; j <= 4; j++)
          {
            Element2d face(TRIG);
            el.GetFace (j, face);
            for (int lpi1 = 1; lpi1 <= 3; lpi1++)
              {
                int lpi2 = lpi1 % 3 + 1;
                int lpi3 = lpi2 % 3 + 1;

                const auto & p1 = Point (el.PNum(lpi1));
                const auto & p2 = Point (el.PNum(lpi2));
                const auto & p3 = Point (el.PNum(lpi3));

                Vec<3> v1(p1, p2);
                Vec<3> v2(p1, p3);
                double cosphi = (v1 * v2) / (v1.Length() * v2.Length());
                double phi = acos (cosphi);
                if (phi > facephimax) facephimax = phi;
                if (phi < facephimin) facephimin = phi;

                if ((180/M_PI) * phi > badellimit)
                  badel = 1;

              }
          }


        VolumeElement(ei).Flags().badel = badel;
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


  int Mesh :: MarkIllegalElements (int domain)
  {
    if(!boundaryedges)
      BuildBoundaryEdges();

    atomic<int> cnt = 0;
    ParallelForRange( Range(volelements), [&] (auto myrange)
    {
      int cnt_local = 0;
      for (auto el : volelements.Range(myrange))
        if ((domain==0 || el.GetIndex().Nr1() == domain) && !LegalTet (el))
          cnt_local++;
      cnt += cnt_local;
    });
    return cnt;
  }

  // #ifdef NONE
  //   void Mesh :: AddIdentification (int pi1, int pi2, int identnr)
  //   {
  //     IVec<2> pair(pi1, pi2);
  //     //  pair.Sort();
  //     identifiedpoints->Set (pair, identnr);
  //     if (identnr > maxidentnr)
  //       maxidentnr = identnr;
  //     timestamp = NextTimeStamp();
  //   }

  //   int Mesh :: GetIdentification (int pi1, int pi2) const
  //   {
  //     IVec<2> pair(pi1, pi2);
  //     if (identifiedpoints->Used (pair))
  //       return identifiedpoints->Get(pair);
  //     else
  //       return 0;
  //   }

  //   int Mesh :: GetIdentificationSym (int pi1, int pi2) const
  //   {
  //     IVec<2> pair(pi1, pi2);
  //     if (identifiedpoints->Used (pair))
  //       return identifiedpoints->Get(pair);

  //     pair = IVec<2> (pi2, pi1);
  //     if (identifiedpoints->Used (pair))
  //       return identifiedpoints->Get(pair);

  //     return 0;
  //   }


  //   void Mesh :: GetIdentificationMap (int identnr, Array<int> & identmap) const
  //   {
  //     int i, j;

  //     identmap.SetSize (GetNP());
  //     for (i = 1; i <= identmap.Size(); i++)
  //       identmap.Elem(i) = 0;

  //     for (i = 1; i <= identifiedpoints->GetNBags(); i++)
  //       for (j = 1; j <= identifiedpoints->GetBagSize(i); j++)
  //    {
  //      IVec<2> i2;
  //      int nr;
  //      identifiedpoints->GetData (i, j, i2, nr);

  //      if (nr == identnr)
  //        {
  //          identmap.Elem(i2[0]) = i2[1];
  //        }
  //    }
  //   }


  //   void Mesh :: GetIdentificationPairs (int identnr, Array<IVec<2>> & identpairs) const
  //   {
  //     int i, j;

  //     identpairs.SetSize(0);

  //     for (i = 1; i <= identifiedpoints->GetNBags(); i++)
  //       for (j = 1; j <= identifiedpoints->GetBagSize(i); j++)
  //    {
  //      IVec<2> i2;
  //      int nr;
  //      identifiedpoints->GetData (i, j, i2, nr);

  //      if (identnr == 0 || nr == identnr)
  //        identpairs.Append (i2);
  //    }
  //   }
  // #endif

  int Mesh::IdentifyPeriodicBoundaries(const string& id_name,
                                       const string &s1,
                                       const Transformation<3> &mapping,
                                       double pointTolerance)
  {
    auto nr = ident->GetNr(id_name);
    ident->SetType(nr, Identifications::PERIODIC);
    // double lami[4];
    set<PointIndex> identified_points;
    if(pointTolerance < 0.)
      {
        netgen::Point<3> pmin, pmax;
        GetBox(pmin, pmax);
        pointTolerance = 1e-8 * (pmax-pmin).Length();
      }
    size_t nse = GetDimension() == 3 ? surfelements.Size() : segments.Size();
    for(auto nr : Range(nse))
      {
        // in 3d these are surface elements, in 2d segments
        SurfaceElementIndex sei = SurfaceElementIndex::FromNr0(nr);
        SegmentIndex segi = SegmentIndex::FromNr0(nr);
        string_view name = GetDimension() == 3 ? GetRegionName(surfelements[sei]) : GetRegionName(segments[segi]);
        if(name != s1)
          continue;

        const auto& pnums = GetDimension() == 3 ? surfelements[sei].PNums() :
          segments[segi].PNums();
        for(const auto& pi : pnums)
          {
            if(identified_points.find(pi) != identified_points.end())
              continue;
            auto pt = (*this)[pi];
            auto mapped_pt = mapping(pt);
            bool found = false;
            for(auto other_pi : Range(points))
              {
                if((mapped_pt - (*this)[other_pi]).Length() < pointTolerance)
                  {
                    identified_points.insert(pi);
                    ident->Add(pi, other_pi, nr);
                    found = true;
                    break;
                  }
              }
            if(!found)
              {
                cout << "point coordinates = " << pt << endl;
                cout << "mapped coordinates = " << mapped_pt << endl;
                throw Exception("Did not find mapped point with nr " + ToString(pi) + ", are you sure your mesh is periodic?");
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
  void Mesh :: AddPointCurvePoint(const netgen::Point<3> & pt) const
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

  netgen::Point<3> & Mesh :: GetPointCurvePoint(int curve, int n) const
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
    /*
    for (const Element & el : VolumeElements())
      for (PointIndex v : el.Vertices())
        if (v > numvertices) numvertices = v;
        
    for (const Element2d & el : SurfaceElements())
      for (PointIndex v : el.Vertices())
        if (v > numvertices) numvertices = v;

    numvertices += 1-PointIndex::BASE;
    */
    numvertices = -1;
    numvertices =
      ParallelReduce (VolumeElements().Size(),
                      [&](size_t nr)
                      {
                        return Max((*this)[ElementIndex::FromNr0(nr)].Vertices()) - IndexBASE<PointIndex>();
                      },
                      [](auto a, auto b) { return a > b ?  a : b; },
                      numvertices);
    numvertices =
      ParallelReduce (SurfaceElements().Size(),
                      [&](size_t nr)
                      {
                        return Max((*this)[SurfaceElementIndex::FromNr0(nr)].Vertices()) - IndexBASE<PointIndex>();
                      },
                      [](auto a, auto b) { return a > b ?  a : b; },
                      numvertices);
    numvertices =
      ParallelReduce (LineSegments().Size(),
                      [&](size_t nr)
                      {
                        return Max((*this)[SegmentIndex::FromNr0(nr)].Vertices()) - IndexBASE<PointIndex>();
                      },
                      [](auto a, auto b) { return a > b ?  a : b; },
                      numvertices);
    numvertices += 1;
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
      for (PointIndex i = mlold+IndexBASE<PointIndex>(); 
           i < np+IndexBASE<PointIndex>(); i++)
        {
          mlbetweennodes[i][0].Invalidate();
          mlbetweennodes[i][1].Invalidate();
        }

    GetIdentifications().SetMaxPointNr (np + PointIndex::BASE-1);
  }


  Table<ElementIndex, PointIndex> Mesh :: CreatePoint2ElementTable(std::optional<TBitArray<PointIndex>> points, int domain) const
  {
    static Timer timer("Mesh::CreatePoint2VolumeElementTable"); RegionTimer rt(timer);
    
    if(points)
      {
        const auto & free_points = *points;
        return ngcore::CreateSortedTable<ElementIndex, PointIndex>( volelements.Range(),
               [&](auto & table, ElementIndex ei)
               {
                 const auto & el = (*this)[ei];
                 if(el.IsDeleted())
                     return;

                 if(domain && el.GetIndex().Nr1() != domain)
                     return;

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
                 if(el.IsDeleted())
                     return;

                 if(domain && el.GetIndex().Nr1() != domain)
                     return;

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


  CompressedTable<SurfaceElementIndex, PointIndex> Mesh :: CreateCompressedPoint2SurfaceElementTable( FaceRegionIndex fi ) const
  {
    static Timer timer("Mesh::CreatePoint2SurfaceElementTable"); RegionTimer rt(timer);

    CompressedTableCreator<SurfaceElementIndex, PointIndex> creator;
    
    if(!fi.IsValid())
      {
        for ( ; !creator.Done(); creator++)
          for (auto sei : SurfaceElements().Range())
            for (auto pi : (*this)[sei].PNums())
              creator.Add(pi, sei);
      }
    else
      {
        Array<SurfaceElementIndex> face_els;
        GetSurfaceElementsOfFace(fi, face_els);

        for ( ; !creator.Done(); creator++)
          for (auto sei : face_els)
            for (auto pi : (*this)[sei].PNums())
              creator.Add(pi, sei);
      }


    auto compressed_table = creator.MoveTable();
    
    for (auto row : compressed_table.GetTable())
      QuickSort (row);
    
    return compressed_table;
  }



  

  bool Mesh :: PureTrigMesh (int faceindex) const
  {
    // if (!faceindex) return !mparam.quad;
    
    if (!faceindex)
      {
        for (SurfaceElementIndex i : SurfaceElements().Range())
          if ((*this)[i].GetNP() != 3)
            return false;
        return true;
      }

    for (SurfaceElementIndex i : SurfaceElements().Range())
      if ((*this)[i].GetIndex().Nr1() == faceindex &&
          (*this)[i].GetNP() != 3)
        return false;
    return true;
  }

  bool Mesh :: PureTetMesh () const
  {
    for (ElementIndex ei : VolumeElements().Range())
      if (VolumeElement(ei).GetNP() != 4)
        return 0;
    return 1;
  }

  void Mesh :: UpdateTopology ()
  {
    static Timer t("Update Topology"); RegionTimer reg(t);
    ComputeNVertices();
    topology.Update();
    static Timer t_call_update_clusters("call update clusters"); t_call_update_clusters.Start();
    clusters->Update();
    t_call_update_clusters.Stop();
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


    for (SurfaceElementIndex sei : SurfaceElements().Range())
      (*this)[sei].SetCurved (GetCurvedElements().IsCurved (sei));
    for (ElementIndex ei : VolumeElements().Range())
      (*this)[ei].SetCurved (GetCurvedElements().IsCurved (ei));
    
    SetNextMajorTimeStamp();
  }

  void Mesh :: BuildCurvedElements (int aorder)
  {
    if (!GetGeometry())
      throw NgException ("don't have a geometry for mesh curving");
    
    GetCurvedElements().BuildCurvedElements (&GetGeometry()->GetRefinement(), aorder, false);


    for (SurfaceElementIndex sei : SurfaceElements().Range())
      (*this)[sei].SetCurved (GetCurvedElements().IsCurved (sei));
    for (ElementIndex ei : VolumeElements().Range())
      (*this)[ei].SetCurved (GetCurvedElements().IsCurved (ei));
    
    SetNextMajorTimeStamp();
  }

  void Mesh :: SplitFacesByAdjacentDomains ()
  {
    UpdateTopology();
    std::map<std::tuple<int, int, int>, int> face_doms_2_new_face;
    int nfaces = FaceDescriptors().Size();
    Array<bool> first_visit(nfaces);
    first_visit = true;

    for (auto sei : Range(SurfaceElements()))
      {
        ElementIndex eli0, eli1;
        GetTopology().GetSurface2VolumeElement(sei, eli0, eli1);
        if(!eli0.IsValid())
          continue;
        auto & sel = (*this)[sei];
        int face = sel.GetIndex().Nr1();
        int domin = (*this)[eli0].GetIndex().Nr1();
        int domout = eli1.IsValid() ? (*this)[eli1].GetIndex().Nr1() : 0;
        if(domin < domout)
          swap(domin, domout);

        auto key = std::make_tuple(face, domin, domout);
        if(face_doms_2_new_face.find(key) == face_doms_2_new_face.end())
          {
            {
              auto & fd = FaceDescriptors()[FaceRegionIndex::FromNr1(face)];
              if(domout == 0 && min(fd.DomainIn(), fd.DomainOut()) > 0)
                continue;
            }
            if(!first_visit[face-1]) {
              nfaces++;
              FaceRegion new_fd = FaceDescriptors()[FaceRegionIndex::FromNr1(face)];
              new_fd.bcprop = nfaces;
              new_fd.domin = domin;
              new_fd.domout = domout;
              AddFaceDescriptor(new_fd);
              SetBCName(nfaces-1, new_fd.GetBCName());
              face_doms_2_new_face[key] = nfaces;
            }
            else {
              face_doms_2_new_face[key] = face;
              auto & fd = FaceDescriptors()[FaceRegionIndex::FromNr1(face)];
              fd.domin = domin;
              fd.domout = domout;
            }
            first_visit[face-1] = false;
          }
          sel.SetIndex(FaceRegionIndex::FromNr1(face_doms_2_new_face[key]));
      }
    SetNextMajorTimeStamp();
    RebuildSurfaceElementLists ();
    CalcSurfacesOfNode();
    UpdateTopology();
  }

  shared_ptr<Mesh> Mesh :: GetSubMesh(string domains, string faces) const
  {
    // Copy the mesh into a new one, then delete unwanted elements
    // Unused points are deleted by the Compress() function at the end
    auto mesh_ptr = make_unique<Mesh>();
    auto & mesh = *mesh_ptr;
    mesh = (*this);

    auto ndomains = GetNDomains();
    auto nfaces = GetNFD();

    TBitArray<PointIndex> keep_point(GetNP());
    BitArray keep_face(nfaces+1);
    BitArray keep_domain(ndomains+1);
    keep_point.Clear();
    keep_face.Clear();
    keep_domain.Clear();

    regex regex_faces(faces);
    regex regex_domains(domains);

    if(dimension == 3) {
      for(auto dom : Range(ndomains))
        if(regex_match(mesh.GetMaterial(dom+1), regex_domains))
          keep_domain.SetBit(dom+1);

      for(auto fi : Range(nfaces))
      {
        auto & fd = mesh.FaceDescriptors()[FaceRegionIndex::FromNr0(fi)];
        if (regex_match(fd.GetBCName(), regex_faces) 
          || keep_domain[fd.DomainIn()] || keep_domain[fd.DomainOut()])
            keep_face.SetBit(fd.BCProperty());
      }
    }
    else {
      for(auto fi : Range(nfaces))
      {
        auto & fd = mesh.FaceDescriptors()[FaceRegionIndex::FromNr0(fi)];
        auto mat = GetMaterial(fd.BCProperty());
        if (regex_match(mat, regex_faces))
            keep_face.SetBit(fd.BCProperty());
      }
    }

    auto filter_elements = [&keep_point](auto & elements, auto & keep_region, auto region_of)
    {
      for (auto && el : elements)
      {
        if(keep_region[region_of(el)])
          for (auto pi : el.PNums())
            keep_point.SetBit(pi);
        else
          el.Delete();
      }
    };

    filter_elements(mesh.VolumeElements(), keep_domain, [](const ElementRef & el) { return el.GetIndex().Nr1(); });
    // keep_face is filled by BCProperty, tested here by descriptor number (they coincide for generated meshes)
    filter_elements(mesh.SurfaceElements(), keep_face, [](const Element2d & el) { return el.GetIndex().Nr1(); });

    // Keep line segments only if all points are kept
    // Check them in reverse order because they are deleted from the end
    auto nsegments = mesh.LineSegments().Size();
    for(auto i : Range(nsegments))
    {
      SegmentIndex segi = SegmentIndex::FromNr0(nsegments-i-1);
      auto seg = mesh[segi];
      bool keep = true;
      for(auto pi : seg.PNums())
        keep &= keep_point[pi];

      if(!keep)
        mesh.LineSegments().DeleteElement(segi);
    }

    // Check in reverse order because they are deleted from the end
    auto npointelements = mesh.pointelements.Size();
    for(auto i : Range(npointelements))
    {
      auto pel = mesh.pointelements[npointelements-i-1];
      if(!keep_point[pel.pnum])
        mesh.pointelements.DeleteElement(npointelements-i-1);
    }

    mesh.Compress();
    return mesh_ptr;
  }

  void Mesh :: SetMaterial (int domnr, const string & mat)
  {
    if (domnr < 1) throw RangeException("Illegal domain number ", domnr, 1, domnr);
    if (dimension == 2)
      {
        while (Regions<2>().Size() < domnr)
          {
            FaceRegion fd(0, 0, 0, 0);
            fd.SetBCProperty(Regions<2>().Size()+1);
            Regions<2>().Append(fd);
          }
        Regions<2>()[FaceRegionIndex::FromNr1(domnr)].SetBCName(mat);
      }
    else if (dimension == 1)
      {
        while (Regions<1>().Size() < domnr)
          Regions<1>().Append(EdgeRegion());
        Regions<1>()[EdgeRegionIndex::FromNr1(domnr)].SetName(mat);
      }
    else
      {
        auto & vols = Regions<3>();
        while (vols.Size() < domnr)
          vols.Append(VolumeRegion(defaultmat));   // set, like the old code
        vols[VolumeRegionIndex::FromNr1(domnr)].SetName(mat);
      }
  }

  string Mesh :: defaultmat = "default";
  string_view Mesh :: defaultmat_sv = "default";  
  const string & Mesh :: GetMaterial (int domnr) const
  {
    if (dimension == 2)
      return (domnr >= 1 && domnr <= Regions<2>().Size()) ? Regions<2>()[FaceRegionIndex::FromNr1(domnr)].GetBCName() : defaultmat;
    if (dimension == 1)
      return (domnr >= 1 && domnr <= Regions<1>().Size()) ? Regions<1>()[EdgeRegionIndex::FromNr1(domnr)].GetName() : defaultmat;
    return *GetMaterialPtr(domnr);
  }

  Array<optional<string>> Mesh :: DomainNames () const
  {
    Array<optional<string>> names;
    if (dimension == 3) return RegionNames<3>();
    if (dimension == 2)
      for (const auto & fd : Regions<2>()) names.Append(fd.GetBCName());
    else
      for (const auto & ed : Regions<1>()) names.Append(ed.GetName());
    return names;
  }

  void Mesh :: SetDomainNames (Array<optional<string>> names)
  {
    if (dimension == 3) { SetRegionNames<3>(names); return; }
    for (int i = 0; i < names.Size(); i++)
      if (names[i]) SetMaterial(i+1, *names[i]);
  }

  void Mesh ::SetNBCNames ( int nbcn )
  {
    if (dimension >= 2) return;   // boundary names live on the descriptors
    Regions<0>() = RegionArray<0>(nbcn);
  }

  void Mesh ::SetBCName ( int bcnr, const string & abcname )
  {
    if (bcnr < 0) throw RangeException("Illegal bc number ", bcnr, 0, bcnr);
    if (dimension == 3)
      {
        while (Regions<2>().Size() <= bcnr)
          {
            FaceRegion fd(0, 0, 0, 0);
            fd.SetBCProperty(Regions<2>().Size()+1);
            Regions<2>().Append(fd);
          }
        Regions<2>()[FaceRegionIndex::FromNr0(bcnr)].SetBCName(abcname);
      }
    else if (dimension == 2)
      {
        while (Regions<1>().Size() <= bcnr)
          Regions<1>().Append(EdgeRegion());
        Regions<1>()[EdgeRegionIndex::FromNr0(bcnr)].SetName(abcname);
      }
    else
      {
        auto & verts = Regions<0>();
        while (verts.Size() <= bcnr)
          verts.Append(VertexRegion(default_bc));
        verts[VertexRegionIndex::FromNr0(bcnr)].SetName(abcname);
      }
  }

  const string & Mesh ::GetBCName ( int bcnr ) const
  {
    return *GetBCNamePtr(bcnr);
  }

  // boundary names keyed by bc number (BCProperty), the format of files and archives:
  // the name of the first face descriptor with that bc number, empty if nothing is named
  Array<string> Mesh :: BCNamesByNumber () const
  {
    Array<string> names;
    if (dimension == 3)
      {
        int nbc = 0;
        bool named = false;
        for (const auto & fd : Regions<2>())
          {
            nbc = max(nbc, fd.BCProperty());
            if (fd.GetBCName() != "default") named = true;
          }
        if (!named) return names;
        names.SetSize(nbc);
        names = "default";
        Array<bool> done(nbc); done = false;
        for (const auto & fd : Regions<2>())
          if (fd.BCProperty() >= 1 && !done[fd.BCProperty()-1])
            { names[fd.BCProperty()-1] = fd.GetBCName(); done[fd.BCProperty()-1] = true; }
      }
    else if (dimension == 2)
      {
        for (const auto & ed : Regions<1>())
          names.Append(ed.GetName());
      }
    else
      for (auto & vd : Regions<0>())
        names.Append(vd.GetName());
    return names;
  }

  EdgeRegion & Mesh :: EnsureEdgeDescriptor (int nr)
  {
    while (Regions<1>().Size() < nr)
      Regions<1>().Append(EdgeRegion());
    return Regions<1>()[EdgeRegionIndex::FromNr1(nr)];
  }

  // cd2names of files and archives: edge names in 3D, vertex names in 2D
  void Mesh :: SetCD2NameCompat (int cd2nr, const string & name)
  {
    if (dimension == 3)
      EnsureEdgeDescriptor(cd2nr).SetName((name != "default" && !name.empty()) ? name : "default");
    else if (dimension == 2)
      SetCD2Name(cd2nr, name);
  }

  void Mesh :: SetCD2Name ( int cd2nr, const string & abcname )
  {
    if (dimension != 2) throw Exception("SetCD2Name names vertices of 2D meshes only");
    auto & verts = Regions<0>();
    while (verts.Size() < cd2nr)
      verts.Append(VertexRegion(cd2_default_name));
    verts[VertexRegionIndex::FromNr1(cd2nr)].SetName(abcname.empty() ? cd2_default_name : abcname);
  }

  string Mesh :: cd2_default_name = "default";
  string Mesh :: default_bc = "default";
  const string & Mesh :: GetCD2Name (int cd2nr) const
  {
    if (dimension == 2 && cd2nr >= 0 && cd2nr < Regions<0>().Size())
      return Regions<0>()[VertexRegionIndex::FromNr0(cd2nr)].GetName();
    return cd2_default_name;
  }

  void Mesh :: SetNCD3Names( int ncd3n )
  {
    Regions<0>() = RegionArray<0>(ncd3n);
  }

  void Mesh :: SetCD3Name ( int cd3nr, const string & abcname )
  {
    (*testout) << "setCD3Name on vertex " << cd3nr-1 << " to " << abcname << endl;
    auto & verts = Regions<0>();
    while (verts.Size() < cd3nr)
      verts.Append(VertexRegion());
    auto & vd = verts[VertexRegionIndex::FromNr1(cd3nr)];
    if (abcname != "default") vd.SetName(abcname); else vd.ResetName();
  }
  
  int Mesh :: AddCD3Name (const string & aname)
  {
    for (auto i : Regions<0>().Range())
      if (Regions<0>()[i].HasName() && Regions<0>()[i].GetName() == aname)
        return i.Nr0();
    return AddRegion(VertexRegion(aname)).Nr0();
  }
  
  string Mesh :: cd3_default_name = "default";
  static string defaultstring  = "default";
  const string & Mesh :: GetCD3Name (int cd3nr) const
  {
    if (cd3nr < 0 || cd3nr >= Regions<0>().Size())
      return defaultstring;
    return Regions<0>()[VertexRegionIndex::FromNr0(cd3nr)].GetName();
  }

  std::string_view Mesh :: GetRegionName (const Segment & el) const
  {
    if (HasEdgeDescriptor(el))
      return Regions<1>()[el.GetIndex()].GetName();
    return defaultmat_sv;
  }

  std::string_view Mesh :: GetRegionName (const Element2d & el) const
  {
    if (HasFaceDescriptor(el))
      return GetFaceDescriptor(el).GetBCName();
    return defaultmat_sv;
  }

  std::string_view Mesh :: GetRegionName (const ElementRef & el) const
  {
    return GetRegionName(3, el.GetIndex().Nr1());
  }
  

  void Mesh :: SetUserData(const char * id, Array<int> & data)
  {
    if(userdata_int.Used(id))
      delete userdata_int[id];

    Array<int> * newdata = new Array<int>(data);

    userdata_int.Set(id,newdata);      
  }
  bool Mesh :: GetUserData(const char * id, Array<int> & data, int shift) const
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
  void Mesh :: SetUserData(const char * id, Array<double> & data)
  {
    if(userdata_double.Used(id))
      delete userdata_double[id];

    Array<double> * newdata = new Array<double>(data);

    userdata_double.Set(id,newdata);      
  }
  bool Mesh :: GetUserData(const char * id, Array<double> & data, int shift) const
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
        << sizeof (netgen::Point<3>) << " + " << sizeof(POINTTYPE) << " = "
        << GetNP() * (sizeof (netgen::Point<3>) + sizeof(POINTTYPE)) << endl;

    ost << GetNSE() << " Surface elements, of size " 
        << sizeof (Element2d) << " = " 
        << GetNSE() * sizeof(Element2d) << endl;

    ost << GetNE() << " Volume elements, of size " 
        << sizeof (Element) << " = " 
        << GetNE() * sizeof(Element) << endl;

    // ost << "surfs on node:";
    // surfacesonnode.PrintMemInfo (cout);

    auto print_ht = [&ost] (const auto & ht, size_t elsize)
    {
      ost << "Hashtable: " << ht.Size()
          << " entries of size " << elsize
          << " = " << ht.Size() * elsize << " bytes."
          << " Used els: " << ht.UsedElements() << endl;
    };

    ost << "boundaryedges: ";
    if (boundaryedges)
      print_ht (*boundaryedges, sizeof(SortedPointIndices<2>) + sizeof(int));

    ost << "surfelementht: ";
    if (surfelementht)
      print_ht (*surfelementht, sizeof(SortedPointIndices<3>) + sizeof(SurfaceElementIndex));
  }

  shared_ptr<Mesh> Mesh :: Mirror ( netgen::Point<3> p_plane, Vec<3> n_plane )
  {
    Mesh & m = *this;
    auto nm_ = make_shared<Mesh>();
    Mesh & nm = *nm_;
    nm = m;

    netgen::Point<3> pmin, pmax;
    GetBox(pmin, pmax);
    auto v = pmax-pmin;
    double eps = v.Length()*1e-8;

    /*
    auto onPlane = [&] (const MeshPoint & p) -> bool
    {
      auto v = p_plane-p;
      auto l = v.Length();
      if(l<eps) return true;

      // auto ret = fabs(v*n_plane)/l;
      return fabs(v*n_plane) < eps;
    };
    */

    /*
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
    */

    Array<PointIndex, PointIndex> point_map(GetNP());
    Array<PointIndex, PointIndex> point_map1(GetNP());

    nm.Points().SetSize(0);
    
    for(auto pi : Range(points))
      {
        auto & p = m[pi];
        
        auto v = p_plane-p;
        auto l = v.Length();

        if(l < eps || fabs(v*n_plane)/l < eps)
          {
            auto npi = nm.AddPoint(p, p.GetLayer(), p.Type());
            point_map[pi] = npi;
            point_map1[pi] = npi;
          }
        else
          {
            auto new_point = p + 2*(v*n_plane)*n_plane;
            point_map1[pi] = nm.AddPoint(p, p.GetLayer(), p.Type());
            point_map[pi] = nm.AddPoint( new_point, p.GetLayer(), p.Type() );
          }
      }
    
    for (auto el : nm.VolumeElements())
      for(auto i : Range(el.GetNP()))
        el[i] = point_map1[el[i]];
    for(auto & el : nm.SurfaceElements())
      for(auto i : Range(el.GetNP()))
        el[i] = point_map1[el[i]];
    for(auto & el : nm.LineSegments())
      for(auto i : Range(el.GetNP()))
        el[i] = point_map1[el[i]];
    
    for (auto el : VolumeElements())
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
        {
          nel.Invert();
          nm.AddSurfaceElement(nel);
        }
    }

    for (auto ei : Range(LineSegments()))
    {
      auto & el = (*this)[ei];
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

    nm.ComputeNVertices();
    return nm_;
  }

  void AddFacesBetweenDomains(Mesh & mesh)
  {
    static Timer timer("AddFacesBetweenDomains"); RegionTimer rt(timer);
    auto & topo = mesh.GetTopology();
    auto p2el = mesh.CreatePoint2ElementTable();

    Array<size_t> els_per_domain(mesh.GetNDomains()+1);
    els_per_domain = 0;

    for(const auto & el : mesh.VolumeElements())
      els_per_domain[el.GetIndex().Nr1()]++;

    std::map<tuple<int,int>, int> doms_2_new_face;

    for(const auto & [facei, fd]: Enumerate(mesh.FaceDescriptors()))
    {
      auto dom0 = fd.DomainIn();
      auto dom1 = fd.DomainOut();
      if(dom0 > dom1)
        swap(dom0, dom1);

      doms_2_new_face[{dom0, dom1}] = facei+1;
    }

    for(auto dom : Range(1, 1+mesh.GetNDomains()))
    {
      if(els_per_domain[dom] == 0)
        continue;

      mesh.UpdateTopology();

      mesh.FindOpenElements(dom);
      for(const auto & openel : mesh.OpenElements())
      {
        std::set<ElementIndex> has_p1, has_p2, has_p3;
        for (auto ei: topo.GetVertexElements(openel[0]))
          has_p1.insert(ei);
        for (auto ei: topo.GetVertexElements(openel[1]))
          has_p2.insert(ei);
        for (auto ei: topo.GetVertexElements(openel[2]))
          has_p3.insert(ei);

        std::set<ElementIndex> has_p12, has_all;
        set_intersection(has_p1.begin(), has_p1.end(),
                         has_p2.begin(), has_p2.end(),
                         inserter(has_p12, has_p12.begin()));
        set_intersection(has_p12.begin(), has_p12.end(),
                         has_p3.begin(), has_p3.end(),
                         inserter(has_all, has_all.begin()));

        ArrayMem<ElementIndex, 5> els;
        for(auto ei : has_all)
          els.Append(ei);

        if(els.Size() == 2 && mesh[els[0]].GetIndex() != mesh[els[1]].GetIndex())
        {
          int dom0 = mesh[els[0]].GetIndex().Nr1();
          int dom1 = mesh[els[1]].GetIndex().Nr1();
          ElementIndex ei0 = els[0];
          if(dom0 > dom1)
          {
            Swap(dom0, dom1);
            ei0 = els[1];
          }

          if(dom1 == dom)
            continue;

          if(doms_2_new_face.count({dom0, dom1}) == 0)
          {
            auto fd = FaceRegion(-1, dom0, dom1, -1);
            auto new_si = mesh.GetNFD()+1;
            fd.SetBCProperty(new_si);
            auto new_face = mesh.AddFaceDescriptor(fd);
            mesh.SetBCName(new_si - 1, "default");
            doms_2_new_face[{dom0, dom1}] = new_face.Nr1();
          }
          for(auto face : topo.GetFaces(ei0)) {
            auto verts = topo.GetFaceVertices(face);
            if(verts.Contains(openel[0]) && verts.Contains(openel[1]) && verts.Contains(openel[2])) {
              Element2d sel(static_cast<int>(verts.Size()));
              sel.SetIndex(FaceRegionIndex::FromNr1(doms_2_new_face[{dom0, dom1}]));

              for(auto j : Range(verts.Size()))
                sel[j] = verts[j];
              auto normal = Cross(mesh[sel[1]]-mesh[sel[0]], mesh[sel[2]]-mesh[sel[0]]);
              Vec<3> surf_center = Vec<3>(Center(mesh[sel[0]] , mesh[sel[1]] , mesh[sel[2]]));
              Vec<3> center(0., 0., 0.);
              for(auto pi : mesh[ei0].PNums())
                center += Vec<3>(mesh[pi]);
              center *= 1.0/mesh[ei0].GetNP();
              if((normal * (center - surf_center)) < 0)
                sel.Invert();
              mesh.AddSurfaceElement(sel);
              break;
            }
          }
        }
      }
    }


  }

}
