#include <set>

#include <mystdlib.h>
#include "meshing.hpp"
#include "debugging.hpp"
#include "boundarylayer.hpp"

namespace netgen
{
  extern const char * tetrules[];
  // extern const char * tetrules2[];
  extern const char * prismrules2[];
  extern const char * pyramidrules[];
  extern const char * pyramidrules2[];
  extern const char * hexrules[];

  struct MeshingData
  {
      int domain;

      // mesh for one domain (contains all adjacent surface elements)
      unique_ptr<Mesh> mesh;

      // maps from local (domain) mesh to global mesh
      Array<PointIndex, PointIndex> pmap;

      // Array<INDEX_2> connected_pairs;

      MeshingParameters mp;

      unique_ptr<Meshing3> meshing;
  };

  // extract surface meshes belonging to individual domains
  Array<MeshingData> DivideMesh(Mesh & mesh, const MeshingParameters & mp)
  {
      static Timer timer("DivideMesh"); RegionTimer rt(timer);

      Array<MeshingData> ret;
      auto num_domains = mesh.GetNDomains();

      if(num_domains==1 || mp.only3D_domain_nr)
      {
          ret.SetSize(1);
          // no need to divide mesh, just fill in meshing data
          ret[0].domain = 1;
          if(mp.only3D_domain_nr)
              ret[0].domain = mp.only3D_domain_nr;

          ret[0].mesh.reset(&mesh); // careful, this unique_ptr must not delete &mesh! (it will be released in MergeMeshes after meshing)
          ret[0].mp = mp;
          return ret;
      }
      ret.SetSize(num_domains);

      Array<Array<PointIndex, PointIndex>> ipmap;
      ipmap.SetSize(num_domains);
      // auto dim = mesh.GetDimension();
      auto num_points = mesh.GetNP();
      auto num_facedescriptors = mesh.GetNFD();


      constexpr PointIndex state0 = IndexBASE<PointIndex>()-1; 
      constexpr PointIndex state1 = state0+1;
      constexpr PointIndex state2 = state0+2;
      
      for(auto i : Range(ret))
      {
          auto & md = ret[i];
          md.domain = i+1;

          md.mp = mp;
          md.mp.maxh = min2 (mp.maxh, mesh.MaxHDomain(md.domain));

          ret[i].mesh = make_unique<Mesh>();
          auto & m = *ret[i].mesh;

          m.SetLocalH(mesh.GetLocalH());

          ipmap[i].SetSize(num_points);
          ipmap[i] = state0;   // 0; // PointIndex::INVALID;
          m.SetDimension( mesh.GetDimension() );
          m.SetGeometry( mesh.GetGeometry() );

          for(auto i : Range(1, num_facedescriptors+1))
              m.AddFaceDescriptor( mesh.GetFaceDescriptor(i) );
      }

      // mark interior edge points
      for(const auto& seg : mesh.LineSegments())
        {
          if(seg.domin > 0 && seg.domin == seg.domout)
            {
              ipmap[seg.domin-1][seg[0]] = state1; // 1;
              ipmap[seg.domin-1][seg[1]] = state1; // 1;
            }
        }

      // mark used points for each domain, add surface elements (with wrong point numbers) to domain mesh
      for(const auto & sel : mesh.SurfaceElements())
      {
        const auto & fd = mesh.GetFaceDescriptor(sel.GetIndex());
        int dom_in  = fd.DomainIn();
        int dom_out = fd.DomainOut();

        for( auto dom : {dom_in, dom_out} )
        {
            if(dom==0)
              continue;

            auto & sels = ret[dom-1].mesh->SurfaceElements();
            for(auto pi : sel.PNums())
              ipmap[dom-1][pi] = state1;  // 1;
            sels.Append(sel);
        }
      }

      // mark used points for already existing volume elements, add them (with wrong point numbers) to domain mesh
      for(const auto & el : mesh.VolumeElements())
      {
        auto dom = el.GetIndex();
        
        auto & els = ret[dom-1].mesh->VolumeElements();
        for(auto pi : el.PNums())
          ipmap[dom-1][pi] = state1; // 1;
        els.Append(el);
      }

      // mark locked/fixed points for each domain TODO: domain bounding box to add only relevant points?
      for(auto pi : mesh.LockedPoints())
        for(auto i : Range(ret))
          ipmap[i][pi] = state2; // 2;

      // add used points to domain mesh, build point mapping
      for(auto i : Range(ret))
      {
          auto & m = *ret[i].mesh;
          auto & pmap = ret[i].pmap;

          for(auto pi : Range(ipmap[i]))
            if(ipmap[i][pi] != state0)
            {
              const auto& mp = mesh[pi];
              auto pi_new = m.AddPoint( mp, mp.GetLayer(), mp.Type() );
              if(ipmap[i][pi] == state2) // 2)
                mesh.AddLockedPoint(pi_new);
              ipmap[i][pi] = pi_new;
              pmap.Append( pi );
            }
      }

      // add segments
      for(auto i : Range(ret))
      {
          auto & imap = ipmap[i];
          auto & m = *ret[i].mesh;
          for(auto seg : mesh.LineSegments())
            if(imap[seg[0]].IsValid() && imap[seg[1]].IsValid())
              {
                  seg[0] = imap[seg[0]];
                  seg[1] = imap[seg[1]];
                  m.AddSegment(seg);
              }
      }

      auto & identifications = mesh.GetIdentifications();

      for(auto i : Range(ret))
      {
          auto & m = *ret[i].mesh;
          auto & imap = ipmap[i];
          auto nmax = identifications.GetMaxNr ();
          auto & m_ident = m.GetIdentifications();

          for (auto & sel : m.SurfaceElements())
            for(auto & pi : sel.PNums())
              pi = imap[pi];

          for (auto & el : m.VolumeElements())
            for(auto & pi : el.PNums())
              pi = imap[pi];

          for(auto n : Range(1,nmax+1))
          {
              NgArray<INDEX_2> pairs;
              identifications.GetPairs(n, pairs);

              for(auto pair : pairs)
              {
                  auto pi0 = imap[pair[0]];
                  auto pi1 = imap[pair[1]];
                  if(!pi0.IsValid() || !pi1.IsValid())
                      continue;

                  m_ident.Add(pi0, pi1, n);
              }
              m_ident.SetType( n, identifications.GetType(n) );
          }
      }
      return ret;
  }

  // Add between identified surface elements (only consider closesurface identifications)
  void FillCloseSurface( MeshingData & md)
  {
      static Timer timer("FillCloseSurface"); RegionTimer rtimer(timer);

      auto & mesh = *md.mesh;
      auto & identifications = mesh.GetIdentifications();
      auto nmax = identifications.GetMaxNr();

      bool have_closesurfaces = false;
      for(auto i : Range(1,nmax+1))
          if(identifications.GetType(i) == Identifications::CLOSESURFACES)
              have_closesurfaces = true;
      if(!have_closesurfaces)
          return;

      idmap_type map;
      std::set<std::tuple<PointIndex,PointIndex,PointIndex>> hex_faces;
      for(auto identnr : Range(1,nmax+1))
      {
          if(identifications.GetType(identnr) != Identifications::CLOSESURFACES)
              continue;

          identifications.GetMap(identnr, map);
          mesh.FindOpenElements(md.domain);

          for(auto & sel : mesh.OpenElements())
          {
              // For quads: check if this open element is already closed by a hex
              // this happens when we have identifications in two directions
              if(sel.GetNP() == 4)
              {
                  Element2d face = sel;
                  face.NormalizeNumbering();
                  if(hex_faces.count({face[0], face[1], face[2]}))
                      continue;
              }
              bool is_mapped = true;
              for(auto pi : sel.PNums())
                  if(!PointIndex(map[pi]).IsValid())
                      is_mapped = false;

              if(!is_mapped)
                  continue;

              // insert prism/hex
              auto np = sel.GetNP();
              Element el(2*np);
              std::set<PointIndex> pis;
              for(auto i : Range(np))
              {
                  el[i] = sel[i];
                  el[i+np] = map[sel[i]];
                  pis.insert(sel[i]);
                  pis.insert(map[sel[i]]);
              }

              // degenerate element (mapped element onto itself, might happen for surface elements connecting two identified faces)
              if(pis.size() < 2*np)
                  continue;

              // check if new element is inside current domain
              auto p0 = mesh[sel[0]];
              Vec<3> n = -Cross(mesh[sel[1]] - p0, mesh[sel[2]] - p0 );

              if(n*(mesh[el[np]]-p0) < 0.0)
                  continue;

              el.SetIndex(md.domain);
              mesh.AddVolumeElement(el);
              if(el.NP()==8)
              {
                  // remember all adjacent faces of the new hex (to skip corresponding openelements accordingly)
                  for(auto facei : Range(1,7))
                  {
                      Element2d face;
                      el.GetFace(facei, face);
                      face.NormalizeNumbering();
                      hex_faces.insert({face[0], face[1], face[2]});
                  }
              }
          }
      }
  }

  void CloseOpenQuads( MeshingData & md)
  {
    static Timer t("CloseOpenQuads"); RegionTimer rt(t);
    auto & mesh = *md.mesh;
    auto domain = md.domain;
    MeshingParameters & mp = md.mp;

    int oldne;
    if (multithread.terminate)
      return;
    
    mesh.CalcSurfacesOfNode();
    mesh.FindOpenElements(domain);
    
    if (!mesh.GetNOpenElements())
      return;

    for (int qstep = 0; qstep <= 3; qstep++)
     {
       if (qstep == 0 && !mp.try_hexes) continue;
       
       if (mesh.HasOpenQuads())
         {
           string rulefile = ngdir;
           
           const char ** rulep = NULL;
           switch (qstep)
             {
             case 0:
               rulep = hexrules;
               break;
             case 1:
               rulep = prismrules2;
               break;
             case 2: // connect pyramid to triangle
               rulep = pyramidrules2;
               break;
             case 3: // connect to vis-a-vis point
               rulep = pyramidrules;
               break;
             }
           
           Meshing3 meshing(rulep); 
           
           MeshingParameters mpquad = mp;
           
           mpquad.giveuptol = mp.giveuptolopenquads;
           mpquad.baseelnp = 4;
           mpquad.starshapeclass = 1000;
           mpquad.check_impossible = qstep == 1;   // for prisms only (air domain in trafo)
           
           
           for (PointIndex pi : mesh.Points().Range())
             meshing.AddPoint (mesh[pi], pi);

           NgArray<INDEX_2> connectednodes;
           for (int nr = 1; nr <= mesh.GetIdentifications().GetMaxNr(); nr++)
             if (mesh.GetIdentifications().GetType(nr) != Identifications::PERIODIC)
               {
                 mesh.GetIdentifications().GetPairs (nr, connectednodes);
                 for (auto pair : connectednodes)
                 {
                   meshing.AddConnectedPair (pair);
                   meshing.AddConnectedPair ({pair[1], pair[0]});
                 }
               }
           
           for (int i = 1; i <= mesh.GetNOpenElements(); i++)
             {
               Element2d hel = mesh.OpenElement(i);
               meshing.AddBoundaryElement (hel);
             }
           
           oldne = mesh.GetNE();
           
           meshing.GenerateMesh (mesh, mpquad);
           
           // for (int i = oldne + 1; i <= mesh.GetNE(); i++)
           for (ElementIndex i : mesh.VolumeElements().Range().Modify(oldne, 0))
             mesh.VolumeElement(i).SetIndex (domain);
           
           (*testout) 
             << "mesh has " << mesh.GetNE() << " prism/pyramid elements" << endl;
           
           mesh.FindOpenElements(domain);
         }
     }
   

   if (mesh.HasOpenQuads())
   {
      if(debugparam.write_mesh_on_error) {
        md.mesh->Save("open_quads_starting_mesh_"+ToString(md.domain)+".vol.gz");
        GetOpenElements(*md.mesh, md.domain)->Save("open_quads_rest_" + ToString(md.domain)+".vol.gz");
        GetOpenElements(*md.mesh, md.domain, true)->Save("open_quads_rest_" + ToString(md.domain)+"_only_quads.vol.gz");
      }
      PrintSysError ("mesh has still open quads");
      throw NgException ("Stop meshing since too many attempts");
      // return MESHING3_GIVEUP;
   }
  }


  void MeshDomain( MeshingData & md)
  {
    auto & mesh = *md.mesh;
    auto domain = md.domain;
    MeshingParameters & mp = md.mp;

    mesh.CalcSurfacesOfNode();
    mesh.FindOpenElements(md.domain);

    md.meshing = make_unique<Meshing3>(nullptr);
    for (PointIndex pi : mesh.Points().Range())
       md.meshing->AddPoint (mesh[pi], pi);

    for (int i = 1; i <= mesh.GetNOpenElements(); i++)
       md.meshing->AddBoundaryElement (mesh.OpenElement(i));

    if (mp.delaunay && mesh.GetNOpenElements())
    {
      int oldne = mesh.GetNE();

      md.meshing->Delaunay (mesh, domain, mp);

      // for (int i = oldne + 1; i <= mesh.GetNE(); i++)
      for (ElementIndex i : mesh.VolumeElements().Range().Modify(oldne, 0))
         mesh.VolumeElement(i).SetIndex (domain);

      PrintMessage (3, mesh.GetNP(), " points, ",
         mesh.GetNE(), " elements");
      mesh.FindOpenElements(domain);
    }

    Box<3> domain_bbox( Box<3>::EMPTY_BOX ); 
   
    for (auto & sel : mesh.SurfaceElements())
     {
       if (sel.IsDeleted() ) continue;

       for (auto pi : sel.PNums())
         domain_bbox.Add (mesh[pi]);
     }
    domain_bbox.Increase (0.01 * domain_bbox.Diam());

    int cntsteps = 0;
    int meshed;
    if (mesh.GetNOpenElements())
     do
       {
         if (multithread.terminate)
           break;
         
         mesh.FindOpenElements(domain);
         PrintMessage (5, mesh.GetNOpenElements(), " open faces");
         // GetOpenElements( mesh, domain )->Save("open_"+ToString(domain)+"_"+ToString(cntsteps)+".vol");
         cntsteps++;


         if (cntsteps > mp.maxoutersteps) 
         {
           if(debugparam.write_mesh_on_error)
           {
             md.mesh->Save("meshing_error_domain_"+ToString(md.domain)+".vol.gz");
             if(mesh.GetNOpenElements())
               GetOpenElements(*md.mesh, md.domain)->Save("meshing_error_rest_" + ToString(md.domain)+".vol.gz");
           }
           throw NgException ("Stop meshing since too many attempts in domain " + ToString(md.domain));
         }

         PrintMessage (1, "start tetmeshing");

         Meshing3 meshing(tetrules);

         Array<PointIndex, PointIndex> glob2loc(mesh.GetNP());
         glob2loc = PointIndex::INVALID;

         for (PointIndex pi : mesh.Points().Range())
           if (domain_bbox.IsIn (mesh[pi]))
             glob2loc[pi] = meshing.AddPoint (mesh[pi], pi);

         for (auto sel : mesh.OpenElements())
           {
             for(auto & pi : sel.PNums())
               pi = glob2loc[pi];
             meshing.AddBoundaryElement (sel);
           }

         int oldne = mesh.GetNE();

         mp.giveuptol = 15 + 10 * cntsteps; 
         mp.sloppy = 5;
         meshing.GenerateMesh (mesh, mp);
         
         for (auto & el : mesh.VolumeElements().Range(oldne, END))
           el.SetIndex (domain);
         

         mesh.CalcSurfacesOfNode();
         mesh.FindOpenElements(domain);
         CheckMesh(mesh, MESHCONST_OPTVOLUME, __FILE__, __LINE__);

         // teterrpow = 2;
         if (mesh.GetNOpenElements() != 0)
         {
            meshed = 0;
            PrintMessage (5, mesh.GetNOpenElements(), " open faces found");

            MeshOptimize3d optmesh(mesh, mp, OPT_REST);

            const char * optstr = "mcmstmcmstmcmstmcm";
            for (size_t j = 1; j <= strlen(optstr); j++)
            {
               CheckMesh(mesh, MESHCONST_OPTVOLUME, __FILE__, __LINE__);
               mesh.FindOpenElements();
               mesh.CalcSurfacesOfNode();
               mesh.FreeOpenElementsEnvironment(2);
               mesh.CalcSurfacesOfNode();

               switch (optstr[j-1])
               {
               case 'c': optmesh.CombineImprove(); break;
               case 'd': optmesh.SplitImprove(); break;
               case 's': optmesh.SwapImprove(); break;
               case 't': optmesh.SwapImprove2(); break;
               case 'm': optmesh.ImproveMesh(); break;
               }	  

            }

            mesh.FindOpenElements(domain);
            PrintMessage (3, "Call remove problem");
            RemoveProblem (mesh, domain);
            mesh.FindOpenElements(domain);
         }
         else
           {
            meshed = 1;
            PrintMessage (1, "Success !");
           }
       }
    while (!meshed);

    PrintMessage (3, "Check subdomain ", domain, " / ", mesh.GetNDomains());

    mesh.FindOpenElements(domain);

    bool res = (mesh.CheckConsistentBoundary() != 0);
    if (res)
    {
      if(debugparam.write_mesh_on_error)
        md.mesh->Save("inconsistent_surface_domain_"+ToString(md.domain)+".vol.gz");
      PrintError ("Surface mesh not consistent");
      throw NgException ("Stop meshing since surface mesh not consistent");
    }
    RemoveIllegalElements (mesh, domain);
    ConformToFreeSegments (mesh, domain);
  }

  void MergeMeshes( Mesh & mesh, Array<MeshingData> & md )
  {
     // todo: optimize: count elements, alloc all memory, copy vol elements in parallel
     static Timer t("MergeMeshes"); RegionTimer rt(t);
     if(md.Size()==1)
     {
         // assume that mesh was never divided, no need to do anything
         if(&mesh != md[0].mesh.get())
             throw Exception("Illegal Mesh pointer in MeshingData");

         md[0].mesh.release();
         return;
     }

     mesh.VolumeElements().DeleteAll();
     mesh.GetIdentifications().GetIdentifiedPoints().DeleteData();

     for(auto & m_ : md)
     {
         auto first_new_pi = m_.pmap.Range().Next();
         auto & m = *m_.mesh;
         Array<PointIndex, PointIndex> pmap(m.Points().Size());
         for(auto pi : Range(IndexBASE<PointIndex>(), first_new_pi))
             pmap[pi] = m_.pmap[pi];

         for (auto pi : Range(first_new_pi, m.Points().Range().Next()))
             pmap[pi] = mesh.AddPoint(m[pi]);


         for ( auto el : m.VolumeElements() )
         {
             for (auto i : Range(el.GetNP()))
                 el[i] = pmap[el[i]];
             el.SetIndex(m_.domain);
             mesh.AddVolumeElement(el);
         }
         // for(const auto& [p1p2, dummy] : m.GetIdentifications().GetIdentifiedPoints())
         // mesh.GetIdentifications().Add(pmap[p1p2[0]], pmap[p1p2[1]], p1p2[2]);
         for(const auto& [p1p2, dummy] : m.GetIdentifications().GetIdentifiedPoints())         
           mesh.GetIdentifications().Add( pmap[ get<0>(p1p2)[0] ], pmap[ get<0>(p1p2)[1]] , get<1>(p1p2) );
         for(auto i : Range(m.GetIdentifications().GetMaxNr()))
           {
             mesh.GetIdentifications().SetType(i+1, m.GetIdentifications().GetType(i+1));
             if(auto name = m.GetIdentifications().GetName(i+1); name != "")
               mesh.GetIdentifications().SetName(i+1, name);
           }
     }
  }

  void MergeMeshes( Mesh & mesh, FlatArray<Mesh> meshes, PointIndex first_new_pi )
  {
     // todo: optimize: count elements, alloc all memory, copy vol elements in parallel
     static Timer t("MergeMeshes"); RegionTimer rt(t);
     for(auto & m : meshes)
     {
         Array<PointIndex, PointIndex> pmap(m.Points().Size());
         for(auto pi : Range(IndexBASE<PointIndex>(), first_new_pi))
             pmap[pi] = pi;

         for (auto pi : Range(first_new_pi, m.Points().Range().Next()))
             pmap[pi] = mesh.AddPoint(m[pi]);


         for ( auto el : m.VolumeElements() )
         {
             for (auto i : Range(el.GetNP()))
                 el[i] = pmap[el[i]];
             mesh.AddVolumeElement(el);
         }
     }
  }

  // extern double teterrpow; 
  MESHING3_RESULT MeshVolume (const MeshingParameters & mp, Mesh& mesh3d)
  {
    static Timer t("MeshVolume"); RegionTimer reg(t);

     mesh3d.Compress();

     if(mesh3d.GetNDomains()==0)
         return MESHING3_OK;

     auto geo = mesh3d.GetGeometry();
     for (auto i : Range(std::min(geo->GetNSolids(), (size_t)mesh3d.GetNDomains())))
       if (auto name = geo->GetSolid(i).properties.name)
         mesh3d.SetMaterial (i+1, *name);

     for (auto bl : mp.boundary_layers)
       GenerateBoundaryLayer(mesh3d, bl);

     if (!mesh3d.HasLocalHFunction())
         mesh3d.CalcLocalH(mp.grading);

     auto md = DivideMesh(mesh3d, mp);

     try
       {
     ParallelFor( md.Range(), [&](int i)
       {
         try {
           if (mp.checkoverlappingboundary)
             if (md[i].mesh->CheckOverlappingBoundary())
             {
               if(debugparam.write_mesh_on_error)
                 md[i].mesh->Save("overlapping_mesh_domain_"+ToString(md[i].domain)+".vol.gz");
               throw NgException ("Stop meshing since boundary mesh is overlapping");
             }

           if(md[i].mesh->GetGeometry()->GetGeomType() == Mesh::GEOM_OCC)
              FillCloseSurface( md[i] );
           CloseOpenQuads( md[i] );
           MeshDomain(md[i]);
         }
         catch (const Exception & e) {
           if(debugparam.write_mesh_on_error)
             md[i].mesh->Save("meshing_error_domain_"+ToString(md[i].domain)+".vol.gz");
           cerr << "Meshing of domain " << i+1 << " failed with error: " << e.what() << endl;
           throw e;
         }
       }, md.Size());
       }
     catch(...)
       {
         MergeMeshes(mesh3d, md);
         return MESHING3_GIVEUP;
       }

     MergeMeshes(mesh3d, md);

     MeshQuality3d (mesh3d);

     return MESHING3_OK;
  }  


  MESHING3_RESULT OptimizeVolume (const MeshingParameters & mp, 
				  Mesh & mesh3d)
    //				  const CSGeometry * geometry)
  {
    static Timer t("OptimizeVolume"); RegionTimer reg(t);
  #ifndef EMSCRIPTEN
    RegionTaskManager rtm(mp.parallel_meshing ? mp.nthreads : 0);
  #endif // EMSCRIPTEN
    const char* savetask = multithread.task;
    multithread.task = "Optimize Volume";
    
    // int i;

    PrintMessage (1, "Volume Optimization");

    /*
      if (!mesh3d.PureTetMesh())
      return MESHING3_OK;
    */

    // (*mycout) << "optstring = " << mp.optimize3d << endl;
    /*
      const char * optstr = globflags.GetStringFlag ("optimize3d", "cmh");
      int optsteps = int (globflags.GetNumFlag ("optsteps3d", 2));
    */

    mesh3d.CalcSurfacesOfNode();
    CheckMesh(mesh3d, MESHCONST_OPTVOLUME, __FILE__, __LINE__);

    MeshOptimize3d optmesh(mesh3d, mp);

    // optimize only bad elements first
    optmesh.SetMinBadness(1000.);
    bool do_split = mp.optimize3d.find('d') != string::npos;
    bool do_swap = mp.optimize3d.find('s') != string::npos;
    bool do_swap2 = mp.optimize3d.find('t') != string::npos;
    for([[maybe_unused]] auto i : Range(mp.optsteps3d))
      {
        CheckMesh(mesh3d, MESHCONST_OPTVOLUME, __FILE__, __LINE__);
        auto [total_badness, max_badness, bad_els] = optmesh.UpdateBadness();
        if(bad_els==0) break;
        if(do_split) optmesh.SplitImprove();
        if(do_swap) optmesh.SwapImprove();
        if(do_swap2) optmesh.SwapImprove2();
      }

    // Now optimize all elements
    optmesh.SetMinBadness(0);

    CheckMesh(mesh3d, MESHCONST_OPTVOLUME, __FILE__, __LINE__);

    for (auto i : Range(mp.optsteps3d))
      {
	if (multithread.terminate)
	  break;

	// teterrpow = mp.opterrpow;
	// for (size_t j = 1; j <= strlen(mp.optimize3d); j++)
        for (auto j : Range(mp.optimize3d.size()))
	  {
            multithread.percent = 100.* (double(j)/mp.optimize3d.size() + i)/mp.optsteps3d;
	    if (multithread.terminate)
	      break;

      CheckMesh(mesh3d, MESHCONST_OPTVOLUME, __FILE__, __LINE__);

	    switch (mp.optimize3d[j])
	      {
	      case 'c': 
          optmesh.SetGoal(OPT_REST);
          optmesh.CombineImprove();
          optmesh.SetGoal(OPT_QUALITY);
          break;
	      case 'd': optmesh.SplitImprove(); break;
	      case 'D': optmesh.SplitImprove2(); break;
	      case 's': optmesh.SwapImprove(); break;
                // case 'u': optmesh.SwapImproveSurface(mesh3d); break;
	      case 't': optmesh.SwapImprove2(); break;
#ifdef SOLIDGEOM
	      case 'm': mesh3d.ImproveMesh(*geometry); break;
	      case 'M': mesh3d.ImproveMesh(*geometry); break;
#else
	      case 'm': mesh3d.ImproveMesh(mp); break;
	      case 'M': mesh3d.ImproveMesh(mp); break;
#endif
	      case 'j': mesh3d.ImproveMeshJacobian(mp); break;
	      }
	  }
	// mesh3d.mglevels = 1;
	MeshQuality3d (mesh3d);
      }

    CheckMesh(mesh3d, MESHCONST_OPTVOLUME, __FILE__, __LINE__);
  
    multithread.task = savetask;
    return MESHING3_OK;
  }


  void ConformToFreeSegments (Mesh & mesh, int domain)
  {
    auto geo = mesh.GetGeometry();
    if(!geo) return;
    auto n_solids = geo->GetNSolids();
    if(n_solids < domain) return;
    if(geo->GetSolid(domain-1).free_edges.Size() == 0)
      return;

    Array<SegmentIndex> free_segs;
    for (auto segi : Range(mesh.LineSegments()))
      if(mesh[segi].domin == domain && mesh[segi].domout == domain)
        free_segs.Append(segi);

    auto get_nonconforming = [&] (const auto & p2el) {
      Array<size_t> nonconforming;

      for (auto segi : free_segs) {
        auto seg = mesh[segi];

        auto has_p0 = p2el[seg[0]];
        bool has_both = false;

        for(auto ei : has_p0) {
          if(mesh[ei].PNums().Contains(seg[1]))
            has_both = true;
        }

        if(!has_both)
          nonconforming.Append(segi);
      }
      return nonconforming;
    };

    auto split_segment = [&] (SegmentIndex segi, const auto & p2el) {
      // Todo: handle curved segments correctly
      auto seg = mesh[segi];
      auto p_new = Center(mesh[seg[0]], mesh[seg[1]]);
      auto pi_new = mesh.AddPoint(p_new);
      auto seg_new0 = seg;
      auto seg_new1 = seg;
      seg_new0[1] = pi_new;
      seg_new1[0] = pi_new;

      mesh[segi][0] = PointIndex::INVALID;
      mesh.AddSegment(seg_new0);
      mesh.AddSegment(seg_new1);

      double lam[3];
      ElementIndex ei = mesh.GetElementOfPoint(p_new, lam, false, domain);
      if(ei == 0) {
        PrintMessage(1, "Could not find volume element with new point");
        return;
      }
      ei -= 1;

      // split tet into 4 new tests, with new point inside
      auto el = mesh[ei];
      if(el.GetNP() != 4) {
        PrintMessage(1, "Only tet elements are supported to split around free segments");
        return;
      }

      if(el.IsDeleted()) {
        PrintMessage(1,"Element to split is already deleted");
        return;
      }

      int pmap[4][4] = {
        {0,1,2,4},
        {1,3,2,4},
        {0,2,3,4},
        {0,3,1,4}
      };

      PointIndex pis[5] = {el[0], el[1], el[2], el[3], pi_new};

      for (auto i : Range(4)) {
        Element el_new;
        el_new = el;
        for (auto j : Range(4))
          el_new[j] = pis[pmap[i][j]];
        mesh.AddVolumeElement(el_new);
      }
      mesh[ei].Delete();
    };

    size_t last_num_bad_segs = -1;
    for ([[maybe_unused]] auto i : Range(10)) {
      auto p2el = mesh.CreatePoint2ElementTable();

      auto bad_segs = get_nonconforming(p2el);
      auto num_bad_segs = bad_segs.Size();

      if(num_bad_segs == 0)
        return;

      PrintMessage(3, "Non-conforming free segments in domain ", domain, ": ", num_bad_segs);

      if(i>=5 || num_bad_segs != last_num_bad_segs) {
        for(auto i : bad_segs)
          split_segment(i, p2el);
        mesh.Compress();
      }

      MeshingParameters dummymp;
      MeshOptimize3d optmesh(mesh, dummymp, OPT_CONFORM);

      for ([[maybe_unused]] auto i : Range(3)) {
        CheckMesh(mesh, MESHCONST_OPTVOLUME, __FILE__, __LINE__);
        optmesh.ImproveMesh();
        CheckMesh(mesh, MESHCONST_OPTVOLUME, __FILE__, __LINE__);
        optmesh.SwapImprove2 ();
        CheckMesh(mesh, MESHCONST_OPTVOLUME, __FILE__, __LINE__);
        optmesh.ImproveMesh();
        CheckMesh(mesh, MESHCONST_OPTVOLUME, __FILE__, __LINE__);
        optmesh.SwapImprove();
        CheckMesh(mesh, MESHCONST_OPTVOLUME, __FILE__, __LINE__);
        optmesh.ImproveMesh();
        CheckMesh(mesh, MESHCONST_OPTVOLUME, __FILE__, __LINE__);
        optmesh.CombineImprove();
        CheckMesh(mesh, MESHCONST_OPTVOLUME, __FILE__, __LINE__);
      }
      last_num_bad_segs = num_bad_segs;
    }

    auto p2el = mesh.CreatePoint2ElementTable();
    auto bad_segs = get_nonconforming(p2el);

    if(bad_segs.Size() > 0) {
      auto bad_seg = mesh[free_segs[bad_segs[0]]];
      if(debugparam.write_mesh_on_error)
        mesh.Save("free_segment_not_conformed_dom_"+ToString(domain)+"_seg_"+ToString(bad_seg[0])+"_"+ToString(bad_seg[1])+".vol.gz");
      throw Exception("Segment not resolved in volume mesh in domain " + ToString(domain)+ ", seg: " + ToString(bad_seg));
    }
  }


  void RemoveIllegalElements (Mesh & mesh3d, int domain)
  {
    static Timer t("RemoveIllegalElements"); RegionTimer reg(t);
    
    // return, if non-pure tet-mesh
    /*
      if (!mesh3d.PureTetMesh())
      return;
    */
    mesh3d.CalcSurfacesOfNode();

    int nillegal = mesh3d.MarkIllegalElements(domain);
    if(nillegal)
      PrintMessage (1, "Remove Illegal Elements");

    int oldn = nillegal;
    int nillegal_min = nillegal;

    MeshingParameters dummymp;
    MeshOptimize3d optmesh(mesh3d, dummymp, OPT_LEGAL);
    int it = 10;
    while (nillegal && (it--) > 0)
      {
	if (multithread.terminate)
	  break;

	PrintMessage (5, nillegal, " illegal tets");
        optmesh.SplitImprove ();

	mesh3d.MarkIllegalElements();  // test
	optmesh.SwapImprove ();
	mesh3d.MarkIllegalElements();  // test
	optmesh.SwapImprove2 ();

	oldn = nillegal;
	nillegal = mesh3d.MarkIllegalElements();
        nillegal_min = min(nillegal_min, nillegal);
        if(nillegal > nillegal_min)
          break;

	if (oldn != nillegal)
	  it = 10;
      }
    PrintMessage (5, nillegal, " illegal tets");
  }
}
