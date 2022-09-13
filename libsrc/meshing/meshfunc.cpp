#include <set>

#include <mystdlib.h>
#include "meshing.hpp"
#include "debugging.hpp"

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
      auto dim = mesh.GetDimension();
      auto num_points = mesh.GetNP();
      auto num_facedescriptors = mesh.GetNFD();

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
          ipmap[i] = PointIndex::INVALID;
          m.SetDimension( mesh.GetDimension() );
          m.SetGeometry( mesh.GetGeometry() );

          for(auto i : Range(1, num_facedescriptors+1))
              m.AddFaceDescriptor( mesh.GetFaceDescriptor(i) );
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
                ipmap[dom-1][pi] = 1;
            sels.Append(sel);
        }
      }

      // mark used points for already existing volume elements, add them (with wrong point numbers) to domain mesh
      for(const auto & el : mesh.VolumeElements())
      {
        auto dom = el.GetIndex();

        auto & els = ret[dom-1].mesh->VolumeElements();
        for(auto pi : el.PNums())
            ipmap[dom-1][pi] = 1;
        els.Append(el);
      }

      // mark locked/fixed points for each domain TODO: domain bounding box to add only relevant points?
      for(auto pi : mesh.LockedPoints())
        for(auto i : Range(ret))
          ipmap[i][pi] = 2;

      // add used points to domain mesh, build point mapping
      for(auto i : Range(ret))
      {
          auto & m = *ret[i].mesh;
          auto & pmap = ret[i].pmap;

          for(auto pi : Range(ipmap[i]))
            if(ipmap[i][pi])
            {
              const auto& mp = mesh[pi];
              auto pi_new = m.AddPoint( mp, mp.GetLayer(), mp.Type() );
              if(ipmap[i][pi] == 2)
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

      NgArray<int, PointIndex::BASE> map;
      std::set<std::tuple<int,int,int>> hex_faces;
      for(auto identnr : Range(1,nmax+1))
      {
          if(identifications.GetType(identnr) != Identifications::CLOSESURFACES)
              continue;

          identifications.GetMap(identnr, map);
          mesh.FindOpenElements(md.domain);

          for(auto & sel : mesh.OpenElements())
          {
              // For quads: check if this open element is already closed by a hex
              // this happends when we have identifications in two directions
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
              std::set<int> pis;
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
           
           for (int i = oldne + 1; i <= mesh.GetNE(); i++)
             mesh.VolumeElement(i).SetIndex (domain);
           
           (*testout) 
             << "mesh has " << mesh.GetNE() << " prism/pyramid elements" << endl;
           
           mesh.FindOpenElements(domain);
         }
     }
   

   if (mesh.HasOpenQuads())
   {
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

      for (int i = oldne + 1; i <= mesh.GetNE(); i++)
         mesh.VolumeElement(i).SetIndex (domain);

      PrintMessage (3, mesh.GetNP(), " points, ",
         mesh.GetNE(), " elements");
   }

   Box<3> domain_bbox( Box<3>::EMPTY_BOX ); 
   
   for (auto & sel : mesh.SurfaceElements())
     {
       if (sel.IsDeleted() ) continue;

       for (auto pi : sel.PNums())
         domain_bbox.Add (mesh[pi]);
     }
   domain_bbox.Increase (0.01 * domain_bbox.Diam());

   mesh.FindOpenElements(domain);

   int cntsteps = 0;
   int meshed;
   if (mesh.GetNOpenElements())
     do
       {
         if (multithread.terminate)
           break;
         
         mesh.FindOpenElements(domain);
         PrintMessage (5, mesh.GetNOpenElements(), " open faces");
         // GetOpenElements( mesh, domain )->Save("open_"+ToString(cntsteps)+".vol");
         cntsteps++;


         if (cntsteps > mp.maxoutersteps) 
            throw NgException ("Stop meshing since too many attempts");

         PrintMessage (1, "start tetmeshing");

         Meshing3 meshing(tetrules);

         Array<PointIndex, PointIndex> glob2loc(mesh.GetNP());
         glob2loc = PointIndex::INVALID;

         for (PointIndex pi : mesh.Points().Range())
           if (domain_bbox.IsIn (mesh[pi]))
               glob2loc[pi] = meshing.AddPoint (mesh[pi], pi);

         for (auto sel : mesh.OpenElements() )
         {
           for(auto & pi : sel.PNums())
               pi = glob2loc[pi];
           meshing.AddBoundaryElement (sel);
         }

         int oldne = mesh.GetNE();

         mp.giveuptol = 15 + 10 * cntsteps; 
         mp.sloppy = 5;
         meshing.GenerateMesh (mesh, mp);
         
         for (ElementIndex ei = oldne; ei < mesh.GetNE(); ei++)
            mesh[ei].SetIndex (domain);
         

         mesh.CalcSurfacesOfNode();
         mesh.FindOpenElements(domain);

         // teterrpow = 2;
         if (mesh.GetNOpenElements() != 0)
         {
            meshed = 0;
            PrintMessage (5, mesh.GetNOpenElements(), " open faces found");

            MeshOptimize3d optmesh(mp);

            const char * optstr = "mcmstmcmstmcmstmcm";
            for (size_t j = 1; j <= strlen(optstr); j++)
            {
               mesh.FindOpenElements();
               mesh.CalcSurfacesOfNode();
               mesh.FreeOpenElementsEnvironment(2);
               mesh.CalcSurfacesOfNode();

               switch (optstr[j-1])
               {
               case 'c': optmesh.CombineImprove(mesh, OPT_REST); break;
               case 'd': optmesh.SplitImprove(mesh, OPT_REST); break;
               case 's': optmesh.SwapImprove(mesh, OPT_REST); break;
               case 't': optmesh.SwapImprove2(mesh, OPT_REST); break;
               case 'm': mesh.ImproveMesh(mp, OPT_REST); break;
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
   
     {
        PrintMessage (3, "Check subdomain ", domain, " / ", mesh.GetNDomains());

        mesh.FindOpenElements(domain);

        bool res = (mesh.CheckConsistentBoundary() != 0);
        if (res)
        {
           PrintError ("Surface mesh not consistent");
           throw NgException ("Stop meshing since surface mesh not consistent");
        }
     }
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

     for(auto & m_ : md)
     {
         auto first_new_pi = m_.pmap.Range().Next();
         auto & m = *m_.mesh;
         Array<PointIndex, PointIndex> pmap(m.Points().Size());
         for(auto pi : Range(PointIndex(PointIndex::BASE), first_new_pi))
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
     }
  }

  void MergeMeshes( Mesh & mesh, FlatArray<Mesh> meshes, PointIndex first_new_pi )
  {
     // todo: optimize: count elements, alloc all memory, copy vol elements in parallel
     static Timer t("MergeMeshes"); RegionTimer rt(t);
     for(auto & m : meshes)
     {
         Array<PointIndex, PointIndex> pmap(m.Points().Size());
         for(auto pi : Range(PointIndex(PointIndex::BASE), first_new_pi))
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

     if (!mesh3d.HasLocalHFunction())
         mesh3d.CalcLocalH(mp.grading);

     auto md = DivideMesh(mesh3d, mp);

     try
       {
     ParallelFor( md.Range(), [&](int i)
       {
         if (mp.checkoverlappingboundary)
           if (md[i].mesh->CheckOverlappingBoundary())
             throw NgException ("Stop meshing since boundary mesh is overlapping");

         if(md[i].mesh->GetGeometry()->GetGeomType() == Mesh::GEOM_OCC)
            FillCloseSurface( md[i] );
         CloseOpenQuads( md[i] );
         MeshDomain(md[i]);
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
    RegionTaskManager rtm(mp.parallel_meshing ? mp.nthreads : 0);
    const char* savetask = multithread.task;
    multithread.task = "Optimize Volume";
    
    int i;

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
    for (auto i : Range(mp.optsteps3d))
      {
	if (multithread.terminate)
	  break;

	MeshOptimize3d optmesh(mp);

	// teterrpow = mp.opterrpow;
	// for (size_t j = 1; j <= strlen(mp.optimize3d); j++)
        for (auto j : Range(mp.optimize3d.size()))
	  {
            multithread.percent = 100.* (double(j)/mp.optimize3d.size() + i)/mp.optsteps3d;
	    if (multithread.terminate)
	      break;

	    switch (mp.optimize3d[j])
	      {
	      case 'c': optmesh.CombineImprove(mesh3d, OPT_REST); break;
	      case 'd': optmesh.SplitImprove(mesh3d); break;
	      case 'D': optmesh.SplitImprove2(mesh3d); break;
	      case 's': optmesh.SwapImprove(mesh3d); break;
                // case 'u': optmesh.SwapImproveSurface(mesh3d); break;
	      case 't': optmesh.SwapImprove2(mesh3d); break;
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
  
    multithread.task = savetask;
    return MESHING3_OK;
  }




  void RemoveIllegalElements (Mesh & mesh3d)
  {
    static Timer t("RemoveIllegalElements"); RegionTimer reg(t);
    
    int it = 10;
    int nillegal, oldn;

    PrintMessage (1, "Remove Illegal Elements");
    // return, if non-pure tet-mesh
    /*
      if (!mesh3d.PureTetMesh())
      return;
    */
    mesh3d.CalcSurfacesOfNode();

    nillegal = mesh3d.MarkIllegalElements();

    MeshingParameters dummymp;
    MeshOptimize3d optmesh(dummymp);
    while (nillegal && (it--) > 0)
      {
	if (multithread.terminate)
	  break;

	PrintMessage (5, nillegal, " illegal tets");
        optmesh.SplitImprove (mesh3d, OPT_LEGAL);

	mesh3d.MarkIllegalElements();  // test
	optmesh.SwapImprove (mesh3d, OPT_LEGAL);
	mesh3d.MarkIllegalElements();  // test
	optmesh.SwapImprove2 (mesh3d, OPT_LEGAL);

	oldn = nillegal;
	nillegal = mesh3d.MarkIllegalElements();

	if (oldn != nillegal)
	  it = 10;
      }
    PrintMessage (5, nillegal, " illegal tets");
  }
}
