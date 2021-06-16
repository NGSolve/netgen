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

      // mesh for one domain (contains all adjacent surface elments)
      Mesh mesh;

      // maps from lokal (domain) mesh to global mesh
      Array<PointIndex, PointIndex> pmap;

      // todo: store (mapped) identifications
      Array<INDEX_2> connected_pairs;

      MeshingParameters mp;

      unique_ptr<Meshing3> meshing;
  };

  // extract surface meshes belonging to individual domains
  Array<MeshingData> DivideMesh(const Mesh & mesh, const MeshingParameters & mp)
  {
      static Timer timer("DivideMesh"); RegionTimer rt(timer);

      Array<MeshingData> ret;
      auto num_domains = mesh.GetNDomains();
      ret.SetSize(num_domains);

      Array<Array<PointIndex, PointIndex>> ipmap;
      ipmap.SetSize(num_domains);
      auto dim = mesh.GetDimension();
      auto num_points = mesh.GetNP();
      auto num_facedescriptors = mesh.GetNFD();

      auto & identifications = mesh.GetIdentifications();

      for(auto i : Range(ret))
      {
          auto & md = ret[i];
          md.domain = i+1;

          if(mp.only3D_domain_nr && mp.only3D_domain_nr !=md.domain)
              continue;

          md.mp = mp;
          md.mp.maxh = min2 (mp.maxh, mesh.MaxHDomain(md.domain));

          auto & m = ret[i].mesh;

          m.SetLocalH(mesh.GetLocalH());

          ipmap[i].SetSize(num_points);
          ipmap[i] = PointIndex::INVALID;
          m.SetDimension( mesh.GetDimension() );

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

            if(mp.only3D_domain_nr && mp.only3D_domain_nr != dom)
                continue;

            auto & sels = ret[dom-1].mesh.SurfaceElements();
            for(auto pi : sel.PNums())
                ipmap[dom-1][pi] = 1;
            sels.Append(sel);
        }
      }

      // add used points to domain mesh, build point mapping
      for(auto i : Range(ret))
      {
          if(mp.only3D_domain_nr && mp.only3D_domain_nr != ret[i].domain)
              continue;

          auto & m = ret[i].mesh;
          auto & pmap = ret[i].pmap;

          for(auto pi : Range(ipmap[i]))
            if(ipmap[i][pi])
            {
              auto pi_new = m.AddPoint( mesh[pi] );
              ipmap[i][pi] = pi_new;
              pmap.Append( pi );
            }
      }

      NgArray<INDEX_2> connectednodes;
      for(auto i : Range(ret))
      {
          auto & imap = ipmap[i];
          if(mp.only3D_domain_nr && mp.only3D_domain_nr != ret[i].domain)
              continue;

          auto & m = ret[i].mesh;
          for (auto & sel : m.SurfaceElements())
            for(auto & pi : sel.PNums())
              pi = imap[pi];

           for (int nr = 1; nr <= identifications.GetMaxNr(); nr++)
             if (identifications.GetType(nr) != Identifications::PERIODIC)
               {
                 identifications.GetPairs (nr, connectednodes);
                 for (auto pair : connectednodes)
                 {
                     auto pi0 = pair[0];
                     auto pi1 = pair[1];
                     if(imap[pi0].IsValid() && imap[pi1].IsValid())
                         ret[i].connected_pairs.Append({imap[pi0], imap[pi1]});
                 }
               }
           // ret[i].mesh.Save("surface_"+ToString(i)+".vol");
      }
      return ret;
  }

  void CloseOpenQuads( MeshingData & md)
  {
    auto & mesh = md.mesh;
    auto domain = md.domain;
    MeshingParameters & mp = md.mp;

    if(mp.only3D_domain_nr && mp.only3D_domain_nr != domain)
        return;
     
    int oldne;
    if (multithread.terminate)
      return;
    
    mesh.CalcSurfacesOfNode();
    mesh.FindOpenElements();
    
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
           
           mpquad.giveuptol = 15;
           mpquad.baseelnp = 4;
           mpquad.starshapeclass = 1000;
           mpquad.check_impossible = qstep == 1;   // for prisms only (air domain in trafo)
           
           
           for (PointIndex pi : mesh.Points().Range())
             meshing.AddPoint (mesh[pi], pi);

           for (auto pair : md.connected_pairs)
               meshing.AddConnectedPair (pair);
           
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
           
           mesh.FindOpenElements();
         }
     }
   

   if (mesh.HasOpenQuads())
   {
      PrintSysError ("mesh has still open quads");
      throw NgException ("Stop meshing since too many attempts");
      // return MESHING3_GIVEUP;
   }
  }

  void PrepareForBlockFillLocalH(MeshingData & md)
  {
    static Timer t("PrepareForBlockFillLocalH"); RegionTimer rt(t);
    md.meshing = make_unique<Meshing3>(nullptr);

    auto & mesh = md.mesh;

    mesh.CalcSurfacesOfNode();
    mesh.FindOpenElements(md.domain);

    for (PointIndex pi : mesh.Points().Range())
       md.meshing->AddPoint (mesh[pi], pi);

    for (int i = 1; i <= mesh.GetNOpenElements(); i++)
       md.meshing->AddBoundaryElement (mesh.OpenElement(i));

    if (mesh.HasLocalHFunction())
        md.meshing->PrepareBlockFillLocalH(mesh, md.mp);
  }


  void MeshDomain( MeshingData & md)
  {
    auto & mesh = md.mesh;
    auto domain = md.domain;
    MeshingParameters & mp = md.mp;

    if(mp.only3D_domain_nr && mp.only3D_domain_nr != domain)
        return;




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
         GetOpenElements( mesh, domain )->Save("open_"+ToString(cntsteps)+".vol");
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

            mesh.FindOpenElements();	      
            PrintMessage (3, "Call remove problem");
            // mesh.Save("before_remove.vol");
            RemoveProblem (mesh, domain);
            // mesh.Save("after_remove.vol");
            mesh.FindOpenElements();
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
           // mesh.Save("output.vol");
           PrintError ("Surface mesh not consistent");
           throw NgException ("Stop meshing since surface mesh not consistent");
        }
     }
   // OptimizeVolume( md.mp, mesh );
  }

  void MeshDomain(Mesh & mesh3d, const MeshingParameters & c_mp, int k, const Identifications & identifications)
  {
    MeshingParameters mp = c_mp; // copy mp to change them here
    NgArray<INDEX_2> connectednodes;
     
    int oldne;
    int meshed;
    if(mp.only3D_domain_nr && mp.only3D_domain_nr !=k)
      return;
    if (multithread.terminate)
      return;
    
    PrintMessage (2, "");
    PrintMessage (1, "Meshing subdomain ", k, " of ", mesh3d.GetNDomains());
    (*testout) << "Meshing subdomain " << k << endl;
    
    mp.maxh = min2 (mp.maxh, mesh3d.MaxHDomain(k));

    mesh3d.CalcSurfacesOfNode();
    mesh3d.FindOpenElements(k);
    
    if (!mesh3d.GetNOpenElements())
      return;
    
    

    Box<3> domain_bbox( Box<3>::EMPTY_BOX ); 
    
    for (SurfaceElementIndex sei = 0; sei < mesh3d.GetNSE(); sei++)
      {
        const Element2d & el = mesh3d[sei];
        if (el.IsDeleted() ) continue;

        if (mesh3d.GetFaceDescriptor(el.GetIndex()).DomainIn() == k ||
            mesh3d.GetFaceDescriptor(el.GetIndex()).DomainOut() == k)
          
          for (int j = 0; j < el.GetNP(); j++)
            domain_bbox.Add (mesh3d[el[j]]);
      }
    domain_bbox.Increase (0.01 * domain_bbox.Diam());
    

    for (int qstep = 0; qstep <= 3; qstep++)
      // for (int qstep = 0; qstep <= 0; qstep++)  // for hex-filling
     {
       if (qstep == 0 && !mp.try_hexes) continue;
       
       // cout << "openquads = " << mesh3d.HasOpenQuads() << endl;
       if (mesh3d.HasOpenQuads())
         {
           string rulefile = ngdir;
           
           const char ** rulep = NULL;
           switch (qstep)
             {
             case 0:
               rulefile = "/Users/joachim/gitlab/netgen/rules/hexa.rls";
               rulep = hexrules;
               break;
             case 1:
               rulefile += "/rules/prisms2.rls";
               rulep = prismrules2;
               break;
             case 2: // connect pyramid to triangle
               rulefile += "/rules/pyramids2.rls";
               rulep = pyramidrules2;
               break;
             case 3: // connect to vis-a-vis point
               rulefile += "/rules/pyramids.rls";
               rulep = pyramidrules;
               break;
             }
           
           // Meshing3 meshing(rulefile);
           Meshing3 meshing(rulep); 
           
           MeshingParameters mpquad = mp;
           
           mpquad.giveuptol = 15;
           mpquad.baseelnp = 4;
           mpquad.starshapeclass = 1000;
           mpquad.check_impossible = qstep == 1;   // for prisms only (air domain in trafo)
           
           
           // for (PointIndex pi = mesh3d.Points().Begin(); pi < mesh3d.Points().End(); pi++)
           for (PointIndex pi : mesh3d.Points().Range())
             meshing.AddPoint (mesh3d[pi], pi);

           /*
           mesh3d.GetIdentifications().GetPairs (0, connectednodes);
           for (int i = 1; i <= connectednodes.Size(); i++)
             meshing.AddConnectedPair (connectednodes.Get(i));
           */
           // for (int nr = 1; nr <= identifications.GetMaxNr(); nr++)
           //   if (identifications.GetType(nr) != Identifications::PERIODIC)
           //     {
           //       identifications.GetPairs (nr, connectednodes);
           //       for (auto pair : connectednodes)
           //         meshing.AddConnectedPair (pair);
           //     }
           
           for (int i = 1; i <= mesh3d.GetNOpenElements(); i++)
             {
               Element2d hel = mesh3d.OpenElement(i);
               meshing.AddBoundaryElement (hel);
             }
           
           oldne = mesh3d.GetNE();
           
           meshing.GenerateMesh (mesh3d, mpquad);
           
           for (int i = oldne + 1; i <= mesh3d.GetNE(); i++)
             mesh3d.VolumeElement(i).SetIndex (k);
           
           (*testout) 
             << "mesh has " << mesh3d.GetNE() << " prism/pyramid elements" << endl;
           
           mesh3d.FindOpenElements(k);
         }
     }
   

   if (mesh3d.HasOpenQuads())
   {
      PrintSysError ("mesh has still open quads");
      throw NgException ("Stop meshing since too many attempts");
      // return MESHING3_GIVEUP;
   }


   if (mp.delaunay && mesh3d.GetNOpenElements())
   {
      Meshing3 meshing((const char**)NULL);

      mesh3d.FindOpenElements(k);

      /*
      for (PointIndex pi = mesh3d.Points().Begin(); pi < mesh3d.Points().End(); pi++)
         meshing.AddPoint (mesh3d[pi], pi);
      */
      for (PointIndex pi : mesh3d.Points().Range())
         meshing.AddPoint (mesh3d[pi], pi);

      for (int i = 1; i <= mesh3d.GetNOpenElements(); i++)
         meshing.AddBoundaryElement (mesh3d.OpenElement(i));

      oldne = mesh3d.GetNE();

      meshing.Delaunay (mesh3d, k, mp);

      for (int i = oldne + 1; i <= mesh3d.GetNE(); i++)
         mesh3d.VolumeElement(i).SetIndex (k);

      PrintMessage (3, mesh3d.GetNP(), " points, ",
         mesh3d.GetNE(), " elements");
   }


   int cntsteps = 0;
   if (mesh3d.GetNOpenElements())
     do
       {
         if (multithread.terminate)
           break;
         
         mesh3d.FindOpenElements(k);
         PrintMessage (5, mesh3d.GetNOpenElements(), " open faces");
         cntsteps++;

         if (cntsteps > mp.maxoutersteps) 
            throw NgException ("Stop meshing since too many attempts");

         string rulefile = ngdir + "/tetra.rls";
         PrintMessage (1, "start tetmeshing");

         //	  Meshing3 meshing(rulefile);
         Meshing3 meshing(tetrules);

         NgArray<int, PointIndex::BASE> glob2loc(mesh3d.GetNP());
         glob2loc = -1;

         // for (PointIndex pi = mesh3d.Points().Begin(); pi < mesh3d.Points().End(); pi++)
         for (PointIndex pi : mesh3d.Points().Range())
           if (domain_bbox.IsIn (mesh3d[pi]))
               glob2loc[pi] = 
               meshing.AddPoint (mesh3d[pi], pi);

         for (int i = 1; i <= mesh3d.GetNOpenElements(); i++)
         {
            Element2d hel = mesh3d.OpenElement(i);
            for (int j = 0; j < hel.GetNP(); j++)
               hel[j] = glob2loc[hel[j]];
            meshing.AddBoundaryElement (hel);
            // meshing.AddBoundaryElement (mesh3d.OpenElement(i));
         }

         oldne = mesh3d.GetNE();

         mp.giveuptol = 15 + 10 * cntsteps; 
         mp.sloppy = 5;
         meshing.GenerateMesh (mesh3d, mp);
         
         for (ElementIndex ei = oldne; ei < mesh3d.GetNE(); ei++)
            mesh3d[ei].SetIndex (k);
         

         mesh3d.CalcSurfacesOfNode();
         mesh3d.FindOpenElements(k);

         // teterrpow = 2;
         if (mesh3d.GetNOpenElements() != 0)
         {
            meshed = 0;
            PrintMessage (5, mesh3d.GetNOpenElements(), " open faces found");

            MeshOptimize3d optmesh(mp);

            const char * optstr = "mcmstmcmstmcmstmcm";
            for (size_t j = 1; j <= strlen(optstr); j++)
            {
               mesh3d.CalcSurfacesOfNode();
               mesh3d.FreeOpenElementsEnvironment(2);
               mesh3d.CalcSurfacesOfNode();

               switch (optstr[j-1])
               {
               case 'c': optmesh.CombineImprove(mesh3d, OPT_REST); break;
               case 'd': optmesh.SplitImprove(mesh3d, OPT_REST); break;
               case 's': optmesh.SwapImprove(mesh3d, OPT_REST); break;
               case 't': optmesh.SwapImprove2(mesh3d, OPT_REST); break;
               case 'm': mesh3d.ImproveMesh(mp, OPT_REST); break;
               }	  

            }

            mesh3d.FindOpenElements(k);	      
            PrintMessage (3, "Call remove problem");
            RemoveProblem (mesh3d, k);
            mesh3d.FindOpenElements(k);
         }
         else
           {
            meshed = 1;
            PrintMessage (1, "Success !");
           }
       }
     while (!meshed);
   
   PrintMessage (1, mesh3d.GetNP(), " points, ",
                 mesh3d.GetNE(), " elements");
     {
       if(mp.only3D_domain_nr && mp.only3D_domain_nr !=k)
	 return;
        PrintMessage (3, "Check subdomain ", k, " / ", mesh3d.GetNDomains());

        mesh3d.FindOpenElements(k);

        bool res = (mesh3d.CheckConsistentBoundary() != 0);
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
     for(auto & m_ : md)
     {
         auto first_new_pi = m_.pmap.Range().Next();
         auto & m = m_.mesh;
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

     if (mp.checkoverlappingboundary)
        if (mesh3d.CheckOverlappingBoundary())
           throw NgException ("Stop meshing since boundary mesh is overlapping");


     if(mesh3d.GetNDomains()==0)
         return MESHING3_OK;

     if (!mesh3d.HasLocalHFunction()) mesh3d.CalcLocalH(mp.grading);

     auto md = DivideMesh(mesh3d, mp);

     ParallelFor( md.Range(), [&](int i)
       {
         CloseOpenQuads( md[i] );
       });

     for(auto & md_ : md)
         PrepareForBlockFillLocalH(md_);

     ParallelFor( md.Range(), [&](int i)
       {
         MeshDomain(md[i]);
       });

     MergeMeshes(mesh3d, md);

     MeshQuality3d (mesh3d);

     return MESHING3_OK;
  }  


  // extern double teterrpow; 
  MESHING3_RESULT MeshVolume_ori (const MeshingParameters & mp, Mesh& mesh3d)
  {
    static Timer t("MeshVolume"); RegionTimer reg(t);

     if (!mesh3d.HasLocalHFunction()) mesh3d.CalcLocalH(mp.grading);

     mesh3d.Compress();

     //  mesh3d.PrintMemInfo (cout);

     if (mp.checkoverlappingboundary)
        if (mesh3d.CheckOverlappingBoundary())
           throw NgException ("Stop meshing since boundary mesh is overlapping");


     if(mesh3d.GetNDomains()==0)
         return MESHING3_OK;

     Array<Mesh> meshes(mesh3d.GetNDomains()-1);
     auto first_new_pi = mesh3d.Points().Range().Next();

     for(auto & m : meshes)
     {
         m = mesh3d;
         m.SetLocalH(mesh3d.GetLocalH());
     }

     ParallelFor(Range(1, mesh3d.GetNDomains()+1), [&](int k)
        {
          if(k==1)
             MeshDomain(mesh3d, mp, k, mesh3d.GetIdentifications());
          else
             MeshDomain(meshes[k-2], mp, k, mesh3d.GetIdentifications());
        });
     MergeMeshes(mesh3d, meshes, first_new_pi);

     MeshQuality3d (mesh3d);

     return MESHING3_OK;
  }  




  /*


  MESHING3_RESULT MeshVolumeOld (MeshingParameters & mp, Mesh& mesh3d)
  {
  int i, k, oldne;


  int meshed;
  int cntsteps; 


  PlotStatistics3d * pstat;
  if (globflags.GetNumFlag("silentflag", 1) <= 2)
  pstat = new XPlotStatistics3d;
  else
  pstat = new TerminalPlotStatistics3d;

  cntsteps = 0;
  do
  {
  cntsteps++;
  if (cntsteps > mp.maxoutersteps) 
  {
  return MESHING3_OUTERSTEPSEXCEEDED;
  }


  int noldp = mesh3d.GetNP();
      
      
  if ( (cntsteps == 1) && globflags.GetDefineFlag ("delaunay"))
  {
  cntsteps ++;

  mesh3d.CalcSurfacesOfNode();


  for (k = 1; k <= mesh3d.GetNDomains(); k++)
  {
  Meshing3 meshing(NULL, pstat);

  mesh3d.FindOpenElements(k);
	      
  for (i = 1; i <= noldp; i++)
  meshing.AddPoint (mesh3d.Point(i), i);
	      
  for (i = 1; i <= mesh3d.GetNOpenElements(); i++)
  {
  if (mesh3d.OpenElement(i).GetIndex() == k)
  meshing.AddBoundaryElement (mesh3d.OpenElement(i));
  }
	      
  oldne = mesh3d.GetNE();
  if (globflags.GetDefineFlag ("blockfill"))
  {
  if (!globflags.GetDefineFlag ("localh"))
  meshing.BlockFill 
  (mesh3d, mp.h * globflags.GetNumFlag ("relblockfillh", 1));
  else
  meshing.BlockFillLocalH (mesh3d);
  }
	      
  MeshingParameters mpd;
  meshing.Delaunay (mesh3d, mpd);

  for (i = oldne + 1; i <= mesh3d.GetNE(); i++)
  mesh3d.VolumeElement(i).SetIndex (k);
  }
  }

  noldp = mesh3d.GetNP();

  mesh3d.CalcSurfacesOfNode();
  mesh3d.FindOpenElements();
  for (k = 1; k <= mesh3d.GetNDomains(); k++)
  {
  Meshing3 meshing(globflags.GetStringFlag ("rules3d", NULL), pstat);
      
  Point3d pmin, pmax;
  mesh3d.GetBox (pmin, pmax, k);
	  
  rot.SetCenter (Center (pmin, pmax));

  for (i = 1; i <= noldp; i++)
  meshing.AddPoint (mesh3d.Point(i), i);

  for (i = 1; i <= mesh3d.GetNOpenElements(); i++)
  {
  if (mesh3d.OpenElement(i).GetIndex() == k)
  meshing.AddBoundaryElement (mesh3d.OpenElement(i));
  }

  oldne = mesh3d.GetNE();


  if ( (cntsteps == 1) && globflags.GetDefineFlag ("blockfill"))
  {
  if (!globflags.GetDefineFlag ("localh"))
  {
  meshing.BlockFill 
  (mesh3d, 
  mp.h * globflags.GetNumFlag ("relblockfillh", 1));
  }
  else
  {
  meshing.BlockFillLocalH (mesh3d);
  }
  }


  mp.giveuptol = int(globflags.GetNumFlag ("giveuptol", 15));

  meshing.GenerateMesh (mesh3d, mp);

  for (i = oldne + 1; i <= mesh3d.GetNE(); i++)
  mesh3d.VolumeElement(i).SetIndex (k);
  }



  mesh3d.CalcSurfacesOfNode();
  mesh3d.FindOpenElements();
      
  teterrpow = 2;
  if (mesh3d.GetNOpenElements() != 0)
  {
  meshed = 0;
  (*mycout) << "Open elements found, old" << endl;
  const char * optstr = "mcmcmcmcm";
  int j;
  for (j = 1; j <= strlen(optstr); j++)
  switch (optstr[j-1])
  {
  case 'c': mesh3d.CombineImprove(); break;
  case 'd': mesh3d.SplitImprove(); break;
  case 's': mesh3d.SwapImprove(); break;
  case 'm': mesh3d.ImproveMesh(2); break;
  }	  
	  
  (*mycout) << "Call remove" << endl;
  RemoveProblem (mesh3d);
  (*mycout) << "Problem removed" << endl;
  }
  else
  meshed = 1;
  }
  while (!meshed);

  MeshQuality3d (mesh3d);

  return MESHING3_OK;
  }  

  */




  /*
  MESHING3_RESULT MeshMixedVolume(MeshingParameters & mp, Mesh& mesh3d)
  {
    int i, j;
    MESHING3_RESULT res;
    Point3d pmin, pmax;

    mp.giveuptol = 10;
    mp.baseelnp = 4;
    mp.starshapeclass = 100;

    //  TerminalPlotStatistics3d pstat;
  
    Meshing3 meshing1("pyramids.rls");
    for (i = 1; i <= mesh3d.GetNP(); i++)
      meshing1.AddPoint (mesh3d.Point(i), i);

    mesh3d.FindOpenElements();
    for (i = 1; i <= mesh3d.GetNOpenElements(); i++)
      if (mesh3d.OpenElement(i).GetIndex() == 1)
	meshing1.AddBoundaryElement (mesh3d.OpenElement(i));

    res = meshing1.GenerateMesh (mesh3d, mp);

    mesh3d.GetBox (pmin, pmax);
    PrintMessage (1, "Mesh pyramids, res = ", res);
    if (res)
      exit (1);


    for (i = 1; i <= mesh3d.GetNE(); i++)
      mesh3d.VolumeElement(i).SetIndex (1);

    // do delaunay
  
    mp.baseelnp = 0;
    mp.starshapeclass = 5;

    Meshing3 meshing2(NULL);
    for (i = 1; i <= mesh3d.GetNP(); i++)
      meshing2.AddPoint (mesh3d.Point(i), i);
    
    mesh3d.FindOpenElements();
    for (i = 1; i <= mesh3d.GetNOpenElements(); i++)
      if (mesh3d.OpenElement(i).GetIndex() == 1)
	meshing2.AddBoundaryElement (mesh3d.OpenElement(i));

    MeshingParameters mpd;
    meshing2.Delaunay (mesh3d, mpd);

    for (i = 1; i <= mesh3d.GetNE(); i++)
      mesh3d.VolumeElement(i).SetIndex (1);


    mp.baseelnp = 0;
    mp.giveuptol = 10;

    for (int trials = 1; trials <= 50; trials++)
      {
	if (multithread.terminate)
	  return MESHING3_TERMINATE;

	Meshing3 meshing3("tetra.rls");
	for (i = 1; i <= mesh3d.GetNP(); i++)
	  meshing3.AddPoint (mesh3d.Point(i), i);
      
	mesh3d.FindOpenElements();
	for (i = 1; i <= mesh3d.GetNOpenElements(); i++)
	  if (mesh3d.OpenElement(i).GetIndex() == 1)
	    meshing3.AddBoundaryElement (mesh3d.OpenElement(i));
      
	if (trials > 1)
	  CheckSurfaceMesh2 (mesh3d);
	res = meshing3.GenerateMesh (mesh3d, mp);
      
	for (i = 1; i <= mesh3d.GetNE(); i++)
	  mesh3d.VolumeElement(i).SetIndex (1);

	if (res == 0) break;



	for (i = 1; i <= mesh3d.GetNE(); i++)
	  {
	    const Element & el = mesh3d.VolumeElement(i);
	    if (el.GetNP() != 4)
	      {
		for (j = 1; j <= el.GetNP(); j++)
		  mesh3d.AddLockedPoint (el.PNum(j));
	      }
	  }

	mesh3d.CalcSurfacesOfNode();
	mesh3d.FindOpenElements();

	MeshOptimize3d optmesh;

	teterrpow = 2;
	const char * optstr = "mcmcmcmcm";
	for (j = 1; j <= strlen(optstr); j++)
	  switch (optstr[j-1])
	    {
	    case 'c': optmesh.CombineImprove(mesh3d, OPT_REST); break;
	    case 'd': optmesh.SplitImprove(mesh3d); break;
	    case 's': optmesh.SwapImprove(mesh3d); break;
	    case 'm': mesh3d.ImproveMesh(); break;
	    }	  
	        
	RemoveProblem (mesh3d);
      }


    PrintMessage (1, "Meshing tets, res = ", res);
    if (res)
      {
	mesh3d.FindOpenElements();
	PrintSysError (1, "Open elements: ", mesh3d.GetNOpenElements());
	exit (1);
      }


  
    for (i = 1; i <= mesh3d.GetNE(); i++)
      {
	const Element & el = mesh3d.VolumeElement(i);
	if (el.GetNP() != 4)
	  {
	    for (j = 1; j <= el.GetNP(); j++)
	      mesh3d.AddLockedPoint (el.PNum(j));
	  }
      }
  
    mesh3d.CalcSurfacesOfNode();
    mesh3d.FindOpenElements();
  
    MeshOptimize3d optmesh;

    teterrpow = 2;
    const char * optstr = "mcmcmcmcm";
    for (j = 1; j <= strlen(optstr); j++)
      switch (optstr[j-1])
	{
	case 'c': optmesh.CombineImprove(mesh3d, OPT_REST); break;
	case 'd': optmesh.SplitImprove(mesh3d); break;
	case 's': optmesh.SwapImprove(mesh3d); break;
	case 'm': mesh3d.ImproveMesh(); break;
	}	  


    return MESHING3_OK;
  }
*/






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
