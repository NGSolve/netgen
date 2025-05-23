#include <mystdlib.h>

#include "meshing.hpp"
#include <opti.hpp>

namespace netgen
{

  class ImprovementRule
  {
  public:
    NgArray<Element2d> oldels;
    NgArray<Element2d> newels;
    NgArray<PointIndices<2>> deledges;
    Array<int,PointIndex> incelsonnode;
    Array<int,PointIndex> reused;
    int bonus;
    int onp;
  };


  void MeshOptimize2d :: GenericImprove ()
  {
    static Timer timer("MeshOptimize2d::GenericImprove"); RegionTimer reg(timer);
    if (!faceindex)
      {
	if (writestatus)
	  PrintMessage (3, "Generic Improve");

	for (faceindex = 1; faceindex <= mesh.GetNFD(); faceindex++)
	  GenericImprove ();
      
	faceindex = 0;
        return;
      }

    // int j, k, l, ri;
    int np = mesh.GetNP();
    int ne = mesh.GetNSE();
    //    SurfaceElementIndex sei;

    
//     for (SurfaceElementIndex sei = 0; sei < ne; sei++)
//       {
// 	const Element2d & el = mesh[sei];
// 	(*testout) << "element " << sei << ": " <<flush;
// 	for(int j=0; j<el.GetNP(); j++)
// 	  (*testout) << el[j] << " " << flush;
// 	(*testout) << "IsDeleted() " << el.IsDeleted()<< endl;
//       }

    bool ok;
    int olddef, newdef;

    NgArray<ImprovementRule*> rules;
    NgArray<SurfaceElementIndex> elmap;
    NgArray<int> elrot;
    Array<PointIndex,PointIndex> pmap;
    Array<PointGeomInfo,PointIndex> pgi;

    int surfnr = mesh.GetFaceDescriptor (faceindex).SurfNr();
  
    ImprovementRule * r1;

    PointIndex p1 = IndexBASE<PointIndex>();
    PointIndex p2 = p1+1;
    PointIndex p3 = p2+1;
    PointIndex p4 = p3+1;
    PointIndex p5 = p4+1;
    PointIndex p6 = p5+1;
    PointIndex p7 = p6+1;
    
    
    // 2 triangles to quad
    r1 = new ImprovementRule;
    r1->oldels.Append (Element2d (p1, p2, p3));
    r1->oldels.Append (Element2d (p3, p2, p4));
    r1->newels.Append (Element2d (p1, p2, p4, p3));
    r1->deledges.Append ( { p2, p3 } );
    r1->onp = 4;
    r1->bonus = 2;
    rules.Append (r1);

    // 2 quad to 1 quad
    r1 = new ImprovementRule;
    r1->oldels.Append (Element2d (p1, p2, p3, p4));
    r1->oldels.Append (Element2d (p4, p3, p2, p5));
    r1->newels.Append (Element2d (p1, p2, p5, p4));
    r1->deledges.Append ( { p2, p3 } );
    r1->deledges.Append ( { p3, p4 } );
    r1->onp = 5;
    r1->bonus = 0;
    rules.Append (r1);

    // swap quads
    r1 = new ImprovementRule;
    r1->oldels.Append (Element2d (p1, p2, p3, p4));
    r1->oldels.Append (Element2d (p3, p2, p5, p6));
    r1->newels.Append (Element2d (p1, p6, p3, p4));
    r1->newels.Append (Element2d (p1, p2, p5, p6));
    r1->deledges.Append ( { p2, p3 } );
    r1->onp = 6;
    r1->bonus = 0;
    rules.Append (r1);

    // three quads to 2
    r1 = new ImprovementRule;
    r1->oldels.Append (Element2d (p1, p2, p3, p4));
    r1->oldels.Append (Element2d (p2, p5, p6, p3));
    r1->oldels.Append (Element2d (p3, p6, p7, p4));
    r1->newels.Append (Element2d (p1, p2, p5, p4));
    r1->newels.Append (Element2d (p4, p5, p6, p7));
    r1->deledges.Append ( { p2, p3 } );
    r1->deledges.Append ( { p3, p4 } );
    r1->deledges.Append ( { p3, p6 } );
    r1->onp = 7;
    r1->bonus = -1;
    rules.Append (r1);

    // quad + 2 connected trigs to quad
    r1 = new ImprovementRule;
    r1->oldels.Append (Element2d (p1, p2, p3, p4));
    r1->oldels.Append (Element2d (p2, p5, p3));
    r1->oldels.Append (Element2d (p3, p5, p4));
    r1->newels.Append (Element2d (p1, p2, p5, p4));
    r1->deledges.Append ( { p2, p3 } ); 
    r1->deledges.Append ( { p3, p4 } );
    r1->deledges.Append ( { p3, p5 } );
    r1->onp = 5;
    r1->bonus = 0;
    rules.Append (r1);

    // quad + 2 non-connected trigs to quad (a and b)
    r1 = new ImprovementRule;
    r1->oldels.Append (Element2d (p1, p2, p3, p4));
    r1->oldels.Append (Element2d (p2, p6, p3));
    r1->oldels.Append (Element2d (p1, p4, p5));
    r1->newels.Append (Element2d (p1, p3, p4, p5));
    r1->newels.Append (Element2d (p1, p2, p6, p3));
    r1->deledges.Append ( { p1, p4 } );
    r1->deledges.Append ( { p2, p3 } );
    r1->onp = 6;
    r1->bonus = 0;
    rules.Append (r1);

    r1 = new ImprovementRule;
    r1->oldels.Append (Element2d (p1, p2, p3, p4));
    r1->oldels.Append (Element2d (p2, p6, p3));
    r1->oldels.Append (Element2d (p1, p4, p5));
    r1->newels.Append (Element2d (p1, p2, p4, p5));
    r1->newels.Append (Element2d (p4, p2, p6, p3));
    r1->deledges.Append ( { p1, p4 } );
    r1->deledges.Append ( { p2, p3 } );
    r1->onp = 6;
    r1->bonus = 0;
    rules.Append (r1);

    // two quad + trig -> one quad + trig
    r1 = new ImprovementRule;
    r1->oldels.Append (Element2d (p1, p2, p3, p4));
    r1->oldels.Append (Element2d (p2, p5, p6, p3));
    r1->oldels.Append (Element2d (p4, p3, p6));
    r1->newels.Append (Element2d (p1, p2, p6, p4));
    r1->newels.Append (Element2d (p2, p5, p6));
    r1->deledges.Append ( { p2, p3 } );
    r1->deledges.Append ( { p3, p4 } );
    r1->deledges.Append ( { p3, p6 } );
    r1->onp = 6;
    r1->bonus = -1;
    rules.Append (r1);

    // swap quad + trig (a and b)
    r1 = new ImprovementRule;
    r1->oldels.Append (Element2d (p1, p2, p3, p4));
    r1->oldels.Append (Element2d (p2, p5, p3));
    r1->newels.Append (Element2d (p2, p5, p3, p4));
    r1->newels.Append (Element2d (p1, p2, p4));
    r1->deledges.Append ( { p2, p3 } );
    r1->onp = 5;
    r1->bonus = 0;
    rules.Append (r1);

    r1 = new ImprovementRule;
    r1->oldels.Append (Element2d (p1, p2, p3, p4));
    r1->oldels.Append (Element2d (p2, p5, p3));
    r1->newels.Append (Element2d (p1, p2, p5, p3));
    r1->newels.Append (Element2d (p1, p3, p4));
    r1->deledges.Append ( { p2, p3 } );
    r1->onp = 5;
    r1->bonus = 0;
    rules.Append (r1);


    // 2 quads to quad + 2 trigs
    r1 = new ImprovementRule;
    r1->oldels.Append (Element2d (p1, p2, p3, p4));
    r1->oldels.Append (Element2d (p3, p2, p5, p6));
    r1->newels.Append (Element2d (p1, p5, p6, p4));
    r1->newels.Append (Element2d (p1, p2, p5));
    r1->newels.Append (Element2d (p4, p6, p3));
    r1->deledges.Append ( { p2, p3 } );
    r1->onp = 6;
    r1->bonus = 0;
    //    rules.Append (r1);




    NgArray<int> mapped(rules.Size());
    NgArray<int> used(rules.Size());
    used = 0;
    mapped = 0;

  
    for (int ri = 0; ri < rules.Size(); ri++)
      {
	ImprovementRule & rule = *rules[ri];
	rule.incelsonnode.SetSize (rule.onp);
	rule.reused.SetSize (rule.onp);

        /*
	for (int j = 0; j < rule.onp; j++)
	  {
	    rule.incelsonnode[j] = 0;
	    rule.reused[j] = 0;
	  }
        */
        rule.incelsonnode = 0;
        rule.reused = 0;

        /*
	for (int j = 1; j <= rule.oldels.Size(); j++)
	  {
	    const Element2d & el = rule.oldels.Elem(j);
	    for (int k = 1; k <= el.GetNP(); k++)
	      rule.incelsonnode.Elem(el.PNum(k))--;
	  }
        */
        for (const Element2d & el : rule.oldels)
          for (PointIndex pi : el.PNums())
            rule.incelsonnode[pi]--;

	for (int j = 1; j <= rule.newels.Size(); j++)
	  {
	    const Element2d & el = rule.newels.Elem(j);
	    for (int k = 1; k <= el.GetNP(); k++)
	      {
		rule.incelsonnode[el.PNum(k)]++;
		rule.reused[el.PNum(k)] = 1;
	      }
	  }
      }



  
    DynamicTable<int,PointIndex> elonnode(np);
    Array<int,PointIndex> nelonnode(np);
    TABLE<SurfaceElementIndex> nbels(ne);

    nelonnode = -4;

    for (SurfaceElementIndex sei = 0; sei < ne; sei++)
      {
	const Element2d & el = mesh[sei];

	if (el.GetIndex() == faceindex && !el.IsDeleted())
	  {
	    for (int j = 0; j < el.GetNP(); j++)
	      elonnode.Add (el[j], sei);
	  }
	if(!el.IsDeleted())
	  {
	    for (int j = 0; j < el.GetNP(); j++)
	      nelonnode[el[j]]++;
	  }
      }

    for (SurfaceElementIndex sei = 0; sei < ne; sei++)
      {
	const Element2d & el = mesh[sei];
	if (el.GetIndex() == faceindex && !el.IsDeleted())
	  {
	    for (int j = 0; j < el.GetNP(); j++)
	      {
		for (int k = 0; k < elonnode[el[j]].Size(); k++)
		  {
		    int nbel = elonnode[el[j]] [k];
		    bool inuse = false;
		    for (int l = 0; l < nbels[sei].Size(); l++)
		      if (nbels[sei][l] == nbel)
			inuse = true;
		    if (!inuse)
		      nbels.Add (sei, nbel);
		  }
	      }
	  }
      }


    for (int ri = 0; ri < rules.Size(); ri++)
      {
	const ImprovementRule & rule = *rules[ri];
      
	elmap.SetSize (rule.oldels.Size());
	elrot.SetSize (rule.oldels.Size());
	pmap.SetSize (rule.onp);
	pgi.SetSize (rule.onp);


	for (SurfaceElementIndex sei = 0; sei < ne; sei++)
	  {
	    if (multithread.terminate)
	      break;
	    if (mesh[sei].IsDeleted()) continue;

	    elmap[0] = sei;
	    NgFlatArray<SurfaceElementIndex> neighbours = nbels[sei];
	    
	    for (elrot[0] = 0; elrot[0] < mesh[sei].GetNP(); elrot[0]++)
	      {
		const Element2d & el0 = mesh[sei];
		const Element2d & rel0 = rule.oldels[0];

		if (el0.GetIndex() != faceindex) continue;
		if (el0.IsDeleted()) continue;
		if (el0.GetNP() != rel0.GetNP()) continue;

		// pmap = PointIndex (-1);
                for (auto & p : pmap) p.Invalidate();
 
		for (int k = 0; k < el0.GetNP(); k++)
		  {
		    pmap[rel0[k]] = el0.PNumMod(k+elrot[0]+1);
		    pgi[rel0[k]] = el0.GeomInfoPiMod(k+elrot[0]+1);
		  }
		
		ok = 1;
		for (int i = 1; i < elmap.Size(); i++)
		  {
		    // try to find a mapping for reference-element i

		    const Element2d & rel = rule.oldels[i];
		    bool possible = 0;

		    for (elmap[i] = 0; elmap[i] < neighbours.Size(); elmap[i]++)
		      {
			const Element2d & el = mesh[neighbours[elmap[i]]];
			if (el.IsDeleted()) continue;
			if (el.GetNP() != rel.GetNP()) continue;

			for (elrot[i] = 0; elrot[i] < rel.GetNP(); elrot[i]++)
			  {
			    possible = 1;

			    for (int k = 0; k < rel.GetNP(); k++)
			      if (pmap[rel[k]].IsValid() && 
				  pmap[rel[k]] != el.PNumMod (k+elrot[i]+1))
				possible = 0;

			    if (possible) 
			      {
				for (int k = 0; k < el.GetNP(); k++)
				  {
				    pmap[rel[k]] = el.PNumMod(k+elrot[i]+1);
				    pgi[rel[k]] = el.GeomInfoPiMod(k+elrot[i]+1);
				  }
				break;
			      }
			  }
			if (possible) break;
		      }

		    if (!possible) 
		      {
			ok = 0;
			break;
		      }

		    elmap[i] = neighbours[elmap[i]];
		  }

		for(int i=0; ok && i<rule.deledges.Size(); i++)
		  {
		    ok = !mesh.IsSegment(pmap[rule.deledges[i][0]],
					 pmap[rule.deledges[i][1]]);
		  }
								    
								    
		
		
		if (!ok) continue;

		mapped[ri]++;

		olddef = 0;
		for (auto j : pmap.Range())
		  olddef += sqr (nelonnode[pmap[j]]);
		olddef += rule.bonus;

		newdef = 0;
		for (auto j : pmap.Range())
		  if (rule.reused[j])
		    newdef += sqr (nelonnode[pmap[j]] + 
				   rule.incelsonnode[j]); 

		if (newdef > olddef)
		  continue;

		// calc metric badness
		double bad1 = 0, bad2 = 0;
		// SelectSurfaceOfPoint (mesh.Point(pmap.Get(1)), pgi.Get(1));
		// auto n = geo.GetNormal(surfnr, mesh.Point(pmap.Get(1)), &pgi.Elem(1));
                auto n = geo.GetNormal(surfnr, mesh.Point(pmap.First()), pgi.Data());
		  
		for (int j = 0; j < rule.oldels.Size(); j++)
		  bad1 += mesh[elmap[j]].CalcJacobianBadness (mesh.Points(), n);
		  
		// check new element:
		for (int j = 1; j <= rule.newels.Size(); j++)
		  {
		    const Element2d & rnel = rule.newels.Get(j);
		    Element2d nel(rnel.GetNP());
		    for (int k = 1; k <= rnel.GetNP(); k++)
		      nel.PNum(k) = pmap[rnel.PNum(k)];

		    bad2 += nel.CalcJacobianBadness (mesh.Points(), n);
		  }

		if (bad2 > 1e3) continue;

		if (newdef == olddef && bad2 > bad1) continue;
		  

		// generate new element:
		for (int j = 1; j <= rule.newels.Size(); j++)
		  {
		    const Element2d & rnel = rule.newels.Get(j);
		    Element2d nel(rnel.GetNP());
		    nel.SetIndex (faceindex);
		    for (int k = 1; k <= rnel.GetNP(); k++)
		      {
			nel.PNum(k) = pmap[rnel.PNum(k)];
			nel.GeomInfoPi(k) = pgi[rnel.PNum(k)];
		      }
		      
		    mesh.AddSurfaceElement(nel);
		  }
		  
		for (int j = 0; j < rule.oldels.Size(); j++)
		  mesh.Delete (elmap[j]);

                for (PointIndex j : pmap.Range())
		  nelonnode[pmap[j]] += rule.incelsonnode[j];

		used[ri]++;
	      }
	  }
      }

    mesh.Compress();

    for (int ri = 0; ri < rules.Size(); ri++)
      {
	PrintMessage (5, "rule ", ri+1, " ",
		      mapped[ri], "/", used[ri], " mapped/used");
      }
  }




}
