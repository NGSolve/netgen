#include <mystdlib.h>

#include "meshing3.hpp"
#include "findip.hpp"
#include "findip2.hpp"

namespace netgen
{

double minother;
double minwithoutother;



MeshingStat3d :: MeshingStat3d ()
{
  cntsucc = cnttrials = cntelem = qualclass = 0;
  vol0 = h = 1;
  problemindex = 1;
}  
  

Meshing3 :: Meshing3 (const string & rulefilename) 
{
  tolfak = 1;

  LoadRules (rulefilename.c_str(), NULL);
  adfront = make_unique<AdFront3>();

  problems.SetSize (rules.Size());
  foundmap.SetSize (rules.Size());
  canuse.SetSize (rules.Size());
  ruleused.SetSize (rules.Size());

  /*
  for (int i = 1; i <= rules.Size(); i++)
    {
      problems.Elem(i) = new char[255];
      foundmap.Elem(i) = 0;
      canuse.Elem(i) = 0;
      ruleused.Elem(i) = 0;
    }
  */

  foundmap = 0;
  canuse = 0;
  ruleused = 0;
}


Meshing3 :: Meshing3 (const char ** rulep)
{
  tolfak = 1;

  LoadRules (NULL, rulep);
  adfront = make_unique<AdFront3>();

  problems.SetSize (rules.Size());
  foundmap.SetSize (rules.Size());
  canuse.SetSize (rules.Size());
  ruleused.SetSize (rules.Size());

  /*
  for (int i = 0; i < rules.Size(); i++)
    {
      problems[i] = new char[255];
      foundmap[i] = 0;
      canuse[i]   = 0;
      ruleused[i] = 0;
    }
  */
  foundmap = 0;
  canuse = 0;
  ruleused = 0;
}

Meshing3 :: ~Meshing3 ()
{
  // delete adfront;
  /*
  for (int i = 0; i < rules.Size(); i++)
    {
      delete [] problems[i];
      delete rules[i];
    }
  */
}



/*
  // was war das ????
static double CalcLocH (const Array<Point<3>> & locpoints,
                        const Array<MiniElement2d> & locfaces,
                        double h)
{
  return h;

  
  int i, j;
  double hi, h1, d, dn, sum, weight, wi;
  Point<3> p0, pc;
  Vec<3> n, v1, v2;

  p0.X() = p0.Y() = p0.Z() = 0;
  for (j = 1; j <= 3; j++)
    {
      p0.X() += locpoints.Get(locfaces.Get(1).PNum(j)).X();
      p0.Y() += locpoints.Get(locfaces.Get(1).PNum(j)).Y();
      p0.Z() += locpoints.Get(locfaces.Get(1).PNum(j)).Z();
    }
  p0.X() /= 3; p0.Y() /= 3; p0.Z() /= 3;
  
  v1 = locpoints.Get(locfaces.Get(1)[1]) -
    locpoints.Get(locfaces.Get(1)[0]);
  v2 = locpoints.Get(locfaces.Get(1)[2]) -
    locpoints.Get(locfaces.Get(1)[0]);

  h1 = v1.Length();
  n = Cross (v2, v1);
  n /= n.Length();

  sum = 0;
  weight = 0;

  for(int i = 1; i <= locfaces.Size(); i++)
    {
      pc.X() = pc.Y() = pc.Z() = 0;
      for (j = 1; j <= 3; j++)
        {
          pc.X() += locpoints.Get(locfaces.Get(i).PNum(j)).X();
          pc.Y() += locpoints.Get(locfaces.Get(i).PNum(j)).Y();
          pc.Z() += locpoints.Get(locfaces.Get(i).PNum(j)).Z();
        }
      pc.X() /= 3; pc.Y() /= 3; pc.Z() /= 3;

      d = Dist (p0, pc);
      dn = n * (pc - p0);
      hi = Dist (locpoints.Get(locfaces.Get(i)[0]),
                 locpoints.Get(locfaces.Get(i)[1]));
                 
      if (dn > -0.2 * h1)
        {
          wi = 1 / (h1 + d);
          wi *= wi;
        }
      else
        wi = 0;

      sum += hi * wi;
      weight += wi;
    }

  return sum/weight;
}
*/

Front3PointIndex Meshing3 :: AddPoint (const Point<3> & p, PointIndex globind)
{
  Front3PointIndex fpi = adfront -> AddPoint (p, globind);
  if (globind >= glob2front.Range().Next())
    {
      PointIndex oldnext = glob2front.Range().Next();
      glob2front.SetSize (globind.Nr0()+1);
      glob2front.Range(oldnext, END) = Front3PointIndex::INVALID;
    }
  glob2front[globind] = fpi;
  return fpi;
}  

void Meshing3 :: AddBoundaryElement (const Element2d & elem)
{
  FrontElement2d mini(elem.GetNP());
  for (int j = 0; j < elem.GetNP(); j++)
    mini[j] = glob2front[elem[j]];
  adfront -> AddFace(mini);
}  

int Meshing3 :: AddConnectedPair (PointIndices<2> apair)
{
  return adfront -> AddConnectedPair ( { glob2front[apair[0]], glob2front[apair[1]] } );
}

MESHING3_RESULT Meshing3 :: 
GenerateMesh (Mesh & mesh, const MeshingParameters & mp)
{
  static Timer t("Meshing3::GenerateMesh"); RegionTimer reg(t);
  // static Timer meshing3_timer_a("Meshing3::GenerateMesh a", 2);
  // static Timer meshing3_timer_b("Meshing3::GenerateMesh b", 2);
  // static Timer meshing3_timer_c("Meshing3::GenerateMesh c", 1);
  // static Timer meshing3_timer_d("Meshing3::GenerateMesh d", 2);
  // static int meshing3_timer = NgProfiler::CreateTimer ("Meshing3::GenerateMesh");
  // static int meshing3_timer_a = NgProfiler::CreateTimer ("Meshing3::GenerateMesh a");
  // static int meshing3_timer_b = NgProfiler::CreateTimer ("Meshing3::GenerateMesh b");
  // static int meshing3_timer_c = NgProfiler::CreateTimer ("Meshing3::GenerateMesh c");
  // static int meshing3_timer_d = NgProfiler::CreateTimer ("Meshing3::GenerateMesh d");
  // RegionTimer reg (meshing3_timer);


  Array<Point<3>, LocalPointIndex> locpoints;      // local points
  Array<MiniElement2d> locfaces;                   // local faces
  Array<Front3PointIndex, LocalPointIndex> pindex;  // mapping from local to front point numbering
  Array<int, LocalPointIndex> allowpoint;         // point is allowed (0/1/2) ?
  Array<int> findex;                             // mapping from local to front face numbering
  //INDEX_2_HASHTABLE<int> connectedpairs(100);    // connecgted pairs for prism meshing

  Array<Point<3>, LocalPointIndex> plainpoints;    // points in reference coordinates
  // Array<int> delpoints;   // points to be deleted
  Array<int> delfaces;    // lines to be deleted
  Array<LocalElement> locelements;       // new generated elements

  int oldnp, oldnf;
  int found;
  referencetransform trans;
  Point<3> inp = Point<3>(0,0,0);
  float err;

  int locfacesplit;             //index for faces in outer area
  
  bool loktestmode = false;

  int uselocalh = mp.uselocalh;

  // int giveuptol = mp.giveuptol; // 
  MeshingStat3d stat;      // statistics
  int plotstat_oldne = -1;

  
  // for star-shaped domain meshing
  Array<MeshPoint, LocalPointIndex> grouppoints;      
  Array<MiniElement2d> groupfaces;
  Array<Front3PointIndex, LocalPointIndex> grouppindex;
  Array<int> groupfindex;
  
  
  float minerr;
  int hasfound;
  [[maybe_unused]] double tetvol;
  // int giveup = 0;

  
  Array<Point<3>> tempnewpoints;
  Array<MiniElement2d> tempnewfaces;
  Array<int> tempdelfaces;
  Array<LocalElement> templocelements;


  stat.h = mp.maxh;

  adfront->SetStartFront (mp.baseelnp);


  found = 0;
  stat.vol0 = adfront -> Volume();
  tetvol = 0;

  stat.qualclass = 1;

  while (1)
    {
      if (multithread.terminate)
        throw NgException ("Meshing stopped");

      // break if advancing front is empty
      if (!mp.baseelnp && adfront->Empty())
        break;

      // break, if advancing front has no elements with
      // mp.baseelnp nodes  
      if (mp.baseelnp && adfront->Empty (mp.baseelnp))
        break;

      locpoints.SetSize(0);
      locfaces.SetSize(0);
      locelements.SetSize(0);
      pindex.SetSize(0);
      findex.SetSize(0);

      ClosedHashTable<IVec<2>,int> connectedpairs(100);  // connected pairs for prism meshing
      
      // select base-element (will be locface[1])
      // and get local environment of radius (safety * h)


      int baseelem = adfront -> SelectBaseElement ();
      if (mp.baseelnp && adfront->GetFace (baseelem).GetNP() != mp.baseelnp)
        {
          adfront->IncrementClass (baseelem);     
          continue;
        }

      const FrontElement2d & bel = adfront->GetFace (baseelem);
      const Point<3> p1 = adfront->GetPoint (bel[0]);
      const Point<3> p2 = adfront->GetPoint (bel[1]);
      const Point<3> p3 = adfront->GetPoint (bel[2]);
      

      Point<3> pmid = Center (p1, p2, p3);

      double his = (Dist (p1, p2) + Dist(p1, p3) + Dist(p2, p3)) / 3;
      double hshould = mesh.GetH (pmid);
      
      if (adfront->GetFace (baseelem).GetNP() == 4)
        hshould = max2 (his, hshould);

      double hmax = (his > hshould) ? his : hshould;
      
      // qualclass should come from baseelem !!!!!
      double hinner = hmax * (1 + stat.qualclass);
      double houter = hmax * (1 + 2 * stat.qualclass);

      // meshing3_timer_a.Start();
      stat.qualclass =
        adfront -> GetLocals (baseelem, locpoints, locfaces, 
                              pindex, findex, connectedpairs,
                              houter, hinner,
                              locfacesplit);
      // meshing3_timer_a.Stop();

      // (*testout) << "locfaces = " << endl << locfaces << endl;


      //loktestmode = 1;
      testmode = loktestmode;  //changed 
      // loktestmode = testmode =  (adfront->GetFace (baseelem).GetNP() == 4) && (rules.Size() == 5);

      loktestmode = stat.qualclass > 5;
      

      if (loktestmode)
        {
          (*testout) << "baseel = " << baseelem << ", ind = " << findex[0] << endl;
          Front3PointIndex pi1 = pindex[locfaces[0][0]];
          Front3PointIndex pi2 = pindex[locfaces[0][1]];
          Front3PointIndex pi3 = pindex[locfaces[0][2]];
          (*testout) << "pi = " << pi1 << ", " << pi2 << ", " << pi3 << endl;
        }


      if (testmode)
        {
          (*testout) << "baseelem = " << baseelem << " qualclass = " << stat.qualclass << endl;
          (*testout) << "locpoints = " << endl << locpoints << endl;
          (*testout) << "connected = " << endl << connectedpairs << endl;
        }



      // loch = CalcLocH (locpoints, locfaces, h);
      
      stat.nff = adfront->GetNF();
      stat.vol = adfront->Volume();
      if (stat.vol < 0) break;

      oldnp = locpoints.Size();
      oldnf = locfaces.Size();


      allowpoint.SetSize(locpoints.Size());
      if (uselocalh && stat.qualclass <= 3)
        //for(int i = 1; i <= allowpoint.Size(); i++)
        for(auto i : allowpoint.Range())
          {
            allowpoint[i] =
              (mesh.GetH (locpoints[i]) > 0.4 * hshould / mp.sloppy) ? 2 : 1;
          }
      else
        allowpoint = 2;


      
      if (stat.qualclass >= mp.starshapeclass &&
          mp.baseelnp != 4)   
        {
          // RegionTimer reg1 (meshing3_timer_b);
          // star-shaped domain removing

          grouppoints.SetSize (0);
          groupfaces.SetSize (0);
          grouppindex.SetSize (0);
          groupfindex.SetSize (0);
          
          adfront -> GetGroup (findex[0], grouppoints, groupfaces, 
                               grouppindex, groupfindex);

          bool onlytri = 1;
          for (auto & f : groupfaces)
            if (f.GetNP() != 3) 
              onlytri = 0;
            
              
          if (onlytri && groupfaces.Size() <= 20 + 2*stat.qualclass &&
              FindInnerPoint (grouppoints, groupfaces, inp) &&
              !adfront->PointInsideGroup(grouppindex, groupfaces))
            {
              (*testout) << "inner point found" << endl;

              for(int i = 0; i < groupfaces.Size(); i++)
                adfront -> DeleteFace (groupfindex[i]);
              
              for (int i = 0; i < groupfaces.Size(); i++)
                for (int j = 1; j <= locfaces.Size(); j++)
                  if (findex[j-1] == groupfindex[i])
                    delfaces.Append (j);
              
              
              delfaces.SetSize (0);
              
              Element newel(TET);
              
              PointIndex npi = mesh.AddPoint (inp);
              newel.SetNP(4);
              newel[3] = npi;
              
              for(int i = 0; i < groupfaces.Size(); i++)
                {
                  for (int j = 0; j < 3; j++)
                    {
                      newel[j] = 
                        adfront->GetGlobalIndex 
                        (grouppindex[groupfaces[i][j]]);
                    }
                  mesh.AddVolumeElement (newel);
                }
              continue;
            }
        }
      
      found = 0;
      hasfound = 0;
      minerr = 1e6;

      //      int optother = 0;

      /*
      for(int i = 1; i <= locfaces.Size(); i++)
        {
          (*testout) << "Face " << i << ": ";
          for (j = 1; j <= locfaces.Get(i).GetNP(); j++)
            (*testout) << pindex.Get(locfaces.Get(i).PNum(j)) << " ";
          (*testout) << endl;
        }
      for(int i = 1; i <= locpoints.Size(); i++)
        {
          (*testout) << "p" << i 
                     << ", gi = " << pindex.Get(i) 
                     << " = " << locpoints.Get(i) << endl;
        }
        */

      minother = 1e10;
      minwithoutother = 1e10;

      bool impossible = 1;

      for (int rotind = 1; rotind <= locfaces[0].GetNP(); rotind++)
        {
          // set transformatino to reference coordinates

          if (locfaces[0].GetNP() == 3)
            {
              trans.Set (locpoints[locfaces[0].PNumMod(1+rotind)],
                         locpoints[locfaces[0].PNumMod(2+rotind)],
                         locpoints[locfaces[0].PNumMod(3+rotind)], hshould);
            }
          else
            {
              trans.Set (locpoints[locfaces[0].PNumMod(1+rotind)],
                         locpoints[locfaces[0].PNumMod(2+rotind)],
                         locpoints[locfaces[0].PNumMod(4+rotind)], hshould);
            }

          // trans.ToPlain (locpoints, plainpoints);

          plainpoints.SetSize (locpoints.Size());
          for (auto i : locpoints.Range())
            trans.ToPlain (locpoints[i], plainpoints[i]);
          
          for (auto i : allowpoint.Range())
            if (plainpoints[i](2) > 0)
              allowpoint[i] = false;

          stat.cnttrials++;


          if (stat.cnttrials % 100 == 0)
            {
              (*testout) << "\n";
              for(int i : canuse.Range()) 
                {
                  (*testout) << foundmap[i] << "/" 
                             << canuse[i] << "/"
                             << ruleused[i] << " map/can/use rule " << rules[i]->Name() << "\n";
                }
              (*testout) << endl;
            }

          // NgProfiler::StartTimer (meshing3_timer_c);
          // meshing3_timer_c.Start();

          found = ApplyRules (plainpoints, allowpoint, 
                              locfaces, locfacesplit, connectedpairs,
                              locelements, delfaces, 
                              stat.qualclass, mp.sloppy, rotind, err);

          if (found >= 0) impossible = 0;
          if (found < 0) found = 0;

          // meshing3_timer_c.Stop();
          // NgProfiler::StopTimer (meshing3_timer_c);    

          if (!found) loktestmode = 0;

          // RegionTimer reg2 (meshing3_timer_d);         
          
          if (loktestmode)
            {
              (*testout) << "plainpoints = " << endl << plainpoints << endl;
              (*testout) << "Applyrules found " << found << endl;
            }

          if (found) stat.cntsucc++;

          locpoints.SetSize (plainpoints.Size());
          /*
          for (int i = oldnp+1; i <= plainpoints.Size(); i++)
            trans.FromPlain (plainpoints.Elem(i), locpoints.Elem(i));
          */
          for (auto i : plainpoints.Range().Modify(oldnp,0))
            trans.FromPlain (plainpoints[i], locpoints[i]);
          


          // avoid meshing from large to small mesh-size
          if (uselocalh && found && stat.qualclass <= 3)
            {
              for (int i = 0; i < locelements.Size(); i++)
                {
                  Point<3> pmin = locpoints[locelements[i][0]];
                  Point<3> pmax = pmin;
                  for (int j = 2; j <= 4; j++)
                    {
                      const Point<3> & hp = locpoints[locelements[i].PNum(j)];
                      SetToMin (pmin, hp);
                      SetToMax (pmax, hp);
                    }

                  if (mesh.GetMinH (pmin, pmax) < 0.4 * hshould / mp.sloppy)
                    found = 0;
                }
            }
          if (found)
            {
              for (int i = 1; i <= locelements.Size(); i++)
                for (int j = 0; j < 4; j++)
                  {
                    const Point<3> & hp = locpoints[locelements[i-1][j]];
                    if (Dist (hp, pmid) > hinner)
                      found = 0;
                  }
            }


          if (found)
            ruleused[found-1]++;
          

          // plotstat->Plot(stat);
          
          if (stat.cntelem != plotstat_oldne)
            {
              plotstat_oldne = stat.cntelem;

              PrintMessageCR (5, "El: ", stat.cntelem,
                              //            << " trials: " << stat.cnttrials
                              " faces: ", stat.nff,
                              " vol = ", float(100 * stat.vol / stat.vol0));
  
              multithread.percent = 100 -  100.0 * stat.vol / stat.vol0;
            }


          if (found && (!hasfound || err < minerr) )
            {
              
              if (testmode)
                {
                  (*testout) << "found is active, 3" << endl;
                  //for(int i = 1; i <= plainpoints.Size(); i++)
                  for(auto i : plainpoints.Range())
                    {
                      (*testout) << "p";
                      if (i < pindex.Range().Next())
                        (*testout) << pindex[i] << ": ";
                      else
                        (*testout) << "new: ";
                      (*testout) << plainpoints[i] << endl;
                    }
                }
              
              
              
              hasfound = found;
              minerr = err;
              
              tempnewpoints.SetSize (0);
              // for(int i = oldnp+1; i <= locpoints.Size(); i++)
              for (auto i : locpoints.Range().Modify(oldnp,0))
                tempnewpoints.Append (locpoints[i]);
              
              tempnewfaces.SetSize (0);
              // for(int i = oldnf+1; i <= locfaces.Size(); i++)
              for (auto i : locfaces.Range().Modify(oldnf,0))              
                tempnewfaces.Append (locfaces[i]);
              
              tempdelfaces.SetSize (0);
              // for(int i = 1; i <= delfaces.Size(); i++)
              for (auto i : delfaces.Range())
                tempdelfaces.Append (delfaces[i]);
              
              templocelements.SetSize (0);
              // for(int i = 1; i <= locelements.Size(); i++)
              for (auto i : locelements.Range())
                templocelements.Append (locelements[i]);

              /*
              optother =
                strcmp (problems[found], "other") == 0;
              */
            }
          
          locpoints.SetSize (oldnp);
          locfaces.SetSize (oldnf);
          delfaces.SetSize (0);
          locelements.SetSize (0);
        }
      
      

      if (hasfound)
        {

          /*
          if (optother)
            (*testout) << "Other is optimal" << endl;

          if (minother < minwithoutother)
            {
              (*testout) << "Other is better, " << minother << " less " << minwithoutother << endl;
            }
            */

          for (int i = 0; i < tempnewpoints.Size(); i++)
            locpoints.Append (tempnewpoints[i]);
          for (int i = 0; i < tempnewfaces.Size(); i++)
            locfaces.Append (tempnewfaces[i]);
          for (int i = 0; i < tempdelfaces.Size(); i++)
            delfaces.Append (tempdelfaces[i]);
          for (int i = 0; i < templocelements.Size(); i++)
            locelements.Append (templocelements[i]);


          if (loktestmode)
            {
              (*testout) << "apply rule" << endl;
              for (auto i : locpoints.Range())
                {
                  (*testout) << "p";
                  if (pindex.Range().Contains(i))
                    (*testout) << pindex[i] << ": ";
                  else
                    (*testout) << "new: ";
                  (*testout) << locpoints[i] << endl;
                }
            }



          pindex.SetSize(locpoints.Size());

          // for (int i = oldnp+1; i <= locpoints.Size(); i++)
          for (auto i : locpoints.Range().Modify(oldnp,0))
            {
              PointIndex globind = mesh.AddPoint (locpoints[i]);
              pindex[i] = adfront -> AddPoint (locpoints[i], globind);
            }

          for (int i = 0; i < locelements.Size(); i++)
            {
              Point<3> * hp1, * hp2, * hp3, * hp4;
              hp1 = &locpoints[locelements[i][0]];
              hp2 = &locpoints[locelements[i][1]];
              hp3 = &locpoints[locelements[i][2]];
              hp4 = &locpoints[locelements[i][3]];
              
              tetvol += (1.0 / 6.0) * ( Cross ( *hp2 - *hp1, *hp3 - *hp1) * (*hp4 - *hp1) );

              const LocalElement & locel = locelements[i];
              Element el(locel.GetNP());
              el.SetType (locel.GetType());
              for (int j = 1; j <= locel.GetNP(); j++)
                el.PNum(j) = adfront -> GetGlobalIndex (pindex[locel.PNum(j)]);

              mesh.AddVolumeElement (el);
              stat.cntelem++;
            }

          for(int i = oldnf; i < locfaces.Size(); i++)
            {
              FrontElement2d frontface(locfaces[i].GetNP());
              for (int j = 0; j < locfaces[i].GetNP(); j++)
                frontface[j] = pindex[locfaces[i][j]];
              adfront->AddFace (frontface);
            }
          
          for (int i = 0; i < delfaces.Size(); i++)
            adfront->DeleteFace (findex[delfaces[i]-1]);
        }
      else
        {
          adfront->IncrementClass (findex[0]);
          if (impossible && mp.check_impossible)
            {
              (*testout) << "skip face since it is impossible" << endl;
              for (int j = 0; j < 100; j++)
                adfront->IncrementClass (findex[0]);
            }
        }

      locelements.SetSize (0);
      // delpoints.SetSize(0);
      delfaces.SetSize(0);

      if (stat.qualclass >= mp.giveuptol)
        break;
    }
  
  PrintMessage (5, "");  // line feed after statistics

  for(int i : ruleused.Range())
    (*testout) << setw(4) << ruleused[i]
               << " times used rule " << rules[i] -> Name() << endl;


  if (!mp.baseelnp && adfront->Empty())
    return MESHING3_OK;

  if (mp.baseelnp && adfront->Empty (mp.baseelnp))
    return MESHING3_OK;

  if (stat.vol < -1e-15)
    return MESHING3_NEGVOL;

  return MESHING3_NEGVOL;
}




enum blocktyp { BLOCKUNDEF, BLOCKINNER, BLOCKBOUND, BLOCKOUTER };

void Meshing3 :: BlockFill (Mesh & mesh, double gh)
{
  static Timer t("Mesing3::BlockFill"); RegionTimer reg(t);
  
  PrintMessage (3, "Block-filling called (obsolete) ");

  int i;
  int n1, n2, n3, n, min1, min2, min3, max1, max2, max3;
  int changed, filled;
  double xmin(0), xmax(0), ymin(0), ymax(0), zmin(0), zmax(0);
  double xminb, xmaxb, yminb, ymaxb, zminb, zmaxb;
  //double rad = 0.7 * gh;
  
  for(int i = 1; i <= adfront->GetNP(); i++)
    {
      const auto & p = adfront->GetPoint(Front3PointIndex::FromNr1(i));
      if (i == 1)
        {
          xmin = xmax = p(0);
          ymin = ymax = p(1);
          zmin = zmax = p(2);
        }
      else
        {
          if (p(0) < xmin) xmin = p(0);
          if (p(0) > xmax) xmax = p(0);
          if (p(1) < ymin) ymin = p(1);
          if (p(1) > ymax) ymax = p(1);
          if (p(2) < zmin) zmin = p(2);
          if (p(2) > zmax) zmax = p(2);
        }
    }
  
  xmin -= 5 * gh;
  ymin -= 5 * gh;
  zmin -= 5 * gh;
  
  n1 = int ((xmax-xmin) / gh + 5);
  n2 = int ((ymax-ymin) / gh + 5);
  n3 = int ((zmax-zmin) / gh + 5);
  n = n1 * n2 * n3;
  
  PrintMessage (5, "n1 = ", n1, " n2 = ", n2, " n3 = ", n3);

  Array<blocktyp> inner(n);
  Array<PointIndex> pointnr(n);
  Array<Front3PointIndex> frontpointnr(n);


  // initialize inner to 1

  for (int i = 0; i < n; i++)
    inner[i] = BLOCKUNDEF;


  // set blocks cutting surfaces to 0

  for(int i = 1; i <= adfront->GetNF(); i++)
    {
      const FrontElement2d & el = adfront->GetFace(i);
      xminb = xmax; xmaxb = xmin;
      yminb = ymax; ymaxb = ymin;
      zminb = zmax; zmaxb = zmin;

      for (int j = 0; j < 3; j++)
        {
          const auto & p = adfront->GetPoint (el[j]);
          if (p(0) < xminb) xminb = p(0);
          if (p(0) > xmaxb) xmaxb = p(0);
          if (p(1) < yminb) yminb = p(1);
          if (p(1) > ymaxb) ymaxb = p(1);
          if (p(2) < zminb) zminb = p(2);
          if (p(2) > zmaxb) zmaxb = p(2);
        }

        

      double filldist = 0.2; // globflags.GetNumFlag ("filldist", 0.4);
      xminb -= filldist * gh;
      xmaxb += filldist * gh;
      yminb -= filldist * gh;
      ymaxb += filldist * gh;
      zminb -= filldist * gh;
      zmaxb += filldist * gh;

      min1 = int ((xminb - xmin) / gh) + 1;
      max1 = int ((xmaxb - xmin) / gh) + 1;
      min2 = int ((yminb - ymin) / gh) + 1;
      max2 = int ((ymaxb - ymin) / gh) + 1;
      min3 = int ((zminb - zmin) / gh) + 1;
      max3 = int ((zmaxb - zmin) / gh) + 1;


      for (int i1 = min1; i1 <= max1; i1++)
        for (int i2 = min2; i2 <= max2; i2++)
          for (int i3 = min3; i3 <= max3; i3++)
            inner[i3 + (i2-1) * n3 + (i1-1) * n2 * n3-1] = BLOCKBOUND;      
    }

  


  while (1)
    {
      int undefi = 0;
      Point<3> undefp;

      for (int i1 = 1; i1 <= n1 && !undefi; i1++)
        for (int i2 = 1; i2 <= n2 && !undefi; i2++)
          for (int i3 = 1; i3 <= n3 && !undefi; i3++)
            {
              i = i3 + (i2-1) * n3 + (i1-1) * n2 * n3;
              if (inner[i-1] == BLOCKUNDEF)
                {
                  undefi = i;
                  undefp(0) = xmin + (i1-0.5) * gh;
                  undefp(1) = ymin + (i2-0.5) * gh;
                  undefp(2) = zmin + (i3-0.5) * gh;
                }
            }
              
      if (!undefi)
        break;

      //      PrintMessage (5, "Test point: ", undefp);
      
      if (adfront -> Inside (undefp))
        {
          //      (*mycout) << "inner" << endl;
          inner[undefi-1] = BLOCKINNER;
        }
      else
        {
          //      (*mycout) << "outer" << endl;
          inner[undefi-1] = BLOCKOUTER;
        }

      do
        {
          changed = 0;
          for (int i1 = 1; i1 <= n1; i1++)
            for (int i2 = 1; i2 <= n2; i2++)
              for (int i3 = 1; i3 <= n3; i3++)
                {
                  i = i3 + (i2-1) * n3 + (i1-1) * n2 * n3;

                  for (int k = 1; k <= 3; k++)
                    {
                      int j = 0;
                      switch (k)
                        {
                        case 1: j = i + n2 * n3; break;
                        case 2: j = i + n3; break;
                        case 3: j = i + 1; break;
                        }
                  
                      if (j > n1 * n2 * n3) continue;

                      if (inner[i-1] == BLOCKOUTER && inner[j-1] == BLOCKUNDEF)
                        {
                          changed = 1;
                          inner[j-1] = BLOCKOUTER;
                        }
                      if (inner[j-1] == BLOCKOUTER && inner[i-1] == BLOCKUNDEF)
                        {
                          changed = 1;
                          inner[i-1] = BLOCKOUTER;
                        }
                      if (inner[i-1] == BLOCKINNER && inner[j-1] == BLOCKUNDEF)
                        {
                          changed = 1;
                          inner[j-1] = BLOCKINNER;
                        }
                      if (inner[j-1] == BLOCKINNER && inner[i-1] == BLOCKUNDEF)
                        {
                          changed = 1;
                          inner[i-1] = BLOCKINNER;
                        }
                    }
                }
        }
      while (changed); 

    }



  filled = 0;
  for(int i = 1; i <= n; i++)
    if (inner[i-1] == BLOCKINNER)
      {
        filled++;
      }
  PrintMessage (5, "Filled blocks: ", filled);

  for (int i = 0; i < n; i++)
    {
      pointnr[i] = PointIndex::INVALID;
      frontpointnr[i] = Front3PointIndex::INVALID;
    }
  
  for (int i1 = 1; i1 <= n1-1; i1++)
    for (int i2 = 1; i2 <= n2-1; i2++)
      for (int i3 = 1; i3 <= n3-1; i3++)
        {
          i = i3 + (i2-1) * n3 + (i1-1) * n2 * n3;
          if (inner[i-1] == BLOCKINNER)
            {
              for (int j1 = i1; j1 <= i1+1; j1++)
                for (int j2 = i2; j2 <= i2+1; j2++)
                  for (int j3 = i3; j3 <= i3+1; j3++)
                    {
                      int j = j3 + (j2-1) * n3 + (j1-1) * n2 * n3;
                      if (!pointnr[j-1].IsValid())
                        {
                          Point<3> hp(xmin + (j1-1) * gh, 
                                     ymin + (j2-1) * gh, 
                                     zmin + (j3-1) * gh);
                          pointnr[j-1] = mesh.AddPoint (hp);
                          frontpointnr[j-1] =
                            AddPoint (hp, pointnr[j-1]);

                        }
                    }
            }
        }


  for (int i1 = 2; i1 <= n1-1; i1++)
    for (int i2 = 2; i2 <= n2-1; i2++)
      for (int i3 = 2; i3 <= n3-1; i3++)
        {
          i = i3 + (i2-1) * n3 + (i1-1) * n2 * n3;
          if (inner[i-1] == BLOCKINNER)
            {
              PointIndex pn[9];
              pn[1] = pointnr[i-1];
              pn[2] = pointnr[i];
              pn[3] = pointnr[i+n3-1];
              pn[4] = pointnr[i+n3];
              pn[5] = pointnr[i+n2*n3-1];
              pn[6] = pointnr[i+n2*n3];
              pn[7] = pointnr[i+n2*n3+n3-1];
              pn[8] = pointnr[i+n2*n3+n3];
              static int elind[][4] =
              {
                { 1, 8, 2, 4 },
                { 1, 8, 4, 3 },
                { 1, 8, 3, 7 },
                { 1, 8, 7, 5 },
                { 1, 8, 5, 6 },
                { 1, 8, 6, 2 }
              };
              for (int j = 0; j < 6; j++)
                {
                  Element el(4);
                  for (int k = 1; k <= 4;  k++)
                    el.PNum(k) = pn[elind[j][k-1]];

                  mesh.AddVolumeElement (el);
                }
            }
        }



  for (int i1 = 2; i1 <= n1-1; i1++)
    for (int i2 = 2; i2 <= n2-1; i2++)
      for (int i3 = 2; i3 <= n3-1; i3++)
        {
          i = i3 + (i2-1) * n3 + (i1-1) * n2 * n3;
          if (inner[i-1] == BLOCKINNER)
            {    
              Front3PointIndex pi1, pi2, pi3, pi4;

              Front3PointIndex pn1 = frontpointnr[i-1];
              Front3PointIndex pn2 = frontpointnr[i];
              Front3PointIndex pn3 = frontpointnr[i+n3-1];
              Front3PointIndex pn4 = frontpointnr[i+n3];
              Front3PointIndex pn5 = frontpointnr[i+n2*n3-1];
              Front3PointIndex pn6 = frontpointnr[i+n2*n3];
              Front3PointIndex pn7 = frontpointnr[i+n2*n3+n3-1];
              Front3PointIndex pn8 = frontpointnr[i+n2*n3+n3];

              for (int k = 1; k <= 6; k++)
                {
                  int j = 0;
                  switch (k)
                    {
                    case 1: // j3 = i3+1
                      j = i + 1;
                      pi1 = pn2;
                      pi2 = pn6;
                      pi3 = pn4;
                      pi4 = pn8;
                      break;
                    case 2: // j3 = i3-1
                      j = i - 1;
                      pi1 = pn1;
                      pi2 = pn3;
                      pi3 = pn5;
                      pi4 = pn7;
                      break;
                    case 3: // j2 = i2+1
                      j = i + n3;
                      pi1 = pn3;
                      pi2 = pn4;
                      pi3 = pn7;
                      pi4 = pn8;
                      break;
                    case 4: // j2 = i2-1
                      j = i - n3;
                      pi1 = pn1;
                      pi2 = pn5;
                      pi3 = pn2;
                      pi4 = pn6;
                      break;
                    case 5: // j1 = i1+1
                      j = i + n3*n2;
                      pi1 = pn5;
                      pi2 = pn7;
                      pi3 = pn6;
                      pi4 = pn8;
                      break;
                    case 6: // j1 = i1-1
                      j = i - n3*n2;
                      pi1 = pn1;
                      pi2 = pn2;
                      pi3 = pn3;
                      pi4 = pn4;
                      break;
                    }

                  if (inner[j-1] == BLOCKBOUND)
                    {
                      FrontElement2d face;
                      face[0] = pi4;
                      face[1] = pi1;
                      face[2] = pi3;
                      adfront->AddFace (face);

                      face[0] = pi1;
                      face[1] = pi4;
                      face[2] = pi2;
                      adfront->AddFace (face);

                    }
                }
            }
        }
}


/*
static const AdFront3 * locadfront;
static int TestInner (const Point<3> & p)
{
  return locadfront->Inside (p);
}
static int TestSameSide (const Point<3> & p1, const Point<3> & p2)
{
  return locadfront->SameSide (p1, p2);
}
*/



void Meshing3 :: BlockFillLocalH (Mesh & mesh, 
                                  const MeshingParameters & mp)
{
  static Timer t("Mesing3::BlockFillLocalH"); RegionTimer reg(t);
  
  double filldist = mp.filldist;
  
  // (*testout) << "blockfill local h" << endl;
  // (*testout) << "rel filldist = " << filldist << endl;
  PrintMessage (3, "blockfill local h");


  Array<Point<3> > npoints;
  
  adfront -> CreateTrees();

  Box<3> bbox ( Box<3>::EMPTY_BOX );
  double maxh = 0;

  for (int i = 1; i <= adfront->GetNF(); i++)
    {
      const FrontElement2d & el = adfront->GetFace(i);
      for (int j = 1; j <= 3; j++)
        {
          const auto & p1 = adfront->GetPoint (el.PNumMod(j));
          const auto & p2 = adfront->GetPoint (el.PNumMod(j+1));

          double hi = Dist (p1, p2);
          if (hi > maxh) maxh = hi;

          bbox.Add (p1);
        }
    }


  Point<3> mpmin = bbox.PMin();
  Point<3> mpmax = bbox.PMax();
  Point<3> mpc = Center (mpmin, mpmax);
  double d = max3(mpmax(0)-mpmin(0), 
                  mpmax(1)-mpmin(1), 
                  mpmax(2)-mpmin(2)) / 2;
  mpmin = mpc - Vec<3> (d, d, d);
  mpmax = mpc + Vec<3> (d, d, d);
  Box3d meshbox (mpmin, mpmax);

  LocalH loch2 (mpmin, mpmax, 1);

  if (mp.maxh < maxh) maxh = mp.maxh;

  auto loch_ptr = mesh.LocalHFunction().Copy(bbox);
  auto & loch = *loch_ptr;

  bool changed;
  static Timer t1("loop1");
  t1.Start();
  do 
    {
      loch.ClearFlags();

      static Timer tbox("adfront-bbox");
      tbox.Start();
      for (int i = 1; i <= adfront->GetNF(); i++)
        {
          const FrontElement2d & el = adfront->GetFace(i);
          
          Box<3> bbox (adfront->GetPoint (el[0]));
          bbox.Add (adfront->GetPoint (el[1]));
          bbox.Add (adfront->GetPoint (el[2]));


          double filld = filldist * bbox.Diam();
          bbox.Increase (filld);
      
          loch.CutBoundary (bbox); // .PMin(), bbox.PMax());
        }
      tbox.Stop();

      //      locadfront = adfront;
      loch.FindInnerBoxes (*adfront, NULL);

      npoints.SetSize(0);
      loch.GetInnerPoints (npoints);

      changed = false;
      for (int i = 0; i < npoints.Size(); i++)
        {
          if (loch.GetH(npoints[i]) > 1.5 * maxh)
            {
              loch.SetH (npoints[i], maxh);
              changed = true;
            }
        }
    }
  while (changed);
  t1.Stop();

  if (debugparam.slowchecks)
    (*testout) << "Blockfill with points: " << endl;
  for (int i = 0; i < npoints.Size(); i++)
    {
      if (meshbox.IsIn (npoints[i]))
        {
          PointIndex gpnum = mesh.AddPoint (npoints[i]);
          adfront->AddPoint (npoints[i], gpnum);

          if (debugparam.slowchecks)
            {
              (*testout) << npoints[i] << endl;
              if (!adfront->Inside(npoints[i]))
                {
                  cout << "add outside point" << endl;
                  (*testout) << "outside" << endl;
                }
            }

        }
    }

  

  // find outer points
  
  static Timer tloch2("build loch2");
  tloch2.Start();
  loch2.ClearFlags();

  for (int i = 1; i <= adfront->GetNF(); i++)
    {
      const FrontElement2d & el = adfront->GetFace(i);
      Point<3> pmin = adfront->GetPoint (el[0]);
      Point<3> pmax = pmin;
      
      for (int j = 2; j <= 3; j++)
        {
          const auto & p = adfront->GetPoint (el.PNum(j));
          SetToMin (pmin, p);
          SetToMax (pmax, p);
        }
      
      loch2.SetH (Center (pmin, pmax), Dist (pmin, pmax));
    }

  for (int i = 1; i <= adfront->GetNF(); i++)
    {
      const FrontElement2d & el = adfront->GetFace(i);
      Point<3> pmin = adfront->GetPoint (el[0]);
      Point<3> pmax = pmin;
      
      for (int j = 2; j <= 3; j++)
        {
          const auto & p = adfront->GetPoint (el.PNum(j));
          SetToMin (pmin, p);
          SetToMax (pmax, p);
        }
      
      double filld = filldist * Dist (pmin, pmax);
      pmin = pmin - Vec<3> (filld, filld, filld);
      pmax = pmax + Vec<3> (filld, filld, filld);
      // loch2.CutBoundary (pmin, pmax);
      loch2.CutBoundary (Box<3> (pmin, pmax)); // pmin, pmax);
    }
  tloch2.Stop();

  // locadfront = adfront;
  loch2.FindInnerBoxes (*adfront, NULL);

  npoints.SetSize(0);
  loch2.GetOuterPoints (npoints);
  
  for (int i = 0; i < npoints.Size(); i++)
    {
      if (meshbox.IsIn (npoints[i]))
        {
          PointIndex gpnum = mesh.AddPoint (npoints[i]);
          adfront->AddPoint (npoints[i], gpnum);
        }
    }  
}

}
