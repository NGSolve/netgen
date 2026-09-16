#include <mystdlib.h>

#include <myadt.hpp>
#include <linalg.hpp>
#include <gprim.hpp>

#include <meshing.hpp>

#include "stlgeom.hpp"

namespace netgen
{

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++  EDGE DATA     ++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


/*
void STLEdgeData :: Write(ofstream& of) const
{
  of // << angle << " "
     << p1 << " "
     << p2 << " "
     << lt << " "
     << rt << " "
    //     << status
     << endl;
}

void STLEdgeData :: Read(ifstream& ifs)
{
  // ifs >> angle;
  ifs >> p1;
  ifs >> p2;
  ifs >> lt;
  ifs >> rt;
  //  ifs >> status;
}


int STLEdgeData :: GetStatus () const
{
  if (topedgenr <= 0 || topedgenr > top->GetNTE()) return 0;
  return top->GetTopEdge (topedgenr).GetStatus(); 
}

void STLEdgeData ::SetStatus (int stat)
{
  if (topedgenr >= 1 && topedgenr <= top->GetNTE())
    top->GetTopEdge (topedgenr).SetStatus(stat); 
}


float STLEdgeData :: CosAngle() const
{
  return top->GetTopEdge (topedgenr).CosAngle(); 
}



void STLEdgeDataList :: ResetAll()
{
  int i;
  for (i = 1; i <= edgedata.Size(); i++)
    {
      edgedata.Elem(i).SetUndefined();
    }
}

void STLEdgeDataList :: ResetCandidates()
{
  int i;
  for (i = 1; i <= edgedata.Size(); i++)
    {
      if (edgedata.Get(i).Candidate())
        {edgedata.Elem(i).SetUndefined();}
    }
}

int STLEdgeDataList :: GetNConfEdges() const
{
  int i;
  int cnt = 0;
  for (i = 1; i <= edgedata.Size(); i++)
    {
      if (edgedata.Get(i).Confirmed()) {cnt++;}
    }
  return cnt;
}

void STLEdgeDataList :: ConfirmCandidates()
{
  int i;
  for (i = 1; i <= edgedata.Size(); i++)
    {
      if (edgedata.Get(i).Candidate())
        {edgedata.Elem(i).SetConfirmed();}
    }
}

int STLEdgeDataList :: GetEdgeNum(int np1, int np2) const
{
  IVec<2> ed = IVec<2>(np1,np2).Sort();
  if (hashtab.Used(ed))
    {
      return hashtab.Get(ed);
    }

//   int i;
//   for (i = 1; i <= Size(); i++)
//     {
//       if ((Get(i).p1 == np1 && Get(i).p2 == np2) ||
//        (Get(i).p2 == np1 && Get(i).p1 == np2))
//      {
//        return i;
//      }
//     }

  return 0;
}

const STLEdgeDataList& STLEdgeDataList :: operator=(const STLEdgeDataList& edl)
{
  int i;
  SetSize(edl.Size());
  for (i = 1; i <= Size(); i++)
    {
      Add(edl.Get(i), i);
    }
  return *this;
} 

void STLEdgeDataList :: Add(const STLEdgeData& ed, int i)
{
  hashtab.Set(IVec<2>(ed.p1,ed.p2).Sort(), i);
  Elem(i) = ed;
  AddEdgePP(ed.p1,i);
  AddEdgePP(ed.p2,i);
}

void STLEdgeDataList :: Write(ofstream& of) const
{
  of.precision(16);
  int i;
  of << Size() << endl;
  
  for (i = 1; i <= Size(); i++)
    {
      Get(i).Write(of);
    }
}

void STLEdgeDataList :: Read(ifstream& ifs)
{
  int i,n;
  ifs >> n;

  SetSize(n);
  STLEdgeData ed;
  for (i = 1; i <= n; i++)
    {
      ed.Read(ifs);
      Add(ed,i);
    }
}

int STLEdgeDataList :: GetNEPPStat(int p, int status) const
{
  int i;
  int cnt = 0;
  for (i = 1; i <= GetNEPP(p); i++)
    {
      if (Get(GetEdgePP(p,i)).GetStatus() == status)
        {
          cnt++;
        }
    }
  return cnt;
}

int STLEdgeDataList :: GetNConfCandEPP(int p) const
{
  int i;
  int cnt = 0;
  for (i = 1; i <= GetNEPP(p); i++)
    {
      if (Get(GetEdgePP(p,i)).ConfCand())
        {
          cnt++;
        }
    }
  return cnt;
}


void STLEdgeDataList :: BuildLineWithEdge(int ep1, int ep2, Array<IVec<2>>& line)
{
  int status = Get(GetEdgeNum(ep1,ep2)).GetStatus();

  int found, pstart, p, en, pnew, ennew;
  int closed = 0;
  int j, i;
  for (j = 1; j <= 2; j++)
    {
      if (j == 1) {p = ep1;}
      if (j == 2) {p = ep2;}

      pstart = p;
      en = GetEdgeNum(ep1,ep2);

      found = 1;
      while (found && !closed)
        {
          found = 0;
          
          if (GetNEPPStat(p,status) == 2)
            {
              for (i = 1; i <= GetNEPP(p); i++)
                {               
                  const STLEdgeData& e = Get(GetEdgePP(p,i));
                  if (GetEdgePP(p,i) != en && e.GetStatus() == status) 
                    {
                      if (e.p1 == p) 
                        {pnew = e.p2;}
                      else 
                        {pnew = e.p1;}

                      ennew = GetEdgePP(p,i);
                    }
                }
              if (pnew == pstart) {closed = 1;}
              else
                {
                  line.Append(IVec<2>(p,pnew));
                  p = pnew;
                  en = ennew;
                  found = 1;
                }
            }
        }
    }
  
}
*/




STLEdgeDataList :: STLEdgeDataList (STLTopology & ageom)
  : geom(ageom)
{
  ;
}

STLEdgeDataList :: ~STLEdgeDataList()
{
  ;
}


void STLEdgeDataList :: Store ()
{
  int i, ne = geom.GetNTE();
  storedstatus.SetSize(ne);
  for (i = 1; i <= ne; i++)
    {
      storedstatus[i-1] = Get(i).GetStatus();
    }
}

void STLEdgeDataList :: Restore ()
{
  int i, ne = geom.GetNTE();
  if (storedstatus.Size() == ne)
    for (i = 1; i <= ne; i++)
      geom.GetTopEdge(i).SetStatus (storedstatus[i-1]);
}


void STLEdgeDataList :: ResetAll()
{
  int i, ne = geom.GetNTE();
  for (i = 1; i <= ne; i++)
    geom.GetTopEdge (i).SetStatus (ED_UNDEFINED);
}

int STLEdgeDataList :: GetNConfEdges() const
{
  int i, ne = geom.GetNTE();
  int cnt = 0;
  for (i = 1; i <= ne; i++)
    if (geom.GetTopEdge (i).GetStatus() == ED_CONFIRMED)
      cnt++;
  return cnt; 
}

void STLEdgeDataList :: ChangeStatus(int status1, int status2)
{
  int i, ne = geom.GetNTE();
  for (i = 1; i <= ne; i++)
    if (geom.GetTopEdge (i).GetStatus() == status1)
      geom.GetTopEdge (i).SetStatus (status2);
}

/*
void STLEdgeDataList :: Add(const STLEdgeData& ed, int i)
{
  IVec<2> edge(ed.p1,ed.p2);
  edge.Sort();
  hashtab.Set(edge, i);
  Elem(i) = ed;
  AddEdgePP(ed.p1,i);
  AddEdgePP(ed.p2,i);
}
*/

void STLEdgeDataList :: Write(ofstream& of) const
{
  
  /*
  of.precision(16);
  int i;
  of << Size() << endl;
  
  for (i = 1; i <= Size(); i++)
    {
      Get(i).Write(of);
    }

  */
  of.precision(16);
  int i, ne = geom.GetNTE();
  //of << GetNConfEdges() << endl;
  of << geom.GetNTE() << endl;

  for (i = 1; i <= ne; i++)
    {
      const STLTopEdge & edge = geom.GetTopEdge(i);
      //if (edge.GetStatus() == ED_CONFIRMED)
      of << edge.GetStatus() << " ";

      const Point<3> & p1 = geom.GetPoint (edge[0]);
      const Point<3> & p2 = geom.GetPoint (edge[1]);
      of << p1(0) << " "
         << p1(1) << " "
         << p1(2) << " "
         << p2(0) << " "
         << p2(1) << " "
         << p2(2) << endl;
    }
  
}

void STLEdgeDataList :: Read(ifstream& ifs)
{
  int i, nce;
  Point<3> p1, p2;
  int pi1, pi2;
  int status, ednum;

  ifs >> nce;
  for (i = 1; i <= nce; i++)
    {
      ifs >> status;
      ifs >> p1(0) >> p1(1) >> p1(2);
      ifs >> p2(0) >> p2(1) >> p2(2);

      pi1 = geom.GetPointNum (p1);
      pi2 = geom.GetPointNum (p2);
      ednum = geom.GetTopEdgeNum (pi1, pi2);


      if (ednum)
        { 
          geom.GetTopEdge(ednum).SetStatus (status);
        //      geom.GetTopEdge (ednum).SetStatus (ED_CONFIRMED);
        }
    }
    /*
  int i,n;
  ifs >> n;

  SetSize(n);
  STLEdgeData ed;
  for (i = 1; i <= n; i++)
    {
      ed.Read(ifs);
      Add(ed,i);
    }
  */
}

int STLEdgeDataList :: GetNEPPStat(int p, int status) const
{
  int i;
  int cnt = 0;
  for (i = 1; i <= GetNEPP(p); i++)
    {
      if (Get(GetEdgePP(p,i)).GetStatus() == status)
        {
          cnt++;
        }
    }
  return cnt;
}

int STLEdgeDataList :: GetNConfCandEPP(int p) const
{
  int i;
  int cnt = 0;
  for (i = 1; i <= GetNEPP(p); i++)
    {
      if (Get(GetEdgePP(p,i)).GetStatus() == ED_CANDIDATE || 
          Get(GetEdgePP(p,i)).GetStatus() == ED_CONFIRMED)
        {
          cnt++;
        }
    }
  return cnt;
}


void STLEdgeDataList :: BuildLineWithEdge(int ep1, int ep2, Array<IVec<2>>& line)
{
  int status = Get(GetEdgeNum(ep1,ep2)).GetStatus();

  int found, pstart, p(0), en, pnew(0), ennew(0);
  int closed = 0;
  int j, i;
  for (j = 1; j <= 2; j++)
    {
      if (j == 1) {p = ep1;}
      if (j == 2) {p = ep2;}

      pstart = p;
      en = GetEdgeNum(ep1,ep2);

      found = 1;
      while (found && !closed)
        {
          found = 0;
          
          if (GetNEPPStat(p,status) == 2)
            {
              for (i = 1; i <= GetNEPP(p); i++)
                {               
                  const STLTopEdge & e = Get(GetEdgePP(p,i));
                  if (GetEdgePP(p,i) != en && e.GetStatus() == status) 
                    {
                      if (e[0] == p) 
                        {pnew = e[1];}
                      else 
                        {pnew = e[0];}

                      ennew = GetEdgePP(p,i);
                    }
                }
              if (pnew == pstart) {closed = 1;}
              else
                {
                  line.Append(IVec<2>(p,pnew));
                  p = pnew;
                  en = ennew;
                  found = 1;
                }
            }
        }
    }
  
}

int Exists(int p1, int p2, const Array<IVec<2>>& line)
{
  int i;
  for (i = 1; i <= line.Size(); i++)
    {
      if ( (line[i-1][0] == p1 && line[i-1][1] == p2) ||
           (line[i-1][0] == p2 && line[i-1][1] == p1) )
        {return 1;}
    }
  return 0;
}

void STLEdgeDataList :: BuildClusterWithEdge(int ep1, int ep2, Array<IVec<2>>& line)
{
  int status = Get(GetEdgeNum(ep1,ep2)).GetStatus();

  int p(0), en;
  int j, i, k;
  int oldend;
  int newend = 1;
  STLPointId pnew;
  int ennew(0);

  int changed = 1;
  while (changed)
    {
      changed = 0;
      for (j = 1; j <= 2; j++)
        {
          oldend = newend;
          newend = line.Size();
          for (k = oldend; k <= line.Size(); k++)
            {
              if (j == 1) p = line[k-1][0];
              if (j == 2) p = line[k-1][1];
              en = GetEdgeNum(line[k-1][0], line[k-1][1]);

              for (i = 1; i <= GetNEPP(p); i++)
                {               
                  pnew = 0;
                  const STLTopEdge & e = Get(GetEdgePP(p,i));
                  if (GetEdgePP(p,i) != en && e.GetStatus() == status) 
                    {
                      if (e[0] == p) 
                        {pnew = e[1];}
                      else 
                        {pnew = e[0];}

                      ennew = GetEdgePP(p,i);
                    }
                  if (pnew && !Exists(p,pnew,line))
                    {
                      changed = 1;
                      line.Append(IVec<2>(p,pnew));
                      p = pnew;
                      en = ennew;
                    }
                }
              
            }
        }

    }

}










//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++   STL LINE    +++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

STLLine :: STLLine(const STLGeometry * ageometry)
  : pts(), lefttrigs(), righttrigs()
{
  geometry = ageometry;
  split = 0;
};

int STLLine :: GetNS() const
{
  if (pts.Size() <= 1) {return 0;}
  return pts.Size()-1;
}
void STLLine :: GetSeg(int nr, int& p1, int& p2) const
{
  p1 = pts[nr-1];
  p2 = pts[nr];
}

int STLLine :: GetLeftTrig(int nr) const 
{
  if (nr > lefttrigs.Size()) {PrintSysError("In STLLine::GetLeftTrig!!!"); return 0;}
  return lefttrigs[nr-1];
};

int STLLine :: GetRightTrig(int nr) const 
{
  if (nr > righttrigs.Size()) {PrintSysError("In STLLine::GetRightTrig!!!"); return 0;}
  return righttrigs[nr-1];
};

double STLLine :: GetSegLen(const Array<Point<3>,STLPointId>& ap, int nr) const
{
  return Dist(ap[PNum(nr)],ap[PNum(nr+1)]);
}

double STLLine :: GetLength(const Array<Point<3>,STLPointId>& ap) const
{
  double len = 0;
  for (int i = 2; i <= pts.Size(); i++)
    len += (ap[pts[i-1]] - ap[pts[i-2]]).Length();
  return len;
}

void STLLine :: GetBoundingBox (const Array<Point<3>,STLPointId> & ap, Box<3> & box) const
{
  box.Set (ap[pts[0]]);
  for (int i = 1; i < pts.Size(); i++)
    box.Add (ap[pts[i]]);
}



Point<3> STLLine :: 
GetPointInDist(const Array<Point<3>,STLPointId>& ap, double dist, int& index) const
{
  if (dist <= 0)
    {
      index = 1;
      return ap[StartP()];
    }
  
  double len = 0;
  int i;
  for (i = 1; i < pts.Size(); i++)
    {
      double seglen = Dist (ap[pts[i-1]],
                            ap[pts[i]]);

      if (len + seglen > dist)
        {
          index = i;
          double relval = (dist - len) / (seglen + 1e-16);
          Vec<3> v (ap[pts[i-1]], ap[pts[i]]);
          return ap[pts[i-1]] + relval * v;
        }

      len += seglen;
    }

  index = pts.Size() - 1;
  return ap[EndP()];
}


/*
double stlgh;
double GetH(const Point<3>& p, double x) 
{
  return stlgh;//+0.5)*(x+0.5);
}
*/
STLLine* STLLine :: Mesh(const Array<Point<3>,STLPointId>& ap, 
                         Array<Point<3>>& mp, double ghi,
                         class Mesh& mesh) const
{
  static Timer timer1a("mesh stl-line 1a");
  static Timer timer1b("mesh stl-line 1b");
  static Timer timer2("mesh stl-line 2");
  static Timer timer3("mesh stl-line 3");

  timer1a.Start();

  STLLine* line = new STLLine(geometry);

  //stlgh = ghi; //uebergangsloesung!!!!
  
  double len = GetLength(ap);
  double inthl = 0; //integral of 1/h
  double dist = 0;
  double h;
  int ind;
  Point<3> p;

  Box<3> bbox;
  GetBoundingBox (ap, bbox);
  double diam = bbox.Diam();

  double minh = mesh.LocalHFunction().GetMinH (bbox.PMin(), bbox.PMax());

  double maxseglen = 0;
  for (int i = 1; i <= GetNS(); i++)
    maxseglen = max2 (maxseglen, GetSegLen (ap, i));
  
  int nph = 10+int(maxseglen / minh); //anzahl der integralauswertungen pro segment

  Array<double> inthi(GetNS()*nph);
  Array<double> curvelen(GetNS()*nph);

  timer1a.Stop();
  timer1b.Start();


  for (int i = 1; i <= GetNS(); i++)
    {
      //double seglen = GetSegLen(ap,i);
      for (int j = 1; j <= nph; j++)
        {
          p = GetPointInDist(ap,dist,ind);
          //h = GetH(p,dist/len);
          h = mesh.GetH(p);

          
          dist += GetSegLen(ap,i)/(double)nph;
          
          inthl += GetSegLen(ap,i)/nph/(h);
          inthi[(i-1)*nph+j-1] = GetSegLen(ap,i)/nph/h;
          curvelen[(i-1)*nph+j-1] = GetSegLen(ap,i)/nph;
        }
    }


  int inthlint = int(inthl+1);

  if ( (inthlint < 3) && (StartP() == EndP()))
    {
      inthlint = 3;
    }
  if ( (inthlint == 1) && ShouldSplit())
    {
      inthlint = 2; 
    }
     
  double fact = inthl/(double)inthlint;
  dist = 0;
  int j = 1;


  p = ap[StartP()];
  int pn = AddPointIfNotExists(mp, p, 1e-10*diam);

  int segn = 1;
  line->AddPoint(pn);
  line->AddLeftTrig(GetLeftTrig(segn));
  line->AddRightTrig(GetRightTrig(segn));
  line->AddDist(dist);

  timer1b.Stop();
  timer2.Start();

  inthl = 0; //restart each meshseg
  for (int i = 1; i <= inthlint; i++)
    {
      while (inthl < 1.000000001 && j <= inthi.Size())
        {
          inthl += inthi[j-1]/fact;
          dist += curvelen[j-1];
          j++;
        }

      //went too far:
      j--;
      double tofar = (inthl - 1)/inthi[j-1];
      inthl -= tofar*inthi[j-1];
      dist -= tofar*curvelen[j-1]*fact;

      if (i == inthlint && fabs(dist - len) >= 1E-8) 
        {
          PrintSysError("meshline failed!!!"); 
        }

      if (i != inthlint) 
        {
          p = GetPointInDist(ap,dist,ind);
          pn = AddPointIfNotExists(mp, p, 1e-10*diam);
          segn = ind;
          line->AddPoint(pn);
          line->AddLeftTrig(GetLeftTrig(segn));
          line->AddRightTrig(GetRightTrig(segn));
          line->AddDist(dist);
        }

      inthl = tofar*inthi[j-1];
      dist += tofar*curvelen[j-1]*fact;
      j++;
    }

  timer2.Stop();
  timer3.Start();


  p = ap[EndP()];
  pn = AddPointIfNotExists(mp, p, 1e-10*diam);
  segn = GetNS();
  line->AddPoint(pn);
  line->AddLeftTrig(GetLeftTrig(segn));
  line->AddRightTrig(GetRightTrig(segn));
  line->AddDist(dist);
  
  for (int ii = 1; ii <= line->GetNS(); ii++)
    {
      int p1, p2;
      line->GetSeg(ii,p1,p2);
    }
  /*  
  (*testout) << "line, " << ap.Get(StartP()) << "-" << ap.Get(EndP())
             << " len = " << Dist (ap.Get(StartP()), ap.Get(EndP())) << endl;
  */

  timer3.Stop();

  return line;
}
}
