#include <mystdlib.h>
#include "meshing.hpp"

#ifdef WIN32
#define COMMASIGN ':'
#else
#define COMMASIGN ','
#endif


namespace netgen
{

extern const char * tetrules[];

void LoadVMatrixLine (istream & ist, DenseMatrix & m, int line)
{
  char ch;
  int pnum;
  float f;
  
  ist >> ch;
  while (ch != '}')
    {
      ist.putback (ch);
      ist >> f;
      ist >> ch;
      ist >> pnum;
      
      if (ch == 'x' || ch == 'X')
	m.Elem(line, 3 * pnum - 2) = f;
      if (ch == 'y' || ch == 'Y')
	m.Elem(line, 3 * pnum - 1) = f;
      if (ch == 'z' || ch == 'Z')
	m.Elem(line, 3 * pnum    ) = f;

      if (ch == 'p' || ch == 'P')
	{
	  m.Elem(line  , 3 * pnum-2) = f;
	  m.Elem(line+1, 3 * pnum-1) = f;
	  m.Elem(line+2, 3 * pnum  ) = f;
	}

      ist >> ch;
      if (ch == COMMASIGN)
	ist >> ch;
    }
}





int vnetrule :: NeighbourTrianglePoint (const threeint & t1, const threeint & t2) const
{
  NgArray<int> tr1(3);
  NgArray<int> tr2(3);
  tr1[0]=t1.i1;
  tr1[1]=t1.i2;
  tr1[2]=t1.i3;
  tr2[0]=t2.i1;
  tr2[1]=t2.i2;
  tr2[2]=t2.i3;


  int ret=0;

  for (int i=1; i<=3; i++)
    {
      for (int j=1; j<=3; j++)
	{
	  if ((tr1[i-1]==tr2[j-1] && tr1[(i%3)]==tr2[(j%3)]) ||
              (tr1[i-1]==tr2[(j%3)] && tr1[(i%3)]==tr2[j-1]))
	    {ret = tr2[(j+1)%3];}
	}      
    }

  return ret;

}

void vnetrule :: LoadRule (istream & ist)
{
  char buf[256];
  char ch, ok;
  Point3d p;
  RuleElement2d face(3);
  int i, j, i1, i2, i3, fs, ii, ii1, ii2, ii3;
  twoint edge;
  DenseMatrix tempoldutonewu(30, 20), 
    tempoldutofreezone(30, 20),
    tempoldutofreezonelimit(30, 20),
    tfz(20, 20),
    tfzl(20, 20);

  tempoldutonewu = 0;
  tempoldutofreezone = 0;
  tfz = 0;
  tfzl = 0;


  noldp = 0;
  noldf = 0;

  ist.get (buf, sizeof(buf), '"');
  ist.get (ch);
  ist.get (buf, sizeof(buf), '"');
  ist.get (ch);

  delete [] name;
  name = new char[strlen (buf) + 1];
  strcpy (name, buf);
  //  (*mycout) << "Rule " << name << " found." << endl;

  do
    {
      ist >> buf;

      if (strcmp (buf, "quality") == 0)

	{
	  ist >> quality;
	}

      else if (strcmp (buf, "flags") == 0)
	{
	  ist >> ch;
	  while (ch != ';')
	    {
	      flags.Append (ch);
	      ist >> ch;
	    }
	}

      else if (strcmp (buf, "mappoints") == 0)
	{
	  ist >> ch;

	  while (ch == '(')
	    {
	      ist >> p.X();
	      ist >> ch;    // ','
	      ist >> p.Y();
	      ist >> ch;    // ','
	      ist >> p.Z();
	      ist >> ch;    // ')'

	      points.Append (p);
	      noldp++;

	      tolerances.SetSize (noldp);
	      tolerances[tolerances.Range().Next()-1] = 1;

	      ist >> ch;
	      while (ch != ';')
		{
		  if (ch == '{')
		    {
		      ist >> tolerances[tolerances.Range().Next()-1];
		      ist >> ch;  // '}'
		    }

		  ist >> ch;
		}

	      ist >> ch;
	    }

	  ist.putback (ch);
	}


      else if (strcmp (buf, "mapfaces") == 0)
	{
	  ist >> ch;

	  while (ch == '(')
	    {
	      face.SetNP(3);
	      ist >> face.PNum(1);
	      ist >> ch;    // ','
	      ist >> face.PNum(2);
	      ist >> ch;    // ','
	      ist >> face.PNum(3);
	      ist >> ch;    // ')' or ','
	      if (ch == COMMASIGN)
		{
		  face.SetNP(4);
		  ist >> face.PNum(4);
		  ist >> ch;    // ')' 
		}
	      faces.Append (face);
	      noldf++;

	      ist >> ch;
	      while (ch != ';')
		{
		  if (ch == 'd')
		    {
		      delfaces.Append (noldf);
		      ist >> ch; // 'e'
		      ist >> ch; // 'l'
		    }

		  ist >> ch;
		}

	      ist >> ch;
	    }

	  ist.putback (ch);
	}

      else if (strcmp (buf, "mapedges") == 0)
	{
	  ist >> ch;

	  while (ch == '(')
	    {
	      ist >> edge.i1;
	      ist >> ch;    // ','
	      ist >> edge.i2;
	      ist >> ch;    // ')'

	      edges.Append (edge);

	      ist >> ch;
	      while (ch != ';')
		{
		  ist >> ch;
		}

	      ist >> ch;
	    }

	  ist.putback (ch);
	}


      else if (strcmp (buf, "newpoints") == 0)
	{
	  ist >> ch;

	  while (ch == '(')
	    {
	      ist >> p.X();
	      ist >> ch;    // ','
	      ist >> p.Y();
	      ist >> ch;    // ','
	      ist >> p.Z();
	      ist >> ch;    // ')'

	      points.Append (p);

	      ist >> ch;
	      while (ch != ';')
		{
		  if (ch == '{')
		    {
		      LoadVMatrixLine (ist, tempoldutonewu,
				       3 * (points.Size()-noldp) - 2);

		      ist >> ch; // '{'
		      LoadVMatrixLine (ist, tempoldutonewu,
				       3 * (points.Size()-noldp) - 1);

		      ist >> ch; // '{'
		      LoadVMatrixLine (ist, tempoldutonewu,
				       3 * (points.Size()-noldp)    );
		    }

		  ist >> ch;
		}

	      ist >> ch;
	    }

	  ist.putback (ch);
	}

      else if (strcmp (buf, "newfaces") == 0)
	{
	  ist >> ch;

	  while (ch == '(')
	    {
	      face.SetNP(3);
	      ist >> face.PNum(1);
	      ist >> ch;    // ','
	      ist >> face.PNum(2);
	      ist >> ch;    // ','
	      ist >> face.PNum(3);
	      ist >> ch;    // ')' or ','
	      if (ch == COMMASIGN)
		{
		  face.SetNP(4);
		  ist >> face.PNum(4);
		  ist >> ch;    // ')' 
		}
	      faces.Append (face);

	      ist >> ch;
	      while (ch != ';')
		{
		  ist >> ch;
		}

	      ist >> ch;
	    }

	  ist.putback (ch);
	}

      else if (strcmp (buf, "freezone") == 0)
	{
	  ist >> ch;
	
	  while (ch == '(')
	    {
	      ist >> p.X();
	      ist >> ch;    // ','
	      ist >> p.Y();
	      ist >> ch;    // ','
	      ist >> p.Z();
	      ist >> ch;    // ')'
	    
	      freezone.Append (p);
	    
	      ist >> ch;
	      while (ch != ';')
		{
		  if (ch == '{')
		    {
		      LoadVMatrixLine (ist, tempoldutofreezone,
				       3 * freezone.Size() - 2);
		    
		      ist >> ch; // '{'
		      LoadVMatrixLine (ist, tempoldutofreezone,
				       3 * freezone.Size() - 1);
		    
		      ist >> ch; // '{'
		      LoadVMatrixLine (ist, tempoldutofreezone,
				       3 * freezone.Size()    );
		    }
		
		  ist >> ch;
		}
	    
	      ist >> ch;
	    }
	
	  ist.putback (ch);
	}
      else if (strcmp (buf, "freezone2") == 0)
	{
	  int k, nfp;

	  nfp = 0;
	  ist >> ch;

	  DenseMatrix hm1(3, 50), hm2(50, 50), hm3(50, 50);
	  hm3 = 0;

	  while (ch == '{')
	    {
	      hm1 = 0;
	      nfp++;
	      LoadVMatrixLine (ist, hm1, 1);

	      for (i = 1; i <= points.Size(); i++)
		tfz.Elem(nfp, i) = hm1.Get(1, 3*i-2);


	      p.X() = p.Y() = p.Z() = 0;
	      for (auto pi : points.Range())
		{
		  int i = pi.Nr1();
		  p.X() += hm1.Get(1, 3*i-2) * points[pi].X();
		  p.Y() += hm1.Get(1, 3*i-2) * points[pi].Y();
		  p.Z() += hm1.Get(1, 3*i-2) * points[pi].Z();
		}
	      freezone.Append (p);
	      freezonelimit.Append (p);
	    
	      hm2 = 0;
	      for (i = 1; i <= 3 * noldp; i++)
		hm2.Elem(i, i) = 1;
	      for (i = 1; i <= 3 * noldp; i++)
		for (j = 1; j <= 3 * (points.Size() - noldp); j++)
		  hm2.Elem(j + 3 * noldp, i) = tempoldutonewu.Get(j, i);
		  
	      for (i = 1; i <= 3; i++)
		for (j = 1; j <= 3 * noldp; j++)
		  {
		    double sum = 0;
		    for (k = 1; k <= 3 * points.Size(); k++)
		      sum += hm1.Get(i, k) * hm2.Get(k, j);
		  
		    hm3.Elem(i + 3 * (nfp-1), j) = sum;
		  }

	      //	    (*testout) << "freepoint: " << p << endl;

	      while (ch != ';')
		ist >> ch; 

	      ist >> ch;
	    }

	  tfzl = tfz;

	  tempoldutofreezone = hm3;
	  tempoldutofreezonelimit = hm3;
	  ist.putback(ch);
	}

      else if (strcmp (buf, "freezonelimit") == 0)
	{
	  int k, nfp;
	  nfp = 0;
	  ist >> ch;

	  DenseMatrix hm1(3, 50), hm2(50, 50), hm3(50, 50);
	  hm3 = 0;

	  while (ch == '{')
	    {
	      hm1 = 0;
	      nfp++;
	      LoadVMatrixLine (ist, hm1, 1);

	      for (i = 1; i <= points.Size(); i++)
		tfzl.Elem(nfp, i) = hm1.Get(1, 3*i-2);


	      p.X() = p.Y() = p.Z() = 0;
	      for (auto pi : points.Range())
		{
		  int i = pi.Nr1();
		  p.X() += hm1.Get(1, 3*i-2) * points[pi].X();
		  p.Y() += hm1.Get(1, 3*i-2) * points[pi].Y();
		  p.Z() += hm1.Get(1, 3*i-2) * points[pi].Z();
		}
	      freezonelimit[nfp-1] = p;
	    
	      hm2 = 0;
	      for (i = 1; i <= 3 * noldp; i++)
		hm2.Elem(i, i) = 1;
	      for (i = 1; i <= 3 * noldp; i++)
		for (j = 1; j <= 3 * (points.Size() - noldp); j++)
		  hm2.Elem(j + 3 * noldp, i) = tempoldutonewu.Get(j, i);
		  
	      for (i = 1; i <= 3; i++)
		for (j = 1; j <= 3 * noldp; j++)
		  {
		    double sum = 0;
		    for (k = 1; k <= 3 * points.Size(); k++)
		      sum += hm1.Get(i, k) * hm2.Get(k, j);
		  
		    hm3.Elem(i + 3 * (nfp-1), j) = sum;
		  }

	      //	    (*testout) << "freepoint: " << p << endl;

	      while (ch != ';')
		ist >> ch; 

	      ist >> ch;
	    }

	  tempoldutofreezonelimit = hm3;
	  ist.putback(ch);
	}

      else if (strcmp (buf, "freeset") == 0)
	{
	  freesets.Append (new NgArray<int>);

	  ist >> ch;

	  while (ch != ';')
	    {
	      ist.putback (ch);
	      ist >> i;
	      freesets.Last()->Append(i);
	      ist >> ch;
	    }
	}

      else if (strcmp (buf, "elements") == 0)
	{
	  ist >> ch;

	  while (ch == '(')
	    {
	      elements.Append (RuleElement(4));

	      //	      elements.Last().SetNP(1);
	      ist >> elements.Last().PNum(1);
	      ist >> ch;    // ','

	      if (ch == COMMASIGN)
		{
		  //		  elements.Last().SetNP(2);
		  ist >> elements.Last().PNum(2);
		  ist >> ch;    // ','
		}
	      if (ch == COMMASIGN)
		{
		  //		  elements.Last().SetNP(3);
		  ist >> elements.Last().PNum(3);
		  ist >> ch;    // ','
		}
	      if (ch == COMMASIGN)
		{
		  //		  elements.Last().SetNP(4);
		  elements.Last().SetType(TET);
		  ist >> elements.Last().PNum(4);
		  ist >> ch;    // ','
		}
	      if (ch == COMMASIGN)
		{
		  //		  elements.Last().SetNP(5);
		  elements.Last().SetType(PYRAMID);
		  ist >> elements.Last().PNum(5);
		  ist >> ch;    // ','
		}
	      if (ch == COMMASIGN)
		{
		  //		  elements.Last().SetNP(6);
		  elements.Last().SetType(PRISM);
		  ist >> elements.Last().PNum(6);
		  ist >> ch;    // ','
		}
              
	      if (ch == COMMASIGN)
		{
		  //		  elements.Last().SetNP(6);
		  elements.Last().SetType(HEX);
		  ist >> elements.Last().PNum(7);
		  ist >> ch;    // ','
		}
	      if (ch == COMMASIGN)
		{
		  //		  elements.Last().SetNP(6);
		  elements.Last().SetType(HEX);
		  ist >> elements.Last().PNum(8);
		  ist >> ch;    // ','
		}

	      /*
	      orientations.Append (fourpoints());
	      orientations.Last().i1 = elements.Last().PNum(1);
	      orientations.Last().i2 = elements.Last().PNum(2);
	      orientations.Last().i3 = elements.Last().PNum(3);
	      orientations.Last().i4 = elements.Last().PNum(4);
	      */

	      ist >> ch;
	      while (ch != ';')
		{
		  ist >> ch;
		}

	      ist >> ch;
	    }

	  ist.putback (ch);
	}

      else if (strcmp (buf, "orientations") == 0)

	{
	  ist >> ch;

	  while (ch == '(')
	    {
	      //        fourint a = fourint();
	      orientations.Append (fourpoints());

	      ist >> orientations.Last().i1;
	      ist >> ch;    // ','
	      ist >> orientations.Last().i2;
	      ist >> ch;    // ','
	      ist >> orientations.Last().i3;
	      ist >> ch;    // ','
	      ist >> orientations.Last().i4;
	      ist >> ch;    // ','


	      ist >> ch;
	      while (ch != ';')
		{
		  ist >> ch;
		}

	      ist >> ch;
	    }

	  ist.putback (ch);
	}


      else if (strcmp (buf, "endrule") != 0)
	{
	  PrintSysError ("Parser3d, unknown token " , buf);
	}
    }
  while (!ist.eof() && strcmp (buf, "endrule") != 0);


  //  (*testout) << endl;
  //  (*testout) << Name() << endl;
  //  (*testout) << "no1 = " << GetNO() << endl;

  oldutonewu.SetSize (3 * (points.Size() - noldp), 3 * noldp);
  oldutonewu = 0;

  for (i = 1; i <= oldutonewu.Height(); i++)
    for (j = 1; j <= oldutonewu.Width(); j++)
      oldutonewu.Elem(i, j) = tempoldutonewu.Elem(i, j);


  /*
    oldutofreezone = new SparseMatrixFlex (3 * freezone.Size(), 3 * noldp);
    oldutofreezonelimit = new SparseMatrixFlex (3 * freezone.Size(), 3 * noldp);

    oldutofreezone -> SetSymmetric(0);
    oldutofreezonelimit -> SetSymmetric(0);
    */

  /*
    oldutofreezone = new DenseMatrix (3 * freezone.Size(), 3 * noldp);
    oldutofreezonelimit = new DenseMatrix (3 * freezone.Size(), 3 * noldp);
  
    for (i = 1; i <= oldutofreezone->Height(); i++)
    for (j = 1; j <= oldutofreezone->Width(); j++)
    //      if (j == 4 || j >= 7)
    {
    if (tempoldutofreezone.Elem(i, j))
    (*oldutofreezone)(i, j) = tempoldutofreezone(i, j);
    if (tempoldutofreezonelimit.Elem(i, j))
    (*oldutofreezonelimit)(i, j) = tempoldutofreezonelimit(i, j);
    }
    */




  oldutofreezone = new DenseMatrix (freezone.Size(), points.Size());
  oldutofreezonelimit = new DenseMatrix (freezone.Size(), points.Size());
  //  oldutofreezone = new SparseMatrixFlex (freezone.Size(), points.Size());
  //  oldutofreezonelimit = new SparseMatrixFlex (freezone.Size(), points.Size());

  for (i = 1; i <= freezone.Size(); i++)
    for (j = 1; j <= points.Size(); j++)
      {
	if (tfz.Elem(i, j))
	  (*oldutofreezone).Elem(i, j) = tfz.Elem(i, j);
	if (tfzl.Elem(i, j))
	  (*oldutofreezonelimit).Elem(i, j) = tfzl.Elem(i, j);
      }
  
  /*
  (*testout) << "Rule " << Name() << endl;
  (*testout) << "oldutofreezone = " << (*oldutofreezone) << endl;
  (*testout) << "oldutofreezonelimit = " << (*oldutofreezonelimit) << endl;
  */

  freezonepi.SetSize (freezone.Size());
  for (i = 1; i <= freezonepi.Size(); i++)
    freezonepi[i-1] = 0;
  for (i = 1; i <= freezone.Size(); i++)
    for (auto pj : points.Range().Modify(0, noldp-points.Size()))
      if (Dist (freezone[i-1], points[pj]) < 1e-8)
	freezonepi[i-1] = pj.Nr1();



  
  for (i = 1; i <= elements.Size(); i++)
    {
      if (elements[i-1].GetNP() == 4)
	{
	  orientations.Append (fourpoints());
	  orientations.Last().i1 = elements[i-1].PNum(1);
	  orientations.Last().i2 = elements[i-1].PNum(2);
	  orientations.Last().i3 = elements[i-1].PNum(3);
	  orientations.Last().i4 = elements[i-1].PNum(4);
	}
      if (elements[i-1].GetNP() == 5)
	{
	  orientations.Append (fourpoints());
	  orientations.Last().i1 = elements[i-1].PNum(1);
	  orientations.Last().i2 = elements[i-1].PNum(2);
	  orientations.Last().i3 = elements[i-1].PNum(3);
	  orientations.Last().i4 = elements[i-1].PNum(5);

	  orientations.Append (fourpoints());
	  orientations.Last().i1 = elements[i-1].PNum(1);
	  orientations.Last().i2 = elements[i-1].PNum(3);
	  orientations.Last().i3 = elements[i-1].PNum(4);
	  orientations.Last().i4 = elements[i-1].PNum(5);
	}
    }



  if (freesets.Size() == 0)
    {
      freesets.Append (new NgArray<int>);
      for (i = 1; i <= freezone.Size(); i++)
	freesets[0]->Append(i);
    }


  //  testout << "Freezone: " << endl;

  //  for (i = 1; i <= freezone.Size(); i++)
  //    (*testout) << "freepoint: " << freezone.Get(i) << endl;
  Vector vp(points.Size()), vfp(freezone.Size());


  if (quality < 100)
    {
      for (int i = 1; i <= 3; i++)
	{
	  for (auto pj : points.Range())
	    vp(pj.Nr0()) = points[pj].X(i);
	  oldutofreezone->Mult(vp, vfp);
	  for (int j = 1; j <= freezone.Size(); j++)
	    freezone[j-1].X(i) = vfp(j-1);
	}
      //      for (i = 1; i <= freezone.Size(); i++)
      //	(*testout) << "freepoint: " << freezone.Get(i) << endl;
    }


  for (fs = 1; fs <= freesets.Size(); fs++)
    {
      freefaces.Append (new NgArray<threeint>);

      NgArray<int> & freeset = *freesets[fs-1];
      NgArray<threeint> & freesetfaces = *freefaces.Last();

      for (ii1 = 1; ii1 <= freeset.Size(); ii1++)
	for (ii2 = 1; ii2 <= freeset.Size(); ii2++)
	  for (ii3 = 1; ii3 <= freeset.Size(); ii3++)
	    if (ii1 < ii2 && ii1 < ii3 && ii2 != ii3)
	      {
		i1 = freeset[ii1-1];
		i2 = freeset[ii2-1];
		i3 = freeset[ii3-1];

		Vec3d v1, v2, n;

		v1 = freezone[i3-1] - freezone[i1-1];
		v2 = freezone[i2-1] - freezone[i1-1];
		n = Cross (v1, v2);
		n /= n.Length();
		//		(*testout) << "i1,2,3 = " << i1 << ", " << i2 << ", " << i3 << endl;
		//		(*testout) << "v1 = " << v1 << " v2 = " << v2 << " n = " << n << endl;
		ok = 1;
		for (ii = 1; ii <= freeset.Size(); ii++)
		  {
		    i = freeset[ii-1];
		    //		    (*testout) << "i = " << i << endl;
		    if (i != i1 && i != i2 && i != i3)
		      if ( (freezone[i-1] - freezone[i1-1]) * n < 0 ) ok = 0;
		  }

		if (ok)
		  {
		    freesetfaces.Append (threeint());
		    freesetfaces.Last().i1 = i1;
		    freesetfaces.Last().i2 = i2;
		    freesetfaces.Last().i3 = i3;
		  }
	      }
    }

  for (fs = 1; fs <= freesets.Size(); fs++)
    {
      freefaceinequ.Append (new DenseMatrix (freefaces[fs-1]->Size(), 4));
    }


  {
    int minn;
    //    NgArray<int> pnearness (noldp);
    pnearness.SetSize (noldp);

    pnearness = INT_MAX/10;

    for (j = 1; j <= GetNP(1); j++)
      pnearness[GetPointNr (1, j)] = 0;

    do
      {
	ok = 1;

	for (i = 1; i <= noldf; i++)
	  {
	    minn = INT_MAX/10;
	    for (j = 1; j <= GetNP(i); j++)
	      minn = min2 (minn, pnearness[GetPointNr (i, j)]);

	    for (j = 1; j <= GetNP(i); j++)
	      if (pnearness[GetPointNr (i, j)] > minn+1)
		{
		  ok = 0;
		  pnearness[GetPointNr (i, j)] = minn+1;
		}
	  }

	for (i = 1; i <= edges.Size(); i++)
	  {
	    RulePointIndex pi1 = RuleP(edges[i-1].i1);
	    RulePointIndex pi2 = RuleP(edges[i-1].i2);

	    if (pnearness[pi1] > pnearness[pi2]+1)
	      {
		ok = 0;
		pnearness[pi1] = pnearness[pi2]+1;
	      }
	    if (pnearness[pi2] > pnearness[pi1]+1)
	      {
		ok = 0;
		pnearness[pi2] = pnearness[pi1]+1;
	      }
	  }
	

	for (i = 1; i <= elements.Size(); i++)
	  if (elements[i-1].GetNP() == 6)  // prism rule
	    {
	      for (j = 1; j <= 3; j++)
		{
		  RulePointIndex pi1 = elements[i-1].PNum(j);
		  RulePointIndex pi2 = elements[i-1].PNum(j+3);

		  if (pnearness[pi1] > pnearness[pi2]+1)
		    {
		      ok = 0;
		      pnearness[pi1] = pnearness[pi2]+1;
		    }
		  if (pnearness[pi2] > pnearness[pi1]+1)
		    {
		      ok = 0;
		      pnearness[pi2] = pnearness[pi1]+1;
		    }
		}
	    }
      }
    while (!ok);

    maxpnearness = 0;
    for (auto pi : pnearness.Range())
      maxpnearness = max2 (maxpnearness, pnearness[pi]);


    fnearness.SetSize (noldf);

    for (i = 1; i <= noldf; i++)
      {
	fnearness[i-1] = 0;
	for (j = 1; j <= GetNP(i); j++)
	  fnearness[i-1] += pnearness[GetPointNr (i, j)];
      }

    // (*testout) << "rule " << name << ", pnear = " << pnearness << endl;
  }

  
  //Table of edges:
  for (fs = 1; fs <= freesets.Size(); fs++)
    {
      freeedges.Append (new NgArray<twoint>);
      
      //      NgArray<int> & freeset = *freesets.Get(fs);
      NgArray<twoint> & freesetedges = *freeedges.Last();
      NgArray<threeint> & freesetfaces = *freefaces[fs-1];
      // int k,l;
      // INDEX ind;
      
      for (int k = 1; k <= freesetfaces.Size(); k++)
	{
          // threeint tr = freesetfaces.Get(k);

	  for (int l = k+1; l <= freesetfaces.Size(); l++)
	    {
	      INDEX ind = NeighbourTrianglePoint(freesetfaces[k-1], freesetfaces[l-1]);
	      if (!ind) continue;

	      INDEX_3 f1(freesetfaces[k-1].i1, 
			 freesetfaces[k-1].i2, 
			 freesetfaces[k-1].i3);
	      INDEX_3 f2(freesetfaces[l-1].i1, 
			 freesetfaces[l-1].i2, 
			 freesetfaces[l-1].i3);
	      IVec<2,RulePointIndex> ed(RulePointIndex::INVALID, RulePointIndex::INVALID);
	      for (int f11 = 1; f11 <= 3; f11++)
		for (int f12 = 1; f12 <= 3; f12++)
		  if (f11 != f12)
		    for (int f21 = 1; f21 <= 3; f21++)
		      for (int f22 = 1; f22 <= 3; f22++)		    
			if (f1.I(f11) == f2.I(f21) && f1.I(f12) == f2.I(f22))
			{
			  ed[0] = RuleP(f1.I(f11));
			  ed[1] = RuleP(f1.I(f12));
			}
	      //	      (*testout) << "ed = " << ed.I(1) << "-" << ed.I(2) << endl;
	      //	      (*testout) << "ind = " << ind << " ed = " << ed << endl;
	      for (int eli = 1; eli <= GetNOldF(); eli++)
		{
		  if (GetNP(eli) == 4)
		    {
		      for (int elr = 1; elr <= 4; elr++)
			{
			  if (GetPointNrMod (eli, elr) == ed[0] &&
			      GetPointNrMod (eli, elr+2) == ed[1])
			    {
			      /*
			      (*testout) << "ed is diagonal of rectangle" << endl;
			      (*testout) << "ed = " << ed.I(1) << "-" << ed.I(2) << endl;
			      (*testout) << "ind = " << ind << endl;
			      */
			      ind = 0;
			    }

			}
		    }
		}

	      if (ind)
		{
		  /*
		  (*testout) << "new edge from face " << k 
			     << " = (" << freesetfaces.Get(k).i1 
			     << ", " << freesetfaces.Get(k).i2 
			     << ", " << freesetfaces.Get(k).i3
			     << "), point " << ind << endl;
			     */
		  freesetedges.Append(twoint(k,ind));
		}
	    }	
	}
    }
    
}





void Meshing3 :: LoadRules (const char * filename, const char ** prules)
{
  char buf[256];
  istream * ist;
  char *tr1 = NULL;

  if (filename)
    {
      PrintMessage (3, "rule-filename = ", filename);
      ist = new ifstream (filename);
    }
  else 
    {
      /* connect tetrules to one string */
      PrintMessage (3, "Use internal rules");
      if (!prules) prules = tetrules;

      const char ** hcp = prules; 
      size_t len = 0;
      while (*hcp)
	{
	  len += strlen (*hcp);
	  hcp++;
	}
      tr1 = new char[len+1];
      tr1[0] = 0;
      hcp = prules; //  tetrules;


      char * tt1 = tr1;
      while (*hcp)
	{
	  strcat (tt1, *hcp);
	  tt1 += strlen (*hcp);	  
	  hcp++;
	}


#ifdef WIN32
      // VC++ 2005 workaround
      for(size_t i=0; i<len; i++)
	if(tr1[i] == ',')
	  tr1[i] = ':';
#endif

      ist = new istringstream (tr1);
    }
  
  if (!ist->good())
    {
      cerr << "Rule description file " << filename << " not found" << endl;
      delete ist;
      exit (1);
    }
    
  while (!ist->eof())
    {
      buf[0] = 0;
      (*ist) >> buf;
	
      if (strcmp (buf, "rule") == 0)
	{
	  // vnetrule * rule = new vnetrule;
          auto rule = make_unique<vnetrule>();
	  rule -> LoadRule(*ist);
	  if (!rule->TestOk())
	    {
	      PrintSysError ("Parser3d: Rule ", rules.Size(), " not ok");
	      exit (1);
	    }
	  rules.Append (std::move(rule));
	}
      else if (strcmp (buf, "tolfak") == 0)
	{
	  (*ist) >> tolfak;
	}
    }
  delete ist;
  delete [] tr1;
}
}
