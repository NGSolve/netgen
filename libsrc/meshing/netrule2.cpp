#include <mystdlib.h>
#include "meshing.hpp"

namespace netgen
{

netrule :: netrule ()
{
  name = new char[1];
  name[0] = char(0);
  quality = 0;
}

netrule ::  ~netrule()
{
  delete [] name;
  for(int i = 0; i < oldutofreearea_i.Size(); i++)
    delete oldutofreearea_i[i];
  for(int i = 0; i < freezone_i.Size(); i++)
    delete freezone_i[i];
}



void netrule :: SetFreeZoneTransformation (const Vector & devp, int tolclass)
{
  double lam1 = 1.0/tolclass;
  double lam2 = 1.-lam1;

  double mem1[100], mem2[100], mem3[100];

  int vs = oldutofreearea.Height();
  FlatVector devfree(vs, mem1);

  int fzs = freezone.Size();
  transfreezone.SetSize (fzs);

  if (tolclass <= oldutofreearea_i.Size())
    {
      oldutofreearea_i[tolclass-1] -> Mult (devp, devfree);

      auto& fzi = *freezone_i[tolclass-1];
      for (int i = 0; i < fzs; i++)
	{
	  transfreezone[i][0] = fzi[i][0] + devfree[2*i];
	  transfreezone[i][1] = fzi[i][1] + devfree[2*i+1];
	}
    }
  else
    {
      FlatVector devfree1(vs, mem2);
      FlatVector devfree2(vs, mem3);

      oldutofreearea.Mult (devp, devfree1);
      oldutofreearealimit.Mult (devp, devfree2);
      devfree.Set2 (lam1, devfree1, lam2, devfree2);

      for (int i = 0; i < fzs; i++)
	{
	  transfreezone[i][0] = lam1 * freezone[i][0] + lam2 * freezonelimit[i][0] + devfree[2*i];
	  transfreezone[i][1] = lam1 * freezone[i][1] + lam2 * freezonelimit[i][1] + devfree[2*i+1];
	}
    }


  if (fzs > 0)
    {
      fzmaxx = fzminx = transfreezone[0][0];
      fzmaxy = fzminy = transfreezone[0][1];
    }

  for (int i = 1; i < fzs; i++)
    {
      if (transfreezone[i][0] > fzmaxx) fzmaxx = transfreezone[i][0];
      if (transfreezone[i][0] < fzminx) fzminx = transfreezone[i][0];
      if (transfreezone[i][1] > fzmaxy) fzmaxy = transfreezone[i][1];
      if (transfreezone[i][1] < fzminy) fzminy = transfreezone[i][1];
    }

  for (int i = 0; i < fzs; i++)
    {
      const auto& p1 = transfreezone[i];
      const auto& p2 = transfreezone[(i+1) % fzs];

      Vec<2> vn = { p2[1] - p1[1], p1[0] - p2[0] };

      double len2 = vn.Length2();

      if (len2 < 1e-10)
	{
	  freesetinequ(i, 0) = 0;
	  freesetinequ(i, 1) = 0;
	  freesetinequ(i, 2) = -1;
	}
      else
	{
	  vn /= sqrt (len2);    // scaling necessary ?

	  freesetinequ(i,0) = vn[0]; 
	  freesetinequ(i,1) = vn[1]; 
	  freesetinequ(i,2) = -(p1[0] * vn[0] + p1[1] * vn[1]);
	}
    }
}


/*
int netrule :: IsInFreeZone2 (const Point2d & p) const
{
  for (int i = 0; i < transfreezone.Size(); i++)
    {
      if (freesetinequ(i, 0) * p.X() + 
	  freesetinequ(i, 1) * p[1] +
	  freesetinequ(i, 2) > 0) return 0;
    }
  return 1;
}
*/

int netrule :: IsLineInFreeZone2 (const Point<2> & p1, const Point<2> & p2) const
{
  if ( (p1[0] > fzmaxx && p2[0] > fzmaxx) ||
       (p1[0] < fzminx && p2[0] < fzminx) ||
       (p1[1] > fzmaxy && p2[1] > fzmaxy) ||
       (p1[1] < fzminy && p2[1] < fzminy) ) return 0;

  for (int i = 1; i <= transfreezone.Size(); i++)
    {
      if (freesetinequ.Get(i, 1) * p1[0] + freesetinequ.Get(i, 2) * p1[1] +
	  freesetinequ.Get(i, 3) > -1e-8 &&    // -1e-6
	  freesetinequ.Get(i, 1) * p2[0] + freesetinequ.Get(i, 2) * p2[1] +
	  freesetinequ.Get(i, 3) > -1e-8       // -1e-6
	  ) return 0;
    }

  double nx =  (p2[1] - p1[1]);
  double ny = -(p2[0] - p1[0]);
  double nl = sqrt (nx * nx + ny * ny);
  if (nl > 1e-8)
    {
      nx /= nl;
      ny /= nl;
      double c = - (p1[0] * nx + p1[1] * ny);

      bool allleft = true;
      bool allright = true;

      for (int i = 1; i <= transfreezone.Size(); i++)
	{
	  bool left  = transfreezone.Get(i)[0] * nx + transfreezone.Get(i)[1] * ny + c <  1e-7;
          bool right = transfreezone.Get(i)[0] * nx + transfreezone.Get(i)[1] * ny + c > -1e-7;
	  if (!left) allleft = false;
	  if (!right) allright = false;
	}
      if (allleft || allright) return false;
    }

  return true;
}

int netrule :: ConvexFreeZone () const
{
  int n = transfreezone.Size();
  for (int i = 1; i <= n; i++)
    {
      const bool counterclockwise = CCW (transfreezone.Get(i), 
					 transfreezone.Get(i % n + 1),
					 transfreezone.Get( (i+1) % n + 1 ),
					 1e-7);
      //(*testout) << "ccw " << counterclockwise << endl << " p1 " << transfreezone.Get(i) << " p2 " << transfreezone.Get(i % n + 1)
      //		 << " p3 " << transfreezone.Get( (i+1) % n + 1 ) << endl;
      if (!counterclockwise )
	return 0;
    }
  return 1;
}


/*
float netrule :: CalcPointDist (int pi, const Point2d & p) const
{
  float dx = p.X() - points.Get(pi).X();
  float dy = p.Y() - points.Get(pi).Y();
  const threefloat * tf = &tolerances.Get(pi);

  return tf->f1 * dx * dx + tf->f2 * dx * dy + tf->f3 * dy * dy;
}
*/

float netrule :: CalcLineError (int li, const Vec<2> & v) const
{
  float dx = v[0] - linevecs.Get(li)[0];
  float dy = v[1] - linevecs.Get(li)[1];

  const threefloat * ltf = &linetolerances.Get(li);
  return ltf->f1 * dx * dx + ltf->f2 * dx * dy + ltf->f3 * dy * dy;
}
} // namespace netgen
