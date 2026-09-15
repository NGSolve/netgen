#include <mystdlib.h>
#include "meshing.hpp"


namespace netgen
{

// A special function for Hermann Landes, Erlangen


void CutOffAndCombine (Mesh & mesh, const Mesh & othermesh)
{
  int i, j;
  int nse = othermesh.GetNSE();
  int onp = othermesh.GetNP();


  PrintMessage (1, "other mesh has ",
		othermesh.GetNP(), " points, ",
		othermesh.GetNSE(), " surface elements.");

  Array<Box3d> otherbounds(nse);  
  Box3d otherbox;

  double maxh = 0;
  for (SurfaceElementIndex i : T_Range<SurfaceElementIndex>(nse))
    {
      const Element2d & sel = othermesh[i];
      sel.GetBox(othermesh.Points(), otherbounds[i.Nr1()-1]);

      double loch = othermesh.GetH (othermesh.Point (sel.PNum(1)));
      otherbounds[i.Nr1()-1].Increase(loch);
      if (loch > maxh) maxh = loch;
    }

  otherbox.SetPoint (othermesh[IndexBASE<PointIndex>()]);
  for (PointIndex pi : othermesh.Points().Range())
    otherbox.AddPoint (othermesh[pi]);
  otherbox.Increase (maxh);

  for (ElementIndex i : mesh.VolumeElements().Range())
    {
      Box3d box;
      int remove = 0;

      const Element & el = mesh[i];
      el.GetBox(mesh.Points(), box);

      if (i.Nr1() % 10000 == 0)
	cout << "+" << flush;

      if (box.Intersect(otherbox))
	{
	  for (j = 1; j <= nse && !remove; j++)
	    if (box.Intersect(otherbounds[j-1]))
	      remove = 1;
	}

      if (remove)
	mesh[i].Delete();
    }
  cout << endl;

  TBitArray<PointIndex> connected(mesh.GetNP());
  connected.Clear();
  for (auto & el : mesh.SurfaceElements())
    {
      for (j = 1; j <= 3; j++)
	connected.SetBit(el.PNum(j));
    }
  
  bool changed;
  do
    {
      changed = 0;
      for (auto & el : mesh.VolumeElements())
	{
	  int has = 0, hasnot = 0;
	  if (el[0].IsValid())
	    {
	      for (j = 0; j < 4; j++)
		{
		  if (connected.Test(el[j]))
		    has = 1;
		  else
		    hasnot = 1;
		}
	      if (has && hasnot)
		{
		  changed = 1;
		  for (j = 0; j < 4; j++)
		    connected.SetBit (el[j]);
		}
	    }
	}
      cout << "." << flush;
    }
  while (changed);
  cout << endl;

  for (auto & el : mesh.VolumeElements())
    {
      int hasnot = 0;
      if (el[0].IsValid())
	{
	  for (j = 0; j < 4; j++)
	    {
	      if (!connected.Test(el[j]))
		hasnot = 1;
	    }
	  if (hasnot)
	    el.Delete();
	}
    }

  mesh.Compress();
  
  mesh.FindOpenElements();
  TBitArray<PointIndex> locked(mesh.GetNP());
  locked.Set();
  for (i = 1; i <= mesh.GetNOpenElements(); i++)
    for (j = 1; j <= 3; j++)
      locked.Clear (mesh.OpenElement(i).PNum(j));

  // for (PointIndex i (1); i <= locked.Size(); i++)
  for (PointIndex i : locked.Range())
    if (locked.Test(i))
      {
	mesh.AddLockedPoint (i);
      }



  
  Array<PointIndex, PointIndex> pmat(onp);
  for (PointIndex pi : othermesh.Points().Range())
    pmat[pi] = mesh.AddPoint (othermesh[pi]);

  int fnum = 
    mesh.AddFaceDescriptor (FaceDescriptor(0,0,1,0));

  for (auto & sel : othermesh.SurfaceElements())
    {
      Element2d tri = sel;
      for (j = 1; j <= 3; j++)
	tri.PNum(j) = pmat[tri.PNum(j)];
      tri.SetIndex(fnum);
      mesh.AddSurfaceElement (tri);
    }

  for (PointIndex pi : pmat.Range())
    mesh.AddLockedPoint (pmat[pi]);

  mesh.CalcSurfacesOfNode();
  mesh.CalcLocalH(0.3);
}




void HelmholtzMesh (Mesh & mesh)
{
  int i;
  double ri, ra, rinf;

  cout << "ri = ";
  cin >> ri;
  cout << "ra = ";
  cin >> ra;
  cout << "rinf = ";
  cin >> rinf;

  double det = ri * ra * rinf - ri * ri * rinf;
  double a = (ri - rinf) / det;
  double b = (ri*ri - ra * rinf) / det;
  for (PointIndex pi : mesh.Points().Range())
    {
      Point<3> & p = mesh[pi];
      double rold = sqrt (sqr(p(0)) + sqr(p(1)) + sqr(p(2)));
      if (rold < ri) continue;

      double rnew = 1 / (a * rold - b);
      double fac = rnew / rold;
      p(0) *= fac;
      p(1) *= fac;
      p(2) *= fac;
    }
}
}
