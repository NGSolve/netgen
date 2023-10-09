//
//
// TECPLOT file by Jawor Georgiew
//
#include <mystdlib.h>

#include <myadt.hpp>
#include <linalg.hpp>
#include <csg.hpp>
#include <meshing.hpp>

#include "writeuser.hpp"


namespace netgen
{

void WriteTecPlotFormat (const Mesh & mesh,
			 const filesystem::path & filename)
{
  auto geom = dynamic_pointer_cast<CSGeometry>(mesh.GetGeometry());
  if(geom == nullptr)
    throw Exception("TecPlot format requires a CSGeometry");

  INDEX i;
  int j, k, e, z;
  Vec<3> n;
  
  INDEX np = mesh.GetNP();
  INDEX ne = mesh.GetNE();
  INDEX nse = mesh.GetNSE();
  
  NgArray<int> sn(np);
  ofstream outfile(filename);
  
  outfile << "TITLE=\" " << filename.string() << "\"" << endl;

  // fill hashtable

  INDEX_3_HASHTABLE<int> face2volelement(ne);

  for (i = 1; i <= ne; i++)
    {
      const Element & el = mesh.VolumeElement(i);
      INDEX_3 i3;
      int l;
      for (j = 1; j <= 4; j++)   // loop over faces of tet
	{
	  l = 0;
	  for (k = 1; k <= 4; k++)
	    if (k != j)
	      {
		l++;
		i3.I(l) = el.PNum(k);
	      }
	  i3.Sort();
	  face2volelement.Set (i3, i);
	}
    }
      
      
  for (j = 1; j <= geom->GetNSurf(); j++)       /* Flaeche Nummer j */
    {
      for (i = 1; i <= np; i++)
	sn.Elem(i) = 0;

      e = 0;
       
      for (i = 1; i <= nse; i++)
	{
	  const Element2d & el = mesh.SurfaceElement(i);
	  if (j ==  mesh.GetFaceDescriptor (el.GetIndex ()).SurfNr())
	    {
	      for (k = 1; k <= 3; k++)
		sn.Elem(el.PNum(k)) = 1;
	      e++;                     /* e= Anzahl der neuen Elemente */
	    }
	}

      z = 0;
      for (i = 1; i <= np; i++)
	if (sn.Elem(i) == 1)
	  sn.Elem(i) = ++z;

      outfile << "ZONE T=\" Surface " << j << " \", N=" << z
	      << ", E=" << e << ", ET=TRIANGLE, F=FEPOINT" << endl;

      for (i = 1; i <= np; i++)
	if (sn.Elem(i) != 0)
	  {
	    n = geom->GetSurface(j) -> GetNormalVector ( mesh.Point(i) );
		
	    outfile << mesh.Point(i)(0) << " " /* Knoten Koordinaten */
		    << mesh.Point(i)(1) << " "
		    << mesh.Point(i)(2) << " "
		    << n(0) << " "
		    << n(1) << " "
		    << n(2) << " "
		    << i     << endl;
	  }
	  

      for (i = 1; i <= nse; i++)
	{
	  const Element2d & el = mesh.SurfaceElement(i);
	  if (j ==  mesh.GetFaceDescriptor(el.GetIndex ()).SurfNr())
	    /* FlaechenKnoten (3) */
	    outfile << sn.Get(el.PNum(1)) << " " 
		    << sn.Get(el.PNum(2)) << " "
		    << sn.Get(el.PNum(3)) << endl;
	      
	  /// Hier soll noch die Ausgabe der Nummer des angrenzenden
	      /// Vol.elements erfolgen !

	      for (k = 1; k <= nse; k++)
		{
		  const Element2d & sel = mesh.SurfaceElement(k);
		  INDEX_3 i3;
		  for (j = 1; j <= 3; j++)
		    i3.I(j) = sel.PNum(j);
		  i3.Sort();
		  
		  //int elind = face2volelement.Get(i3);
		}
	}
    }
}

static RegisterUserFormat reg_tecplot ("TecPlot Format", {".mesh"}, nullopt, WriteTecPlotFormat);

}
