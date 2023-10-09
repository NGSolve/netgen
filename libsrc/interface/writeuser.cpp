//
//  Write user dependent output file
//

#include <mystdlib.h>

#include <myadt.hpp>
#include <linalg.hpp>
#include <csg.hpp>
#include <geometry2d.hpp>
#include <meshing.hpp>

#include "writeuser.hpp"
#include "../general/gzstream.h"


namespace netgen
{
  extern MeshingParameters mparam;

  Array<UserFormatRegister::UserFormatEntry> UserFormatRegister::entries;
  std::map<string, int> UserFormatRegister::format_to_entry_index;

  void RegisterUserFormats (NgArray<const char*> & names,
			    NgArray<const char*> & extensions)
			    
{
  for (const auto & entry : UserFormatRegister::entries)
    {
      names.Append (entry.format.c_str());
      extensions.Append (entry.extensions[0].c_str());
    }
}
  
bool WriteUserFormat (const string & format,
		      const Mesh & mesh,
		      const filesystem::path & filename)
{
  if(!UserFormatRegister::HaveFormat(format))
    return true;

  const auto & entry = UserFormatRegister::Get(format);
  if(!entry.write)
    return true;

  (*entry.write)(mesh, filename);
  return false;
}



/*
 *  Neutral mesh format
 *  points, elements, surface elements
 */

void WriteNeutralFormat (const Mesh & mesh,
			 const filesystem::path & filename)
{
  cout << "write neutral, new" << endl;
  int np = mesh.GetNP();
  int ne = mesh.GetNE();
  int nse = mesh.GetNSE();
  int nseg = mesh.GetNSeg();
  int i, j;

  int inverttets = mparam.inverttets;
  int invertsurf = mparam.inverttrigs;

  ofstream outfile (filename);

  outfile.precision(6);
  outfile.setf (ios::fixed, ios::floatfield);
  outfile.setf (ios::showpoint);

  outfile << np << "\n";

  for (i = 1; i <= np; i++)
    {
      const Point3d & p = mesh.Point(i);

      outfile.width(10);
      outfile << p.X() << " ";
      outfile.width(9);
      outfile << p.Y() << " ";
      if (mesh.GetDimension() == 3)
	{
	  outfile.width(9);
	  outfile << p.Z();
	  }
      outfile << "\n";
    }

  if (mesh.GetDimension() == 3)
    {
      outfile << ne << "\n";
      for (i = 1; i <= ne; i++)
	{
	  Element el = mesh.VolumeElement(i);
	  if (inverttets)
	    el.Invert();
	  outfile.width(4);
	  outfile << el.GetIndex() << "  ";
	  for (j = 1; j <= el.GetNP(); j++)
	    {
	      outfile << " ";
	      outfile.width(8);
	      outfile << el.PNum(j);
	    }
	  outfile << "\n";
	}
    }

  outfile << nse << "\n";
  for (i = 1; i <= nse; i++)
    {
      Element2d el = mesh.SurfaceElement(i);
      if (invertsurf)
	el.Invert();
      outfile.width(4);
      outfile << mesh.GetFaceDescriptor (el.GetIndex()).BCProperty() << "    ";
      for (j = 1; j <= el.GetNP(); j++)
	{
	  outfile << " ";
	  outfile.width(8);
	  outfile << el.PNum(j);
	}
      outfile << "\n";
    }


  if (mesh.GetDimension() == 2)
    {
      outfile << nseg << "\n";
      for (int i = 1; i <= nseg; i++)
	{
	  const Segment & seg = mesh.LineSegment(i);
	  outfile.width(4);
	  outfile << seg.si << "    ";

          for (int j = 0; j < seg.GetNP(); j++)
            {
              outfile << " ";
              outfile.width(8);
              outfile << seg[j];
            }
          /*
	  outfile << " ";
	  outfile.width(8);
	  outfile << seg[0];
	  outfile << " ";
	  outfile.width(8);
	  outfile << seg[1];
          if (seg[2] != -1)
            {
              outfile.width(8);
              outfile << seg[2];
            }
          */
	  outfile << "\n";
	}
    }
}









void WriteSurfaceFormat (const Mesh & mesh,
			 const filesystem::path & filename)
{
  // surface mesh
  int i, j;

  cout << "Write Surface Mesh" << endl;

  ofstream outfile (filename);

  outfile << "surfacemesh" << endl;

  outfile << mesh.GetNP() << endl;
  for (i = 1; i <= mesh.GetNP(); i++)
    {
      for (j = 0; j < 3; j++)
	{
	  outfile.width(10);
	  outfile << mesh.Point(i)(j) << " ";
	}
      outfile << endl;
    }
  outfile << mesh.GetNSE() << endl;
  for (i = 1; i <= mesh.GetNSE(); i++)
    {
      for (j = 1; j <= 3; j++)
	{
	  outfile.width(8);
	  outfile << mesh.SurfaceElement(i).PNum(j);
	}
      outfile << endl;
    }
}





/*
 *  save surface mesh as STL file
 */

void WriteSTLFormat (const Mesh & mesh,
		     const filesystem::path & filename)
{
  cout << "\nWrite STL Surface Mesh" << endl;

  auto ext = filename.extension();
  ostream *outfile;

  if(ext == ".gz")
	  outfile = new ogzstream(filename);
  else
	  outfile = new ofstream(filename);

  int i;

  outfile->precision(10);

  *outfile << "solid" << endl;

  for (i = 1; i <= mesh.GetNSE(); i++)
    {
      *outfile << "facet normal ";
      const Point3d& p1 = mesh.Point(mesh.SurfaceElement(i).PNum(1));
      const Point3d& p2 = mesh.Point(mesh.SurfaceElement(i).PNum(2));
      const Point3d& p3 = mesh.Point(mesh.SurfaceElement(i).PNum(3));

      Vec3d normal = Cross(p2-p1,p3-p1);
      if (normal.Length() != 0)
	{
	  normal /= (normal.Length());
	}

      *outfile << normal.X() << " " << normal.Y() << " " << normal.Z() << "\n";
      *outfile << "outer loop\n";

      *outfile << "vertex " << p1.X() << " " << p1.Y() << " " << p1.Z() << "\n";
      *outfile << "vertex " << p2.X() << " " << p2.Y() << " " << p2.Z() << "\n";
      *outfile << "vertex " << p3.X() << " " << p3.Y() << " " << p3.Z() << "\n";

      *outfile << "endloop\n";
      *outfile << "endfacet\n";
    }
  *outfile << "endsolid" << endl;
}





/*
 *  Philippose - 16 August 2010
 *  Save surface mesh as STL file
 *  with a separate solid definition
 *  for each face
 *  - This helps in splitting up the
 *    STL into named boundary faces
 *    when using a third-party mesher
 */
void WriteSTLExtFormat (const Mesh & mesh,
		     const filesystem::path & filename)
{
  cout << "\nWrite STL Surface Mesh (with separated boundary faces)" << endl;

  auto ext = filename.extension();
  ostream *outfile;

  if(ext == ".gz")
	  outfile = new ogzstream(filename);
  else
	  outfile = new ofstream(filename);

  outfile->precision(10);

  int numBCs = 0;

  NgArray<int> faceBCs;
  TABLE<int> faceBCMapping;

  faceBCs.SetSize(mesh.GetNFD());
  faceBCMapping.SetSize(mesh.GetNFD());

  faceBCs = -1;

  // Collect the BC numbers used in the mesh
  for(int faceNr = 1; faceNr <= mesh.GetNFD(); faceNr++)
  {
	  int bcNum = mesh.GetFaceDescriptor(faceNr).BCProperty();

	  if(faceBCs.Pos(bcNum) < 0)
	  {
        numBCs++;
		  faceBCs.Set(numBCs,bcNum);
        faceBCMapping.Add1(numBCs,faceNr);        
	  }
     else
     {
        faceBCMapping.Add1(faceBCs.Pos(bcNum)+1,faceNr);
     }
  }

  faceBCs.SetSize(numBCs);
  faceBCMapping.ChangeSize(numBCs);

  // Now actually write the data to file
  for(int bcInd = 1; bcInd <= faceBCs.Size(); bcInd++)
  {
      *outfile << "solid Boundary_" << faceBCs.Elem(bcInd) << "\n";

      for(int faceNr = 1;faceNr <= faceBCMapping.EntrySize(bcInd); faceNr++)
      {
        Array<SurfaceElementIndex> faceSei;
          mesh.GetSurfaceElementsOfFace(faceBCMapping.Get(bcInd,faceNr),faceSei);

          for (int i = 0; i < faceSei.Size(); i++)
          {
        	  *outfile << "facet normal ";
        	  const Point3d& p1 = mesh.Point(mesh[faceSei[i]].PNum(1));
        	  const Point3d& p2 = mesh.Point(mesh[faceSei[i]].PNum(2));
        	  const Point3d& p3 = mesh.Point(mesh[faceSei[i]].PNum(3));

        	  Vec3d normal = Cross(p2-p1,p3-p1);
        	  if (normal.Length() != 0)
        	  {
        		  normal /= (normal.Length());
        	  }

        	  *outfile << normal.X() << " " << normal.Y() << " " << normal.Z() << "\n";
        	  *outfile << "outer loop\n";

        	  *outfile << "vertex " << p1.X() << " " << p1.Y() << " " << p1.Z() << "\n";
        	  *outfile << "vertex " << p2.X() << " " << p2.Y() << " " << p2.Z() << "\n";
        	  *outfile << "vertex " << p3.X() << " " << p3.Y() << " " << p3.Z() << "\n";

        	  *outfile << "endloop\n";
        	  *outfile << "endfacet\n";
          }
      }
      *outfile << "endsolid Boundary_" << faceBCs.Elem(bcInd) << "\n";
  }
}




/*
 *
 *  write surface mesh as VRML file
 *
 */

void WriteVRMLFormat (const Mesh & mesh,
		      bool faces,
		      const filesystem::path & filename)
{

  if (faces)

    {
      // Output in VRML, IndexedFaceSet is used
      // Bartosz Sawicki <sawickib@ee.pw.edu.pl>

      int np = mesh.GetNP();
      int nse = mesh.GetNSE();
      int i, j;

      ofstream outfile (filename);

      outfile.precision(6);
      outfile.setf (ios::fixed, ios::floatfield);
      outfile.setf (ios::showpoint);

      outfile << "#VRML V2.0 utf8 \n"
	         "Background {\n"
		 "    skyColor [1 1 1]\n"
     		 "    groundColor [1 1 1]\n"
		 "}\n"
		 "Group{ children [\n"
		 "Shape{ \n"
		 "appearance Appearance { material Material { }} \n"
                 "geometry IndexedFaceSet { \n"
                 "coord Coordinate { point [ \n";


      for (i = 1; i <= np; i++)
        {
          const Point3d & p = mesh.Point(i);
          outfile.width(10);
          outfile << p.X() << " ";
          outfile << p.Y() << " ";
          outfile << p.Z() << " \n";
	}

      outfile << "  ] } \n"
                 "coordIndex [ \n";

      for (i = 1; i <= nse; i++)
	{
	  const Element2d & el = mesh.SurfaceElement(i);

	  for (j = 1; j <= 3; j++)
	    {
	      outfile.width(8);
	      outfile << el.PNum(j)-1;
	    }
	  outfile << " -1 \n";
	}

      outfile << "  ] \n";

      //define number and RGB definitions of colors
      outfile << "color Color { color [1 0 0, 0 1 0, 0 0 1, 1 1 0]} \n"
                 "colorIndex [\n";

      for (i = 1; i <= nse; i++)
	{
	  outfile << mesh.GetFaceDescriptor(mesh.SurfaceElement(i).GetIndex ()).BCProperty();
          outfile << endl;
	}

      outfile << " ] \n"
                 "colorPerVertex FALSE \n"
                 "creaseAngle 0 \n"
		 "solid FALSE \n"
                 "ccw FALSE \n"
		 "convex TRUE \n"
                 "} } # end of Shape\n"
		 "] }\n";

    } /* end of VRMLFACES */


  else

    {
        // Output in VRML, IndexedLineSet is used
	// Bartosz Sawicki <sawickib@ee.pw.edu.pl>

      int np = mesh.GetNP();
      int nse = mesh.GetNSE();
      int i, j;

      ofstream outfile (filename);

      outfile.precision(6);
      outfile.setf (ios::fixed, ios::floatfield);
      outfile.setf (ios::showpoint);

      outfile << "#VRML V2.0 utf8 \n"
	         "Background {\n"
		 "    skyColor [1 1 1]\n"
     		 "    groundColor [1 1 1]\n"
		 "}\n"
		 "Group{ children [\n"
	         "Shape{ \n"
		 "appearance Appearance { material Material { }} \n"
                 "geometry IndexedLineSet { \n"
                 "coord Coordinate { point [ \n";


      for (i = 1; i <= np; i++)
        {
          const Point3d & p = mesh.Point(i);
          outfile.width(10);
          outfile << p.X() << " ";
          outfile << p.Y() << " ";
          outfile << p.Z() << " \n";
	}

      outfile << "  ] } \n"
                 "coordIndex [ \n";

      for (i = 1; i <= nse; i++)
	{
	  const Element2d & el = mesh.SurfaceElement(i);

	  for (j = 1; j <= 3; j++)
	    {
	      outfile.width(8);
	      outfile << el.PNum(j)-1;
	    }
	  outfile.width(8);
	  outfile << el.PNum(1)-1;
	  outfile << " -1 \n";
	}

      outfile << "  ] \n";

/* Uncomment if you want color mesh
      outfile << "color Color { color [1 1 1, 0 1 0, 0 0 1, 1 1 0]} \n"
                 "colorIndex [\n";

      for (i = 1; i <= nse; i++)
	{
	  outfile << mesh.GetFaceDescriptor(mesh.SurfaceElement(i).GetIndex ()).BCProperty();
          outfile << endl;
	}

      outfile << " ] \n"
*/
      outfile << "colorPerVertex FALSE \n"
                 "} } #end of Shape\n"
		 "] } \n";

    }

}

void WriteVRMLFormatLineset (const Mesh & mesh, const filesystem::path & filename)
{
  WriteVRMLFormat(mesh, false, filename);
}

void WriteVRMLFormatFaceset (const Mesh & mesh, const filesystem::path & filename)
{
  WriteVRMLFormat(mesh, true, filename);
}






/*
 * FEPP .. a finite element package developed at University Linz, Austria
 */
void WriteFEPPFormat (const Mesh & mesh,
		      const filesystem::path & filename)
{

  ofstream outfile (filename);

  if (mesh.GetDimension() == 3)

    {

      // output for FEPP

      int np = mesh.GetNP();
      int ne = mesh.GetNE();
      int nse = mesh.GetNSE();
      // int ns = mesh.GetNFD();
      int i, j;

      outfile.precision(5);
      outfile.setf (ios::fixed, ios::floatfield);
      outfile.setf (ios::showpoint);

      outfile << "volumemesh4" << endl;
      outfile << nse << endl;
      for (i = 1; i <= nse; i++)
	{
	  const Element2d & el = mesh.SurfaceElement(i);

	  //	  int facenr = mesh.facedecoding.Get(el.GetIndex()).surfnr;
	  outfile.width(4);
	  outfile << el.GetIndex() << " ";
	  outfile.width(4);
	  //	  outfile << mesh.GetFaceDescriptor(el.GetIndex()).BCProperty() << " ";
	  outfile << mesh.GetFaceDescriptor(el.GetIndex()).BCProperty() << " ";
	  outfile.width(4);
	  outfile << el.GetNP() << "    ";
	  for (j = 1; j <= el.GetNP(); j++)
	    {
	      outfile.width(8);
	      outfile << el.PNum(j);
	    }
	  outfile << "\n";
	}


      outfile << ne << "\n";
      for (i = 1; i <= ne; i++)
	{
	  const Element & el = mesh.VolumeElement(i);
	  outfile.width(4);
	  outfile << el.GetIndex() << " ";
	  outfile.width(4);
	  outfile << el.GetNP() << " ";
	  for (j = 1; j <= el.GetNP(); j++)
	    {
	      outfile.width(8);
	      outfile << el.PNum(j);
	    }
	  outfile << "\n";
	}

      outfile << np << "\n";
      for (i = 1; i <= np; i++)
	{
	  const Point3d & p = mesh.Point(i);

	  outfile.width(10);
	  outfile << p.X() << " ";
	  outfile.width(9);
	  outfile << p.Y() << " ";
	  outfile.width(9);
	  outfile << p.Z() << "\n";
	}

      /*
      if (typ == WRITE_FEPPML)
	{
	  int nbn =  mesh.mlbetweennodes.Size();
	  outfile << nbn << "\n";
	  for (i = 1; i <= nbn; i++)
	    outfile << mesh.mlbetweennodes.Get(i).I1() << " "
		    << mesh.mlbetweennodes.Get(i).I2() << "\n";


	  //	  int ncon = mesh.connectedtonode.Size();
	  //	  outfile << ncon << "\n";
	  //	  for (i = 1; i <= ncon; i++)
	  //	    outfile << i << " " << mesh.connectedtonode.Get(i) << endl;
	}
      */

      /*
      // write CSG surfaces
      if (&geom && geom.GetNSurf() >= ns)
	{
	  outfile << ns << endl;
	  for (i = 1; i <= ns; i++)
	    geom.GetSurface(mesh.GetFaceDescriptor(i).SurfNr())->Print(outfile);
	}
      else
      */
	outfile << "0" << endl;
    }


  else

    { // 2D fepp format

      ;
      /*
      extern SplineGeometry2d * geometry2d;
      if (geometry2d)
	Save2DMesh (mesh, &geometry2d->GetSplines(), outfile);
      else
	Save2DMesh (mesh, 0, outfile);
      */
    }
}






/*
 *  Edge element mesh format
 *  points, elements, edges
 */

void WriteEdgeElementFormat (const Mesh & mesh,
			     const filesystem::path & filename)
{
  cout << "write edge element format" << endl;

  const MeshTopology * top = &mesh.GetTopology();
  int npoints = mesh.GetNP();
  int nelements = mesh.GetNE();
  int nsurfelem = mesh.GetNSE();
  int nedges = top->GetNEdges();
  int i, j;

  int inverttets = mparam.inverttets;
  int invertsurf = mparam.inverttrigs;
  NgArray<int> edges;

  ofstream outfile (filename);

  outfile.precision(6);
  outfile.setf (ios::fixed, ios::floatfield);
  outfile.setf (ios::showpoint);


  // vertices with coordinates
  outfile << npoints << "\n";
  for (i = 1; i <= npoints; i++)
    {
      const Point3d & p = mesh.Point(i);

      outfile.width(10);
      outfile << p.X() << " ";
      outfile.width(9);
      outfile << p.Y() << " ";
      outfile.width(9);
      outfile << p.Z() << "\n";
    }

  // element - edge - list
  outfile << nelements << " " << nedges << "\n";
  for (i = 1; i <= nelements; i++)
    {
      Element el = mesh.VolumeElement(i);
      if (inverttets)
      	el.Invert();
      outfile.width(4);
      outfile << el.GetIndex() << "  ";
      outfile.width(8);
      outfile << el.GetNP();
      for (j = 1; j <= el.GetNP(); j++)
	{
	  outfile << " ";
	  outfile.width(8);
	  outfile << el.PNum(j);
	}

      // top->GetElementEdges(i,edges);
      auto eledges = top->GetEdges(ElementIndex(i-1));
      outfile << endl << "      ";
      outfile.width(8);
      outfile << eledges.Size();
      for (j=1; j <= eledges.Size(); j++)
	{
	  outfile << " ";
	  outfile.width(8);
	  outfile << eledges[j-1]+1;
	}
      outfile << "\n";

      // orientation:
      top->GetElementEdgeOrientations(i,edges);
      outfile << "              ";
      for (j=1; j <= edges.Size(); j++)
	{
	  outfile << " ";
	  outfile.width(8);
	  outfile << edges[j-1];
	}
      outfile << "\n";
    }

  // surface element - edge - list (with boundary conditions)
  outfile << nsurfelem << "\n";
  for (i = 1; i <= nsurfelem; i++)
    {
      Element2d el = mesh.SurfaceElement(i);
      if (invertsurf)
	el.Invert();
      outfile.width(4);
      outfile << mesh.GetFaceDescriptor (el.GetIndex()).BCProperty() << "  ";
      outfile.width(8);
      outfile << el.GetNP();
      for (j = 1; j <= el.GetNP(); j++)
	{
	  outfile << " ";
	  outfile.width(8);
	  outfile << el.PNum(j);
	}

      top->GetSurfaceElementEdges(i,edges);
      outfile << endl << "      ";
      outfile.width(8);
      outfile << edges.Size();
      for (j=1; j <= edges.Size(); j++)
	{
	  outfile << " ";
	  outfile.width(8);
	  outfile << edges[j-1];
	}
      outfile << "\n";
    }


  // int v1, v2;
  // edge - vertex - list
  outfile << nedges << "\n";
  for (i=1; i <= nedges; i++)
    {
      // top->GetEdgeVertices(i,v1,v2);
      auto [v1,v2] = top->GetEdgeVertices(i-1);
      outfile.width(4);
      outfile << v1;
      outfile << " ";
      outfile.width(8);
      outfile << v2 << endl;
    }
}

static RegisterUserFormat reg_neutral("Neutral Format", {".mesh"}, ReadFile, WriteNeutralFormat);
static RegisterUserFormat reg_surface("Surface Mesh Format", {".mesh", ".surf"} ,ReadFile, WriteSurfaceFormat);
static RegisterUserFormat reg_stl ("STL Format", {".stl", ".stlb"}, ReadFile, WriteSTLFormat);
static RegisterUserFormat reg_stlx ("STL Extended Format", {".stl", ".stlb"}, nullopt, WriteSTLExtFormat);
static RegisterUserFormat reg_vrml ("VRML Format", {".*"}, nullopt, WriteVRMLFormatLineset);
static RegisterUserFormat reg_vrml_faces ("VRML Format Faces", {".*"}, nullopt, WriteVRMLFormatFaceset);
static RegisterUserFormat reg_fepp ("Fepp Format", {"*"}, nullopt, WriteFEPPFormat);
static RegisterUserFormat reg_edgeelement ("EdgeElement Format", {"*"}, nullopt, WriteEdgeElementFormat);

}

