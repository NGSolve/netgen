//
//  Write Fluent file
//  Johannes Gerstmayr, University Linz
//

#include <mystdlib.h>

#include <myadt.hpp>
#include <linalg.hpp>
#include <csg.hpp>
#include <meshing.hpp>

#include "writeuser.hpp"

namespace netgen
{


void WriteFluentFormat (const Mesh & mesh,
                        const filesystem::path & filename)

{
  cout << "start writing fluent export" << endl;
      
  int np = mesh.GetNP();
  int ne = mesh.GetNE();
  int nse = mesh.GetNSE();
  int i, j;

  ofstream outfile (filename);
  char str[100];

  outfile.precision(6);
  //outfile.setf (ios::fixed, ios::floatfield);
  //outfile.setf (ios::showpoint);
      
  outfile << "(0 \"Exported file from NETGEN \")" << endl;
  outfile << "(0 \"Dimension:\")" << endl;
  outfile << "(2 3)" << endl << endl;

  outfile << "(0 \"Nodes:\")" << endl;

  //number of nodes:
  snprintf (str, size(str), "(10 (0 1 %x 1))",np); //hexadecimal!!!
  outfile << str << endl;

  //nodes of zone 1:
  snprintf (str, size(str), "(10 (7 1 %x 1)(",np); //hexadecimal!!!
  outfile << str << endl;
  for (PointIndex pi : mesh.Points().Range())
    {
      const Point<3> & p = mesh[pi];

      //outfile.width(10);
      outfile << p(0) << " ";
      outfile << p(1) << " ";
      outfile << p(2) << "\n";
    }
  outfile << "))" << endl << endl;

  //write faces with elements

  outfile << "(0 \"Faces:\")" << endl;

  Element2d face, face2;
  int /* i2, */ j2;
  Array<PointIndices<3>> surfaceelp;
  Array<int> surfaceeli;
  Array<ElementIndex> locels;

  //no cells=no tets
  //no faces=2*tets

  int noverbface = 2*ne-nse/2;
      
  snprintf (str, size(str), "(13 (0 1 %x 0))",(noverbface+nse)); //hexadecimal!!!
  outfile << str << endl;
      
  snprintf (str, size(str), "(13 (4 1 %x 2 3)(",noverbface); //hexadecimal!!!
  outfile << str << endl;

  const_cast<Mesh&> (mesh).BuildElementSearchTree(3);

  for (ElementIndex i : T_Range<ElementIndex>(ne))
    {
      if (ne > 2000)
        {
          if (i.Nr1()%2000 == 0)
            {
              cout << (double)i.Nr1()/(double)ne*100. << "%" << endl;
            }
        }

      Element el = mesh[i];
      //if (inverttets)
      //  el.Invert();
          
      //outfile << el.GetIndex() << "    ";
      if (el.GetNP() != 4) {cout << "only tet-meshes supported in write fluent!" << endl;}
          
      //faces:
          
      Box3d box;
      el.GetBox(mesh.Points(), box);
      box.IncreaseRel(1e-6);

      mesh.GetIntersectingVolEls(box.PMin(),box.PMax(),locels);
      // int nel = locels.Size();
      // int locind;

      //cout << "nel=" << nel << endl;

      for (j = 1; j <= el.GetNFaces(); j++)
        {
          el.GetFace(j, face);
          face.Invert();
          int eli2 = 0;
          int stopsig = 0;
              
          for (auto locind : locels)
            {
              Element el2 = mesh[locind];
              //if (inverttets)
              //  el2.Invert();

              for (j2 = 1; j2 <= el2.GetNFaces(); j2++)
                {
                  el2.GetFace(j2, face2);

                  if (face2.HasFace(face)) {eli2 = locind.Nr1(); stopsig = 1; break;}
                }
              if (stopsig) break;
            }
              
          if (eli2==i.Nr1()) cout << "error in WRITE_FLUENT!!!" << endl;
              
          if (eli2 > i.Nr1()) //don't write faces two times!
            {
              //i: left cell, eli: right cell
              outfile << hex << face[1] << " "
                << hex << face[0] << " "
                << hex << face[2] << " "
                << hex << i.Nr1()  << " "
                << hex << eli2 << "\n";
            }
          if (eli2 == 0) 
            {
              surfaceelp.Append(PointIndices<3>(face[1],face[0],face[2]));
              surfaceeli.Append(i.Nr1());
            }
        }
    }
  outfile << "))" << endl;
      
  snprintf (str, size(str), "(13 (2 %x %x 3 3)(",(noverbface+1),noverbface+nse); //hexadecimal!!!
  outfile << str << endl;

  for (i = 1; i <= surfaceelp.Size(); i++)
    {
      outfile << hex << surfaceelp[i-1][0].Nr1() << " "
              << hex << surfaceelp[i-1][1].Nr1() << " "
              << hex << surfaceelp[i-1][2].Nr1() << " "
              << hex << surfaceeli[i-1] << " " << 0 << "\n";
    }

  outfile << "))" << endl << endl;

  outfile << "(0 \"Cells:\")" << endl;
      
  snprintf (str, size(str), "(12 (0 1 %x 0))",ne); //hexadecimal!!!
  outfile << str << endl;

  snprintf (str, size(str), "(12 (1 1 %x 1 2))",ne); //hexadecimal!!!
  outfile << str << endl << endl;




  outfile << "(0 \"Zones:\")\n"
          << "(45 (1 fluid fluid)())\n"
    //      << "(45 (2 velocity-inlet velocity_inlet.1)())\n"
    //      << "(45 (3 pressure-outlet pressure_outlet.2)())\n"
          << "(45 (2 wall wall)())\n"
          << "(45 (4 interior default-interior)())\n" << endl;

  cout << "done" << endl;
}

static RegisterUserFormat reg_fluent ("Fluent Format", {".mesh"}, nullopt, WriteFluentFormat);
}
