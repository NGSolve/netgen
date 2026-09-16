//
//  Write dolfin file
//
//  by
//  Kent-Andre Mardal <kent-and@simula.no>


#include <mystdlib.h>

#include <myadt.hpp>
#include <linalg.hpp>
#include <csg.hpp>
#include <meshing.hpp>

namespace netgen
{

#include "writeuser.hpp"



  void WriteDolfinFormat (const Mesh & mesh, const filesystem::path & filename)
  {
    cout << "start writing dolfin export" << endl;

    int np = mesh.GetNP();
    int ne = mesh.GetNE();
    // int nse = mesh.GetNSE();
    int nsd = mesh.GetDimension(); 
    // int invertsurf = mparam.inverttrigs;
    // int i, j;

    ofstream outfile (filename);

    // char str[100];
    outfile.precision(8);
    outfile.setf (ios::fixed, ios::floatfield);
    outfile.setf (ios::showpoint);

    if ( nsd == 3) {

      outfile << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" <<endl; 
      outfile << ""<<endl; 

      outfile << "<dolfin xmlns:dolfin=\"http://www.phi.chalmers.se/dolfin/\">"<<endl;
      outfile << "  <mesh celltype=\"tetrahedron\" dim=\"3\">" <<endl; 
      outfile << "      <vertices size=\""<<np<<"\">"<<endl; 
      for (PointIndex pi : mesh.Points().Range()) {
        const Point<3> & p = mesh[pi];
        outfile << "      <vertex index=\""<<pi.Nr0()<<"\" x=\""<<p(0)<<"\" y=\""<<p(1)<<"\" z=\""<<p(2)<<"\"/>"<<endl; 
      }
      outfile << "      </vertices>"<<endl; 



      outfile << "      <cells size=\""<<ne<<"\">"<<endl; 
      for (ElementIndex i : T_Range<ElementIndex>(ne)) {
        const Element & el = mesh[i];

        outfile << "      <tetrahedron index=\""<<i.Nr1()-1<<"\" v0=\""<<el[0]-1<<"\" v1=\""<<el[1]-1<<"\" v2=\""<<el[2]-1<<"\" v3=\""<<el[3]-1<<"\"/>"<<endl; 
      }
      outfile << "      </cells>"<<endl; 
    }
    outfile << "   </mesh>"<<endl; 
    outfile << "</dolfin>"<<endl; 

    cout << "done writing dolfin export" << endl;
  }
}
