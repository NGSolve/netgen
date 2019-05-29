/*************************************
 * Write Gmsh file
 * First issue the 04/26/2004 by Paul CARRICO (paul.carrico@free.fr)
 * At the moment, the GMSH format is available for
 * linear tetrahedron elements i.e. in 3D
 * (based on Neutral Format)
 *
 * Second issue the 05/05/2004 by Paul CARRICO
 * Thanks to Joachim Schoeberl for the correction of a minor bug
 * the 2 initial Gmsh Format (i.e. volume format and surface format) are group together)
 * in only one file
 **************************************/

#include <mystdlib.h>

#include <myadt.hpp>
#include <linalg.hpp>
#include <csg.hpp>
#include <meshing.hpp>

namespace netgen
{
#include "writeuser.hpp"

  extern MeshingParameters mparam;


/*
 *  VTK mesh format
 *  points, elements, surface elements
 */

void WriteVtkFormat (const Mesh & mesh,
                      const string & filename)
{
  ofstream outfile (filename.c_str());
  outfile.precision(6);
  outfile.setf (ios::fixed, ios::floatfield);
  outfile.setf (ios::showpoint);

  int np = mesh.GetNP();  /// number of point
  int ne = mesh.GetNE();  /// number of element
  int nse = mesh.GetNSE();  /// number of surface element (BC)

  outfile << "# vtk DataFile Version 2.0\n";
  outfile << "Created with netgen\n";
  outfile << "ASCII\n";
  outfile << "DATASET UNSTRUCTURED_GRID\n";
  
  outfile << "POINTS " << np << " double\n";
  for (int i=0; i<np; i++)
  {
  auto & p = mesh.Point(i+1);
    outfile << p[0] << " " << p[1] << " " << p[2] << "\n";
  }

  std::vector<int> types;
  std::vector<int> domains;
  if (ne > 0)
  {
    unsigned int size = 0;
    for (int i=0; i<ne; i++)
      size += mesh.VolumeElement(i+1).GetNV() + 1; // only save "linear" corners

    outfile << "CELLS " << ne << " " << size << "\n";
    for (int i=0; i<ne; i++)
    {
      auto& el = mesh.VolumeElement(i+1);
      domains.push_back(el.GetIndex());
      switch (el.GetType())
      {
      case TET:
      case TET10: // reorder to follow VTK convention & zero based indices
        outfile << 4 << " " << el[0]-1 << " " << el[1]-1 << " " << el[3]-1 << " " << el[2]-1 << "\n";
        types.push_back(10);
        break;
      case PRISM: // reorder to follow VTK convention & zero based indices
        outfile << 6 << " "
          << el[0]-1 << " " << el[2]-1 << " " << el[1]-1 << " " 
          << el[3]-1 << " " << el[5]-1 << " " << el[4]-1 << "\n";
        types.push_back(13);
        break;
      default:
        throw ngcore::Exception("Unexpected element type");
        break;
      }
    }
  }
  else
  {
    unsigned int size = 0;
    for (int i=0; i<nse; i++)
      size += mesh.SurfaceElement(i+1).GetNV() + 1;

    outfile << "CELLS " << nse << " " << size << "\n";
    for (int i=0; i<nse; i++)
    {
      auto& el = mesh.SurfaceElement(i+1);
      domains.push_back(el.GetIndex());
      switch (el.GetType())
      {
      case TRIG:
      case TRIG6:
        outfile << 3 << " " << el[0]-1 << " " << el[1]-1 << " " << el[2]-1 << "\n";
        types.push_back(5);
        break;
      case QUAD:
        outfile << 4 << " " << el[0]-1 << " " << el[1]-1 << " " << el[2]-1 << " " << el[3]-1 << "\n";
        types.push_back(9);
        break;
      default:
        throw ngcore::Exception("Unexpected element type");
      break;
      }
    }
  }

  outfile << "CELL_TYPES " << types.size() << "\n";
  for (auto type_id: types)
  {
    outfile << type_id << "\n";
  }

  outfile << "CELL_DATA " << domains.size() << "\n";
  outfile << "SCALARS scalars int 1\n";
  outfile << "LOOKUP_TABLE default\n";
  for (auto id: domains)
  {
    outfile << id << "\n";
  }

  outfile.close();
}
}


