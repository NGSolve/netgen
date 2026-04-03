//
//  Write Abaqus file
//
//

#include <mystdlib.h>

#include <myadt.hpp>
#include <linalg.hpp>
#include <csg.hpp>
#include <meshing.hpp>

#include "writeuser.hpp"

namespace netgen
{

using std::vector;

struct AbaqusElementType
{
  const char * name;
  const vector<int> permutation;

  AbaqusElementType(const char * name, const vector<int> & permutation)
    : name(name), permutation(permutation)
  {}
};

static inline const AbaqusElementType & GetAbaqusType(int dim, int num_nodes)
{
  // maps num_nodes to AbaqusElementType for each dimension
  typedef std::map<int, AbaqusElementType> AbaqusElementTypes;
  static const std::map<int, AbaqusElementType> abaqus_eltypes[3] =
  {
    // 1D
    AbaqusElementTypes{
      {2, AbaqusElementType{"T2D2", vector{0,1}}},
    },
    // 2D
    AbaqusElementTypes{
      {3, AbaqusElementType{"CPS3", vector{0,1,2}}},
      {6, AbaqusElementType{"CPS6", vector{0,1,2,5,6,4}}},      
    },
    // 3D
    AbaqusElementTypes{
      {4, AbaqusElementType{"C3D4", vector{0,1,3,2}}},
      {10, AbaqusElementType{"C3D10", vector{0,1,3,2,4,8,6,5,7,9}}},
    }
  };

  const auto & eltypes = abaqus_eltypes[dim-1];
  if (eltypes.count(num_nodes) > 0)
    return eltypes.at(num_nodes);
  else
    throw Exception("unsupported " + ToString(dim)+"d Element type with " + ToString(num_nodes) + " nodes");
}

static void WritePoints ( const Mesh & mesh, ostream & out )
{
  out << "*Node" << endl;
  for(auto pi : mesh.Points().Range() )
  {
    out << pi+1-IndexBASE<PointIndex>() << ", ";
    auto p = mesh[pi];
    out << p[0] << ", " << p[1] << ", " << p[2] << '\n';
  }
}

template <typename ElIndex>
static void WriteElement(ostream & out, const Mesh& mesh, ElIndex ei, const vector<int> & permutation, int & el_counter)
{
  el_counter++;
  auto el = mesh[ei];
  out << el_counter;
  for(auto i : Range(el.PNums()))
    out << ", " << el[permutation[i]]+1-IndexBASE<PointIndex>();
  out << '\n';
}

template <typename ElIndex, typename Elements>
static void WriteElements ( ostream & out, const Mesh & mesh, int dim, const Elements & el_range, int & el_counter)
{
  // map index, num_nodes to elements
  std::map<std::tuple<int, int>, Array<ElIndex>> elset_map;

  for(auto ei : el_range)
    {
      const auto & el = mesh[ei];
      int index = 0;
      if constexpr(std::is_same_v<ElIndex,SegmentIndex>)
        index = el.edgenr;
      else
        index = el.GetIndex();
      elset_map[{index, el.GetNP()}].Append(ei);
    }

  for(auto & [key, elems] : elset_map)
    {
      auto [index, num_nodes] = key;
      auto name = mesh.GetRegionName(elems[0]);
      if (name == "") name = "default";
      PrintMessage (5, index, ": ", name);
      const auto & eltype = GetAbaqusType(dim, num_nodes) ;
      out << "*Element, type=" << eltype.name << ", ELSET=" << name << endl;
      for(auto ei : elems)
        WriteElement(out, mesh, ei, eltype.permutation, el_counter);
    }
}

void WriteAbaqusFormat (const Mesh & mesh,
			const filesystem::path & filename)

{
  PrintMessage (1, "Write Abaqus Mesh");

  ofstream outfile (filename);

  outfile << "*Heading" << endl;
  outfile << " " << filename << endl;

  outfile.precision(8);

  int element_counter = 0;
  WritePoints(mesh, outfile);
  if(mesh.GetDimension() < 3)
    WriteElements<SegmentIndex>(outfile, mesh, 1, mesh.LineSegments().Range(), element_counter);
  WriteElements<SurfaceElementIndex>(outfile, mesh, 2, mesh.SurfaceElements().Range(), element_counter);
  WriteElements<ElementIndex>(outfile, mesh, 3, mesh.VolumeElements().Range(), element_counter);

  // Write identifications (untested!)
  if (mesh.GetIdentifications().GetMaxNr())
    {
      const auto np = mesh.GetNP();
      // periodic identification, implementation for
      // Helmut J. Boehm, TU Vienna
	  
      auto mpcfilename = filename;
      if (filename.extension() == ".inp")
        mpcfilename.replace_extension(".mpc");
      else
        mpcfilename.concat(".mpc");
	  
      ofstream mpc (mpcfilename);

      int masternode(0);

      NgArray<INDEX_2> pairs;
      NgBitArray master(np), help(np);
      master.Set();
      for (int i = 1; i <= 3; i++)
	{
	  mesh.GetIdentifications().GetPairs (i, pairs);
	  help.Clear();
	  for (int j = 1; j <= pairs.Size(); j++)
	    {
	      help.Set (pairs.Get(j).I1());
	    }
	  master.And (help);
	}
      for (int i = 1; i <= np; i++)
	if (master.Test(i))
	  masternode = i;

      cout << "masternode = " << masternode << " = "
	   << mesh.Point(masternode) << endl;
      NgArray<int> minions(3);
      for (int i = 1; i <= 3; i++)
	{
	  mesh.GetIdentifications().GetPairs (i, pairs);
	  for (int j = 1; j <= pairs.Size(); j++)
	    {
	      if (pairs.Get(j).I1() == masternode)
		minions.Elem(i) = pairs.Get(j).I2();
	    }
	  cout << "minion(" << i << ") = " << minions.Get(i)
	       << " = " << mesh.Point(minions.Get(i)) << endl;
	}
	  
	  
      outfile << "**\n"
	      << "*NSET,NSET=CTENODS\n"
	      << minions.Get(1) << ", " 
	      << minions.Get(2) << ", " 
	      << minions.Get(3) << endl;

	  
      outfile << "**\n"
	      << "**POINT_fixed\n"
	      << "**\n"
	      << "*BOUNDARY, OP=NEW\n";
      for (int j = 1; j <= 3; j++)
	outfile << masternode << ", " << j << ",,    0.\n";

      outfile << "**\n"
	      << "*BOUNDARY, OP=NEW\n";
      for (int j = 1; j <= 3; j++)
	{
	  Vec3d v(mesh.Point(masternode), mesh.Point(minions.Get(j)));
	  double vlen = v.Length();
	  int dir = 0;
	  if (fabs (v.X()) > 0.9 * vlen) dir = 2;
	  if (fabs (v.Y()) > 0.9 * vlen) dir = 3;
	  if (fabs (v.Z()) > 0.9 * vlen) dir = 1;
	  if (!dir)
	    cout << "ERROR: Problem with rigid body constraints" << endl;
	  outfile << minions.Get(j) << ", " << dir << ",,    0.\n";
	}

      outfile << "**\n"
	      << "*EQUATION, INPUT=" << mpcfilename << endl;
	  

      NgBitArray eliminated(np);
      eliminated.Clear();
      for (int i = 1; i <= mesh.GetIdentifications().GetMaxNr(); i++)
	{
	  mesh.GetIdentifications().GetPairs (i, pairs);
	  if (!pairs.Size())
	    continue;
	      
	  for (int j = 1; j <= pairs.Size(); j++)
	    if (pairs.Get(j).I1() != masternode && 
		!eliminated.Test(pairs.Get(j).I2()))
	      {
		eliminated.Set (pairs.Get(j).I2());
		for (int k = 1; k <= 3; k++)
		  {
		    mpc << "4" << "\n";
		    mpc << pairs.Get(j).I2() << "," << k << ", -1.0, ";
		    mpc << pairs.Get(j).I1() << "," << k << ", 1.0, ";
		    mpc << minions.Get(i) << "," << k << ", 1.0, ";
		    mpc << masternode << "," << k << ", -1.0 \n";
		  }
	      }
	}
    }


  PrintMessage(1, "done");
}

static RegisterUserFormat reg_abaqus ("Abaqus Format", {".mesh"}, nullopt, WriteAbaqusFormat);
}
