// Write Chemnitz file format


#include <mystdlib.h>

#include <myadt.hpp>

#include <linalg.hpp>
#include <csg.hpp>
#include <meshing.hpp>

#include "writeuser.hpp"

namespace netgen
{

  class POINT3D
  {
  public:
    POINT3D () { };
    double x, y, z;
  };

  class VOLELEMENT
  {
  public:
    int domnr;
    PointIndex p1, p2, p3, p4;
    int faces[4];

    VOLELEMENT () 
    { for (int i = 0; i < 4; i++) faces[i] = 0; }
  };
  
  class SURFELEMENT
  {
  public:
    SURFELEMENT () { };
    int snr;
    PointIndex p1, p2, p3;
  };
  

  class FACE
  {
  public:
    PointIndex p1, p2, p3;
    int edges[3];

    FACE () 
    { for (int i = 0; i < 3; i++) edges[i] = 0; }
  };

  class EDGE
  {
  public:
    EDGE () { };
    PointIndex p1, p2;
  };

  static Array<POINT3D> points;
  static Array<VOLELEMENT> volelements;
  static Array<SURFELEMENT> surfelements;

  static Array<FACE> faces;
  static Array<EDGE> edges;


  void ReadFile (char * filename)
  {
    int n;
    ifstream infile(filename);
    char reco[100];
  
  
    infile >> reco;  // file format recognition
  
    infile >> n;   // number of surface elements
    cout << n << " Surface elements" << endl;
  
    for (int i = 1; i <= n; i++)
      {
        SURFELEMENT sel;
        infile >> sel.snr >> sel.p1 >> sel.p2 >> sel.p3;
        surfelements.Append (sel);
      }
    
    infile >> n;   // number of volume elements
    cout << n << " Volume elements" << endl;
  
    for (int i = 1; i <= n; i++)
      {
        VOLELEMENT el;
        infile >> el.p1 >> el.p2 >> el.p3 >> el.p4;
        volelements.Append (el);
      }
    
    infile >> n;   // number of points 
    cout << n << " Points" << endl;
  
    for (int i = 1; i <= n; i++)
      {
        POINT3D p;
        infile >> p.x >> p.y >> p.z;
        points.Append (p);
      }
  }
  
  

  void ReadFileMesh (const Mesh & mesh)
  {
    int n = mesh.GetNSE();   // number of surface elements
    cout << n << " Surface elements" << endl;
  
    for (int i = 1; i <= n; i++)
      {
        SURFELEMENT sel;
        const Element2d & el = mesh.SurfaceElement(i);
        sel.snr = el.GetIndex();
        sel.p1 = el.PNum(1);
        sel.p2 = el.PNum(2);
        sel.p3 = el.PNum(3);
        surfelements.Append (sel);
      }
    
    n = mesh.GetNE();   // number of volume elements
    cout << n << " Volume elements" << endl;
  
    for (int i = 1; i <= n; i++)
      {
        VOLELEMENT el;
        const Element & nel = mesh.VolumeElement(i);
        el.p1 = nel.PNum(1);
        el.p2 = nel.PNum(2);
        el.p3 = nel.PNum(3);
        el.p4 = nel.PNum(4);
        //      infile >> el.p1 >> el.p2 >> el.p3 >> el.p4;
        volelements.Append (el);
      }
    
    n = mesh.GetNP();   // number of points 
    cout << n << " Points" << endl;
  
    for (PointIndex pi : mesh.Points().Range())
      {
        POINT3D p;
        const auto & mp = mesh[pi];
        p.x = mp(0);
        p.y = mp(1);
        p.z = mp(2);
        //      infile >> p.x >> p.y >> p.z;
        points.Append (p);
      }
  }
  



  void Convert ()
  {
    INDEX_3_HASHTABLE<int> faceindex(volelements.Size()/5 + 1);
    INDEX_2_HASHTABLE<int> edgeindex(volelements.Size()/5 + 1);

    // face j of a tet is the one opposite to its point j
    static const int facepoints[4][3] = { {1,2,3}, {0,2,3}, {0,1,3}, {0,1,2} };

    for (int i = 1; i <= volelements.Size(); i++)
      {
        const auto & vel = volelements[i-1];
        PointIndex vp[4] = { vel.p1, vel.p2, vel.p3, vel.p4 };

        for (int j = 0; j < 4; j++)
          {
            SortedPointIndices<3> i3 (vp[facepoints[j][0]],
                                      vp[facepoints[j][1]],
                                      vp[facepoints[j][2]]);
            int facei;
            if (faceindex.Used (i3))
              facei = faceindex.Get(i3);
            else
              {
                FACE fa;
                auto [fp1, fp2, fp3] = i3;
                fa.p1 = fp1; fa.p2 = fp2; fa.p3 = fp3;
                faces.Append (fa);
                facei = faces.Size();
                faceindex.Set (i3, facei);
              }

            volelements[i-1].faces[j] = facei;
          }
      }

    // edge j of a face is the one opposite to its point j
    static const int edgepoints[3][2] = { {1,2}, {0,2}, {0,1} };

    for (int i = 1; i <= faces.Size(); i++)
      {
        PointIndex fp[3] = { faces[i-1].p1, faces[i-1].p2, faces[i-1].p3 };

        for (int j = 0; j < 3; j++)
          {
            SortedPointIndices<2> i2 (fp[edgepoints[j][0]], fp[edgepoints[j][1]]);
            int edgei;
            if (edgeindex.Used (i2))
              edgei = edgeindex.Get(i2);
            else
              {
                EDGE ed;
                auto [ep1, ep2] = i2;
                ed.p1 = ep1; ed.p2 = ep2;
                edges.Append (ed);
                edgei = edges.Size();
                edgeindex.Set (i2, edgei);
              }

            faces[i-1].edges[j] = edgei;
          }
      }
  }


  void WriteFile (ostream & outfile)
  {
    outfile 
      << "#VERSION: 1.0" << endl
      << "#PROGRAM: NETGEN" << endl
      << "#EQN_TYPE: POISSON" << endl
      << "#DIMENSION: 3D" << endl
      << "#DEG_OF_FREE: 1" << endl
      << "#DESCRIPTION: I don't know" << endl
      << "##RENUM: not done" << endl
      << "#USER: Kleinzen" << endl
      << "DATE: 10.06.1996" << endl;
  
    outfile << "#HEADER:   8" << endl
            << points.Size() << "  " << edges.Size() << "  " 
            << faces.Size() << "  " << volelements.Size() << "  0  0  0  0" << endl;
  
    outfile << "#VERTEX:   " << points.Size() << endl;
    for (int i = 1; i <= points.Size(); i++)
      outfile << "  " << i << "  " << points[i-1].x << "  " << points[i-1].y 
              << "  " << points[i-1].z << endl;
    	
    outfile << "#EDGE:  " << edges.Size() << endl;
    for (int i = 1; i <= edges.Size(); i++)
      outfile << "  " << i << "  1  " 
              << edges[i-1].p1 << "  " 
              << edges[i-1].p2 
              << "  0" << endl;
    
    outfile << "#FACE:  " << faces.Size() << endl;  
    for (int i = 1; i <= faces.Size(); i++)
      outfile << "  " << i << "  1  3  " 
              << faces[i-1].edges[0] << "  " 
              << faces[i-1].edges[1] << "  " 
              << faces[i-1].edges[2] << endl;
    	
    outfile << "#SOLID:  " << volelements.Size() << endl;
    for (int i = 1; i <= volelements.Size(); i++)
      outfile << "  " << i << "  1  4  " 
              << volelements[i-1].faces[0] << "  "
              << volelements[i-1].faces[1] << "  "
              << volelements[i-1].faces[2] << "  "
              << volelements[i-1].faces[3] << endl;
    	
    outfile << "#END_OF_DATA" << endl;
  }
    

  void WriteUserChemnitz (const Mesh & mesh,
                          const filesystem::path & filename)
  {
    ofstream outfile (filename);

    ReadFileMesh (mesh);
    Convert ();
  
    WriteFile (outfile);
    cout << "Wrote Chemnitz standard file" << endl;
  }
static RegisterUserFormat reg_chemnitz ("Chemnitz Format", {"*"}, nullopt, WriteUserChemnitz );
}
