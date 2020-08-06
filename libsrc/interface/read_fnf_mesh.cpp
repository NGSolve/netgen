
//
//  Read Pro/ENGINEER neutral format
//

#include <mystdlib.h>

#include <myadt.hpp>
#include <linalg.hpp>
#include <csg.hpp>
#include <meshing.hpp>
#include <sys/stat.h>


namespace netgen
{
#include "writeuser.hpp"

  bool ReadLine (istream & in, string & buf)
  {
    do
      {
        buf = "";
        
        while (in.good())
          {
            char ch = in.get();
            if (ch == '\n') break;
            if (ch == '\r') break;
            if (ch == '\\')
              {
                // while (iswhite (ch = in.get() )
                ch = in.get();   // '\n'   CR
                ch = in.get();   // '\n'   LF
              }
            else
              buf += ch;
          }
      }
    while (in.good() && (buf == "" || buf[0] == '#'));
    
    return in.good();
  }
  



  
  class LoadType
  {
  public:
    int id;
    string name;
    string placement;
    string valuetype;
    NgArray<double> places;
  };



  
  void ReadFNFFormat (Mesh & mesh, 
                      const string & filename)
  {
    ifstream fin (filename.c_str());

    string buf;

    mesh.SetDimension (3);

    while (ReadLine (fin, buf))
      {
        stringstream sbuf(buf);
        string start_sect, token; char ch;
        
        sbuf >> start_sect;

        if (start_sect == "%START_SECT")
          {
            sbuf >> ch >> token;
            
            if (token == "HEADER")
              {
                while (1)
                  {
                    ReadLine (fin, buf);
                    stringstream sbuf(buf);
                    string token;
                    
                    sbuf >> token;
                    
                    if (token == "%TITLE")
                      {
                        char ch;
                        string name;
                        sbuf >> ch >> name;
                        PrintMessage(3, "Title: ", name);
                      }
                    else if (token == "%STATISTICS")
                      {
                        ;
                      }
                    else if (token == "%END_SECT")
                      {
                        break;
                      }
                    else
                      {
                        PrintMessage(1, "SECTION HEADER, unknown field: ",buf);
                      }
                  }                
              }
 

            else if (token == "ELEM_TYPES")
              {
                while (1)
                  {
                    ReadLine (fin, buf);
                    stringstream sbuf(buf);
                    string token;
                    
                    sbuf >> token;
                    
                    if (token == "%ELEM_TYPE")
                      {
			int nr;
			string def;
                        char ch;
                        sbuf >> nr >> def >> ch;
			if (def == "DEF")
			  {
			    string classname, type;
			    sbuf >> classname >> type;
			    if (classname != "SOLID" || type != "TETRA")
			      cerr << "Element not supported: " << buf << endl;
			  }
                      }
                    else if (token == "%END_SECT")
                      {
                        break;
                      }
                    else
                      {
                        PrintMessage(1, "SECTION ELEM_TYPE, unknown field: ", buf);
                      }
                  }                
              }
 

            else if (token == "COORD_SYSTEMS")
              {
                while (1)
                  {
                    ReadLine (fin, buf);
                    stringstream sbuf(buf);
                    string token;
                    
                    sbuf >> token;
                    
                    if (token == "%END_SECT")
                      {
                        break;
                      }
                    else
                      {
                        // cout << "COORD_SYSTEMS, unknown field: " << buf << endl;
                      }
                  }                
              }

            else if (token == "ANALYSIS")
              {
                // ignore this section
                while (1)
                  {
                    ReadLine (fin, buf);
                    stringstream sbuf(buf);
                    string token;
                    sbuf >> token;
                    if (token == "%END_SECT")
                      {
                        break;
                      }
                  }                
              }
 
            else if (token == "MATERIALS")
              {
		*testout << "parse materials" << endl;
                NgArray<double> young_modulus, poisson_ratio, mass_density;

                while (1)
                  {
                    ReadLine (fin, buf);
                    stringstream sbuf(buf);
                    string token;
                    
                    sbuf >> token;
                    
                    if (token == "%MATERIAL")
                      {
                        int nr;
                        string prop;
                        char ch;
                        double val;

                        sbuf >> nr >> prop >> ch;
                        if (prop == "DEF")
                          {
                            string name;
                            sbuf >> name;
                            mesh.SetMaterial(nr, name);
                          }
                        else
                          {
                            sbuf >> val;
			    *testout << "prop = " << prop << ", val = " << val << endl;
                            if (prop == "YOUNG_MODULUS")
                              young_modulus.Append (val);
                            else if  (prop == "POISSON_RATIO")
                              poisson_ratio.Append (val);
                            else if  (prop == "MASS_DENSITY")
                              mass_density.Append (val);
                          }
                      }
                    else if (token == "%END_SECT")
                      {
                        mesh.SetUserData ("YOUNG_MODULUS", young_modulus);
                        mesh.SetUserData ("POISSON_RATIO", poisson_ratio);
                        mesh.SetUserData ("MASS_DENSITY", mass_density);
			*testout << "young = " << young_modulus << endl;
			*testout << "poisson = " << poisson_ratio << endl;
                        break;
                      }
                    else
                      {
                        PrintMessage(1, "SECTION MATERIALS, unknown field: ", buf);
                      }
                  }
              }

            
            else if (token == "MESH")
              {
                while (1)
                  {
                    ReadLine (fin, buf);
                    stringstream sbuf(buf);
                    string token;
                    sbuf >> token;
                    if (token == "%NODE")
                      {
                        string st;
                        char ch;
                        int nr, ks_id;
                        double x,y,z;
                        sbuf >> nr >> st >> ch >> x >> y >> z >> ks_id;
                        mesh.AddPoint (Point3d (x,y,z) );
                      }
                    else if (token == "%ELEM")
                      {
                        string elemid, def;
                        char ch;
                        int elnr, typid, matid;
                        string propid;
                        sbuf >> elnr >> def >> ch;
                        sbuf >> typid >> matid >> propid;
                        NgArray<int> pnums;
                        while (1)
                          {
                            int pn;
                            sbuf >> pn;
                            if (!sbuf.good()) break;
                            pnums.Append (pn);
                          }
                        int pe2ng [] = { 0, 1, 2, 3, 4, 7, 5, 6,  8, 9 };
                        Element el(pnums.Size());
                        for (int j = 0; j < pnums.Size(); j++)
                          el[pe2ng[j]] = pnums[j];
                        el.SetIndex (matid);
                        mesh.AddVolumeElement (el);
                      }
                    else if (token == "%END_SECT")
                      {
                        break;
                      }
                    else
                      {
                        PrintMessage(1, "SECTION MESH, unknown: ", buf);
                      }
                  }
              }
            else if (token == "MESH_TOPOLOGY")
              {
                while (1)
                  {
                    ReadLine (fin, buf);
                    stringstream sbuf(buf);
                    string token, kw;
                    int nr;
                    char ch;

                    sbuf >> token;
                    if (token == "%EDGE")
                      {
                        sbuf >> nr >> kw >> ch;
                        if (kw == "NODES")
                          {
                            NgArray<int> enums;
                            while (1)
                              {
                                int en;
                                sbuf >> en;
                                if (!sbuf.good()) break;
                                enums.Append (en);
                              }
                            for (int j = 0; j+2 < enums.Size(); j+=2)
                              {
                                Segment seg;
                                seg[0] = enums[j];
                                seg[1] = enums[j+2];
                                seg[2] = enums[j+1];
                                seg.edgenr = nr;
                                mesh.AddSegment (seg);
                              }
                          }
                      }
                    else if (token == "%SURFACE")
                      {
                        sbuf >> nr >> kw >> ch;
                        if (kw == "FACES")
                          {
                            NgArray<int> fnums;
                            while (1)
                              {
                                int fn;
                                sbuf >> fn;
                                if (!sbuf.good()) break;
                                fnums.Append (fn);
                              }                            

                            FaceDescriptor fd(-1, -1, -1, -1);
                            fd.SetBCProperty (nr);
			    *testout << "add fd " << mesh.GetNFD() << ", nr = " << nr << endl;
                            mesh.AddFaceDescriptor (fd);
                              
                            for (int j = 0; j < fnums.Size(); j += 2)
                              {
                                int elnr = fnums[j];
                                int fnr = fnums[j+1];
                                
                                const Element & el = mesh.VolumeElement (elnr);
                                if(j == 0)
                                  mesh.GetFaceDescriptor(nr).SetDomainIn(el.GetIndex());
                                Element2d el2d;
                                el.GetFace (fnr, el2d);
                                el2d.SetIndex (nr);
                                  
                                mesh.AddSurfaceElement (el2d);
                              }
                          }
                      }
                    else if (token == "%END_SECT")
                      {
                        break;
                      }
                    else
                      {
                        PrintMessage(1, "SECTION MESH, unknown: ", buf);
                      }
                  }
              }



 
            else if (token == "LOADS")
              {
                NgArray<LoadType*> loadtypes;

                while (1)
                  {
                    ReadLine (fin, buf);
                    stringstream sbuf(buf);
                    string token;
                    
                    sbuf >> token;
                    
                    if (token == "%LOAD_TYPE")
                      {
                        string def;
                        char ch;

                        LoadType * lt = new LoadType;
                        sbuf >> lt->id >> def >> ch >> lt->name >> lt->placement >> lt->valuetype;
                        
                        if (lt->name == "DISPLACEMENT")
                          PrintMessage(3, "loadtype DISPLACEMENT found");

                        if (lt->placement != "FACE" && lt->placement != "EDGE" && lt->placement != "NODE")
                          PrintMessage(1, "unsupported placement ", lt->placement);

                        loadtypes.Append (lt);
                      }

                    else if (token == "%LOAD")
                      {
                        int id;
                        string def;
                        char ch;
                        int placement;
                        int load_type_id, con_case_id;
                        sbuf >> id >> def >> ch;
                        
                        if (def == "DEF")
                          {
                            sbuf >> load_type_id >> con_case_id;
                          }
                        if (def == "VAL")
                          {
                            sbuf >> placement;
                            for (int i = 0; i < loadtypes.Size(); i++)
                              if (load_type_id == loadtypes[i]->id)
                                loadtypes[i]->places.Append (placement);
                          }
                      }
                    else if (token == "%CON_CASE")
                      { ; }
                    
                    else if (token == "%END_SECT")
                      {
                        for (int i = 0; i < loadtypes.Size(); i++)
                          {
                            stringstream str;
                            str << loadtypes[i]->places;
                            if (loadtypes[i]->placement == "FACE" && loadtypes[i]->name == "DISPLACEMENT")
                              {
                                mesh.SetUserData ("CONSTRAINT_DISP_FACE", loadtypes[i]->places);
                                PrintMessage(3, "constrained faces: ", str.str());
                              }
                            if (loadtypes[i]->placement == "EDGE" && loadtypes[i]->name == "DISPLACEMENT")
                              {
                                mesh.SetUserData ("CONSTRAINT_DISP_EDGE", loadtypes[i]->places);
                                PrintMessage(3,"constrained edges: ", str.str());
                              }
                            if (loadtypes[i]->placement == "NODE" && loadtypes[i]->name == "DISPLACEMENT")
                              {
                                mesh.SetUserData ("CONSTRAINT_DISP_NODE", loadtypes[i]->places);
                                PrintMessage(3, "constrained nodes: ", str.str());
                              }
                          }
                        break;
                      }
                    else
                      {
                        PrintMessage(1, "SECTION LOADS, unknown field: ", buf);
                      }
                  }
              }



            else
              {
                PrintMessage(1, "unknown section ", token);
              }
          }
        else
          PrintMessage(3, "parse line: (", buf, ")");
      }
  }
}
