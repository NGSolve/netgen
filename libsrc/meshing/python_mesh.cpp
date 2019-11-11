#ifdef NG_PYTHON

#include <../general/ngpython.hpp>
#include <core/python_ngcore.hpp>
#include "python_mesh.hpp"

#include <mystdlib.h>
#include "meshing.hpp"
// #include <csg.hpp>
// #include <geometry2d.hpp>
#include <../interface/writeuser.hpp>


using namespace netgen;

extern const char *ngscript[];

namespace netgen
{
  extern bool netgen_executable_started;
  extern shared_ptr<NetgenGeometry> ng_geometry;
  extern void Optimize2d (Mesh & mesh, MeshingParameters & mp);

#ifdef PARALLEL
  /** we need allreduce in python-wrapped communicators **/
  template <typename T>
  inline T MyMPI_AllReduceNG (T d, const MPI_Op & op /* = MPI_SUM */, MPI_Comm comm)
  {
    T global_d;
    MPI_Allreduce ( &d, &global_d, 1, MyGetMPIType<T>(), op, comm);
    return global_d;
  }
#else
  // enum { MPI_SUM = 0, MPI_MIN = 1, MPI_MAX = 2 };
  // typedef int MPI_Op;
  template <typename T>
  inline T MyMPI_AllReduceNG (T d, const MPI_Op & op /* = MPI_SUM */, MPI_Comm comm)
  { return d; }
#endif
}


void TranslateException (const NgException & ex)
{
  string err = string("Netgen exception: ")+ex.What();
  PyErr_SetString(PyExc_RuntimeError, err.c_str());
}

static Transformation<3> global_trafo(Vec<3> (0,0,0));

DLL_HEADER void ExportNetgenMeshing(py::module &m) 
{
  py::register_exception<NgException>(m, "NgException");
  m.attr("_netgen_executable_started") = py::cast(netgen::netgen_executable_started);
  string script;
  const char ** hcp = ngscript;
  while (*hcp)
      script += *hcp++;

  m.attr("_ngscript") = py::cast(script);

  m.def("_GetStatus", []()
        {
          MyStr s; double percent;
          GetStatus(s, percent);
          return py::make_tuple(s.c_str(), percent);
        });
  m.def("_PushStatus", [](string s) { PushStatus(MyStr(s)); });
  m.def("_SetThreadPercentage", [](double percent) { SetThreadPercent(percent); });

  py::class_<NgMPI_Comm> (m, "MPI_Comm")
    .def_property_readonly ("rank", &NgMPI_Comm::Rank)
    .def_property_readonly ("size", &NgMPI_Comm::Size)
    .def("Barrier", &NgMPI_Comm::Barrier)
    
#ifdef PARALLEL
    .def("WTime", [](NgMPI_Comm  & c) { return MPI_Wtime(); })
#else
    .def("WTime", [](NgMPI_Comm  & c) { return -1.0; })
#endif
    .def("Sum", [](NgMPI_Comm  & c, double x) { return MyMPI_AllReduceNG(x, MPI_SUM, c); })
    .def("Min", [](NgMPI_Comm  & c, double x) { return MyMPI_AllReduceNG(x, MPI_MIN, c); })
    .def("Max", [](NgMPI_Comm  & c, double x) { return MyMPI_AllReduceNG(x, MPI_MAX, c); })
    .def("Sum", [](NgMPI_Comm  & c, int x) { return MyMPI_AllReduceNG(x, MPI_SUM, c); })
    .def("Min", [](NgMPI_Comm  & c, int x) { return MyMPI_AllReduceNG(x, MPI_MIN, c); })
    .def("Max", [](NgMPI_Comm  & c, int x) { return MyMPI_AllReduceNG(x, MPI_MAX, c); })
    .def("Sum", [](NgMPI_Comm  & c, size_t x) { return MyMPI_AllReduceNG(x, MPI_SUM, c); })
    .def("Min", [](NgMPI_Comm  & c, size_t x) { return MyMPI_AllReduceNG(x, MPI_MIN, c); })
    .def("Max", [](NgMPI_Comm  & c, size_t x) { return MyMPI_AllReduceNG(x, MPI_MAX, c); })
    .def("SubComm", [](NgMPI_Comm & c, std::vector<int> proc_list) {
        Array<int> procs(proc_list.size());
        for (int i = 0; i < procs.Size(); i++)
          { procs[i] = proc_list[i]; }
        if (!procs.Contains(c.Rank()))
          { throw Exception("rank "+ToString(c.Rank())+" not in subcomm"); }
	return c.SubCommunicator(procs);
      }, py::arg("procs"));
  ;



  
  py::class_<NGDummyArgument>(m, "NGDummyArgument")
    .def("__bool__", []( NGDummyArgument &self ) { return false; } )
    ;
  
  py::class_<Point<2>> (m, "Point2d")
    .def(py::init<double,double>())
    .def ("__str__", &ToString<Point<2>>)
    .def(py::self-py::self)
    .def(py::self+Vec<2>())
    .def(py::self-Vec<2>())
    .def("__getitem__", [](Point<2>& self, int index) { return self[index]; })
    ;

  py::class_<Point<3>> (m, "Point3d")
    .def(py::init<double,double,double>())
    .def ("__str__", &ToString<Point<3>>)
    .def(py::self-py::self)
    .def(py::self+Vec<3>())
    .def(py::self-Vec<3>())
    .def("__getitem__", [](Point<2>& self, int index) { return self[index]; })
    ;

  m.def ("Pnt", FunctionPointer
           ([](double x, double y, double z) { return global_trafo(Point<3>(x,y,z)); }));
  m.def ("Pnt", FunctionPointer
           ([](double x, double y) { return Point<2>(x,y); }));

  /*
    // duplicated functions ????
  m.def ("Pnt", FunctionPointer
           ([](double x, double y, double z) { return Point<3>(x,y,z); }));
  m.def ("Pnt", FunctionPointer
           ([](double x, double y) { return Point<2>(x,y); }));
  */

  py::class_<Vec<2>> (m, "Vec2d")
    .def(py::init<double,double>())
    .def ("__str__", &ToString<Vec<3>>)
    .def(py::self+py::self)
    .def(py::self-py::self)
    .def(-py::self)
    .def(double()*py::self)
    .def("Norm", &Vec<2>::Length)
    .def("__getitem__", [](Vec<2>& vec, int index) { return vec[index]; })
    .def("__len__", [](Vec<2>& /*unused*/) { return 2; })
    ;

  py::class_<Vec<3>> (m, "Vec3d")
    .def(py::init<double,double,double>())
    .def ("__str__", &ToString<Vec<3>>)
    .def(py::self+py::self)
    .def(py::self-py::self)
    .def(-py::self)
    .def(double()*py::self)
    .def("Norm", &Vec<3>::Length)
    .def("__getitem__", [](Vec<3>& vec, int index) { return vec[index]; })
    .def("__len__", [](Vec<3>& /*unused*/) { return 3; })
    ;

  m.def ("Vec", FunctionPointer
           ([] (double x, double y, double z) { return global_trafo(Vec<3>(x,y,z)); }));
  m.def ("Vec", FunctionPointer
           ([] (double x, double y) { return Vec<2>(x,y); }));

  py::class_<Transformation<3>> (m, "Trafo")
    .def(py::init<Vec<3>>(), "a translation")
    .def(py::init<Point<3>,Vec<3>,double>(), "a rotation given by point on axes, direction of axes, angle")
    .def("__mul__", [](Transformation<3> a, Transformation<3> b)->Transformation<3>
         { Transformation<3> res; res.Combine(a,b); return res; })
    .def("__call__", [] (Transformation<3> trafo, Point<3> p) { return trafo(p); })
    ;

  m.def ("GetTransformation", [] () { return global_trafo; });
  m.def ("SetTransformation", [] (Transformation<3> trafo) { global_trafo = trafo; });
  m.def ("SetTransformation", 
         [](int dir, double angle)
         {
           if (dir > 0)
             global_trafo.SetAxisRotation (dir, angle*M_PI/180);
           else
             global_trafo = Transformation<3> (Vec<3>(0,0,0));
         },
         py::arg("dir")=int(0), py::arg("angle")=int(0));
  m.def ("SetTransformation", 
         [](Point<3> p0, Vec<3> ex, Vec<3> ey, Vec<3> ez)
            {
              Point<3> pnts[4];
              pnts[0] = p0;
              pnts[1] = p0 + ex;
              pnts[2] = p0 + ey;
              pnts[3] = p0 + ez;
              global_trafo = Transformation<3> (pnts);
            },
         py::arg("p0"), py::arg("ex"), py::arg("ey"), py::arg("ez"));


  
  py::class_<PointIndex>(m, "PointId")
    .def(py::init<int>())
    .def("__repr__", &ToString<PointIndex>)
    .def("__str__", &ToString<PointIndex>)
    .def_property_readonly("nr", &PointIndex::operator int)
    .def("__eq__" , FunctionPointer( [](PointIndex &self, PointIndex &other)
                  { return static_cast<int>(self)==static_cast<int>(other); }) )
    .def("__hash__" , FunctionPointer( [](PointIndex &self ) { return static_cast<int>(self); }) )
    ;

  py::class_<ElementIndex>(m, "ElementId3D")
    .def(py::init<int>())
    .def("__repr__", &ToString<ElementIndex>)
    .def("__str__", &ToString<ElementIndex>)
    .def_property_readonly("nr", &ElementIndex::operator int)
    .def("__eq__" , FunctionPointer( [](ElementIndex &self, ElementIndex &other)
                  { return static_cast<int>(self)==static_cast<int>(other); }) )
    .def("__hash__" , FunctionPointer( [](ElementIndex &self ) { return static_cast<int>(self); }) )
    ;


  py::class_<SurfaceElementIndex>(m, "ElementId2D")
    .def(py::init<int>())
    .def("__repr__", &ToString<SurfaceElementIndex>)
    .def("__str__", &ToString<SurfaceElementIndex>)
    .def_property_readonly("nr", &SurfaceElementIndex::operator int)
    .def("__eq__" , FunctionPointer( [](SurfaceElementIndex &self, SurfaceElementIndex &other)
                  { return static_cast<int>(self)==static_cast<int>(other); }) )
    .def("__hash__" , FunctionPointer( [](SurfaceElementIndex &self ) { return static_cast<int>(self); }) )
    ;

  py::class_<SegmentIndex>(m, "ElementId1D")
    .def(py::init<int>())
    .def("__repr__", &ToString<SegmentIndex>)
    .def("__str__", &ToString<SegmentIndex>)
    .def_property_readonly("nr", &SegmentIndex::operator int)
    .def("__eq__" , FunctionPointer( [](SegmentIndex &self, SegmentIndex &other)
                  { return static_cast<int>(self)==static_cast<int>(other); }) )
    .def("__hash__" , FunctionPointer( [](SegmentIndex &self ) { return static_cast<int>(self); }) )
    ;



  /*  
  py::class_<Point<3>> ("Point")
    .def(py::init<double,double,double>())
    ;
  */

  py::class_<MeshPoint /* ,py::bases<Point<3>> */ >(m, "MeshPoint")
    .def(py::init<Point<3>>())
    .def("__str__", &ToString<MeshPoint>)
    .def("__repr__", &ToString<MeshPoint>)
    .def_property_readonly("p", [](const MeshPoint & self)
                           {
                             py::list l;
                             l.append ( py::cast(self[0]) );
                             l.append ( py::cast(self[1]) );
                             l.append ( py::cast(self[2]) );
                             return py::tuple(l);
                           })
    .def("__getitem__", [](const MeshPoint & self, int index) {
	  if(index<0 || index>2)
              throw py::index_error();
	  return self[index];
	})
    .def("__setitem__", [](MeshPoint & self, int index, double val) {
	  if(index<0 || index>2)
              throw py::index_error();
	  self(index) = val;
	})
    ;

  py::class_<Element>(m, "Element3D")
    .def(py::init([](int index, std::vector<PointIndex> vertices)
                  {
                    int np = vertices.size();
                    ELEMENT_TYPE et;
                    switch (np)
                      {
                      case 4: et = TET; break;
                      case 5: et = PYRAMID; break;
                      case 6: et = PRISM; break;
                      case 8: et = HEX; break;
                      case 10: et = TET10; break;
                      case 13: et = PYRAMID13; break;
                      case 15: et = PRISM15; break;
                      case 20: et = HEX20; break;
                      default:
                        throw Exception ("no Element3D with " + ToString(np) +
                                         " points");
                      }

                    auto newel = new Element(et);
                    for(int i=0; i<np; i++)
                      (*newel)[i] = vertices[i];
                    newel->SetIndex(index);
                    return newel;
                  }),
          py::arg("index")=1,py::arg("vertices"),
         "create volume element"
         )
    .def("__repr__", &ToString<Element>)
    .def_property("index", &Element::GetIndex, &Element::SetIndex)
    .def_property("curved", &Element::IsCurved, &Element::SetCurved)    
    .def_property_readonly("vertices", 
                  FunctionPointer ([](const Element & self) -> py::list
                                   {
                                     py::list li;
                                     for (int i = 0; i < self.GetNV(); i++)
                                       li.append (py::cast(self[i]));
                                     return li;
                                   }))
    .def_property_readonly("points", 
                  FunctionPointer ([](const Element & self) -> py::list
                                   {
                                     py::list li;
                                     for (int i = 0; i < self.GetNP(); i++)
                                       li.append (py::cast(self[i]));
                                     return li;
                                   }))
    ;

  py::class_<Element2d>(m, "Element2D")
    .def(py::init ([](int index, std::vector<PointIndex> vertices)
                   {
                     Element2d * newel = nullptr;
                     if (vertices.size() == 3)
                       {
                         newel = new Element2d(TRIG);
                         for (int i = 0; i < 3; i++)
                           (*newel)[i] = vertices[i];
                         newel->SetIndex(index);
                       }
                     else if (vertices.size() == 4)
                       {
                         newel = new Element2d(QUAD);
                         for (int i = 0; i < 4; i++)
                           (*newel)[i] = vertices[i];
                         newel->SetIndex(index);
                       }
                     else if (vertices.size() == 6)
                       {
                         newel = new Element2d(TRIG6);
                         for(int i = 0; i<6; i++)
                           (*newel)[i] = vertices[i];
                         newel->SetIndex(index);
                       }
                     else if (vertices.size() == 8)
                       {
                         newel = new Element2d(QUAD8);
                         for(int i = 0; i<8; i++)
                           (*newel)[i] = vertices[i];
                         newel->SetIndex(index);
                       }
                     else 
                       throw NgException("Inconsistent number of vertices in Element2D");
                     return newel;
                   }),
                   py::arg("index")=1,py::arg("vertices"),
         "create surface element"
         )
    .def_property("index", &Element2d::GetIndex, &Element2d::SetIndex)
    .def_property("curved", &Element2d::IsCurved, &Element2d::SetCurved)
    .def_property_readonly("vertices",
                  FunctionPointer([](const Element2d & self) -> py::list
                                  {
                                    py::list li;
                                    for (int i = 0; i < self.GetNV(); i++)
                                      li.append(py::cast(self[i]));
                                    return li;
                                  }))
    .def_property_readonly("points", 
                  FunctionPointer ([](const Element2d & self) -> py::list
                                   {
                                     py::list li;
                                     for (int i = 0; i < self.GetNP(); i++)
                                       li.append (py::cast(self[i]));
                                     return li;
                                   }))
    ;

  py::class_<Segment>(m, "Element1D")
    .def(py::init([](py::list vertices, py::list surfaces, int index, int edgenr)
                  {
                    Segment * newel = new Segment();
                    for (int i = 0; i < 2; i++)
                      (*newel)[i] = py::extract<PointIndex>(vertices[i])();
                    newel -> si = index;
                    newel -> edgenr = edgenr;
                    newel -> epgeominfo[0].edgenr = edgenr;
                    newel -> epgeominfo[1].edgenr = edgenr;
                    // needed for codim2 in 3d
                    newel -> edgenr = index;
                    if (len(surfaces))
                      {
                        newel->surfnr1 = py::extract<int>(surfaces[0])();
                        newel->surfnr2 = py::extract<int>(surfaces[1])();
                      }
                    return newel;
                  }),
          py::arg("vertices"),
           py::arg("surfaces")=py::list(),
           py::arg("index")=1,
           py::arg("edgenr")=1,
         "create segment element"
         )
    .def("__repr__", &ToString<Segment>)
    .def_property_readonly("vertices", 
                  FunctionPointer ([](const Segment & self) -> py::list
                                   {
                                     py::list li;
                                     for (int i = 0; i < 2; i++)
                                       li.append (py::cast(self[i]));
                                     return li;
                                   }))
    .def_property_readonly("points", 
                  FunctionPointer ([](const Segment & self) -> py::list
                                   {
                                     py::list li;
                                     for (int i = 0; i < self.GetNP(); i++)
                                       li.append (py::cast(self[i]));
                                     return li;
                                   }))
    .def_property_readonly("surfaces", 
                  FunctionPointer ([](const Segment & self) -> py::list
                                   {
                                     py::list li;
                                     li.append (py::cast(self.surfnr1));
                                     li.append (py::cast(self.surfnr2));
                                     return li;
                                   }))
    .def_property_readonly("index", FunctionPointer([](const Segment &self) -> size_t
		  {
		    return self.si;
		  }))
    .def_property_readonly("edgenr", FunctionPointer([](const Segment & self) -> size_t
						     {
						       return self.edgenr;
						     }))
    ;


  py::class_<Element0d>(m, "Element0D")
    .def(py::init([](PointIndex vertex, int index)
                  {
                    Element0d * instance = new Element0d;
                    instance->pnum = vertex;
                    instance->index = index;
                    return instance;
                  }),
         py::arg("vertex"),
         py::arg("index")=1,
         "create point element"
         )
    .def("__repr__", &ToString<Element0d>)
    .def_property_readonly("vertices", 
                  FunctionPointer ([](const Element0d & self) -> py::list
                                   {
                                     py::list li;
                                     li.append (py::cast(self.pnum));
                                     return li;
                                   }))
    ;
  
  
  


  py::class_<FaceDescriptor>(m, "FaceDescriptor")
    .def(py::init<const FaceDescriptor&>())
    .def(py::init([](int surfnr, int domin, int domout, int bc)
                  {
                    FaceDescriptor * instance = new FaceDescriptor();
                    instance->SetSurfNr(surfnr);
                    instance->SetDomainIn(domin);
                    instance->SetDomainOut(domout);
                    instance->SetBCProperty(bc);
                             return instance;
                  }),
         py::arg("surfnr")=1, 
         py::arg("domin")=1,
         py::arg("domout")=py::int_(0),
         py::arg("bc")=py::int_(0),
         "create facedescriptor")
    .def("__str__", &ToString<FaceDescriptor>)
    .def("__repr__", &ToString<FaceDescriptor>)
    .def_property("surfnr", &FaceDescriptor::SurfNr, &FaceDescriptor::SetSurfNr)
    .def_property("domin", &FaceDescriptor::DomainIn, &FaceDescriptor::SetDomainIn)
    .def_property("domout", &FaceDescriptor::DomainOut, &FaceDescriptor::SetDomainOut)
    .def_property("bc", &FaceDescriptor::BCProperty, &FaceDescriptor::SetBCProperty)
    .def_property("bcname",
                  [](FaceDescriptor & self) -> string { return self.GetBCName(); },
                  [](FaceDescriptor & self, string name) { self.SetBCName(new string(name)); } // memleak
                  )
    .def("SetSurfaceColor", [](FaceDescriptor & self, py::list color )
          {
            Vec3d c;
            c.X() = py::extract<double>(color[0])();
            c.Y() = py::extract<double>(color[1])();
            c.Z() = py::extract<double>(color[2])();
            self.SetSurfColour(c);
          })
    ;

  

  ExportArray<Element,size_t>(m);
  ExportArray<Element2d,SurfaceElementIndex>(m);
  ExportArray<Segment,size_t>(m);
  ExportArray<Element0d>(m);
  ExportArray<MeshPoint,PointIndex>(m);
  ExportArray<FaceDescriptor>(m);

  py::implicitly_convertible< int, PointIndex>();

  py::class_<NetgenGeometry, shared_ptr<NetgenGeometry>> (m, "NetgenGeometry", py::dynamic_attr())
             ;
  
  py::class_<Mesh,shared_ptr<Mesh>>(m, "Mesh")
    // .def(py::init<>("create empty mesh"))

    .def(py::init( [] (int dim, NgMPI_Comm comm)
                   {
                     auto mesh = make_shared<Mesh>();
		     mesh->SetCommunicator(comm);
                     mesh -> SetDimension(dim);
                     SetGlobalMesh(mesh);  // for visualization
                     mesh -> SetGeometry (nullptr);
                     return mesh;
                   } ),
         py::arg("dim")=3, py::arg("comm")=NgMPI_Comm{}
         )
    .def(NGSPickle<Mesh>())
    .def_property_readonly("comm", [](const Mesh & amesh) -> NgMPI_Comm
			   { return amesh.GetCommunicator(); },
                           "MPI-communicator the Mesh lives in")
    /*
    .def("__init__",
         [](Mesh *instance, int dim)
                           {
                             new (instance) Mesh();
                             instance->SetDimension(dim);
                           },
           py::arg("dim")=3
          )
    */
    
    .def_property_readonly("_timestamp", &Mesh::GetTimeStamp)
    .def_property_readonly("ne", [](Mesh& m) { return m.GetNE(); })
    .def("Distribute", [](shared_ptr<Mesh> self, NgMPI_Comm comm) {
	self->SetCommunicator(comm);
	if(comm.Size()==1) return self;
	// if(MyMPI_GetNTasks(comm)==2) throw NgException("Sorry, cannot handle communicators with NP=2!");
	// cout << " rank " << MyMPI_GetId(comm) << " of " << MyMPI_GetNTasks(comm) << " called Distribute " << endl;
	if(comm.Rank()==0) self->Distribute();
	else self->SendRecvMesh();
	return self;
      }, py::arg("comm"))
    .def_static("Receive", [](NgMPI_Comm comm) -> shared_ptr<Mesh> {
        auto mesh = make_shared<Mesh>();
        mesh->SetCommunicator(comm);
        mesh->SendRecvMesh();
        return mesh;
      }, py::arg("comm"))
    .def("Load",  FunctionPointer 
	 ([](shared_ptr<Mesh> self, const string & filename)
	  {

	    auto comm = self->GetCommunicator();
	    int id = comm.Rank();
	    int ntasks = comm.Size();
	    auto & mesh = self;

	    {
	      ifstream infile(filename.c_str());
	      if(!infile.good())
		throw NgException(string("Error opening file ") + filename);
	    }

	    if ( filename.find(".vol") == string::npos )
	      {
		if(ntasks>1)
		  throw NgException("Not sure what to do with this?? Does this work with MPI??");
		mesh->SetCommunicator(comm);
		ReadFile(*mesh,filename.c_str());
		//mesh->SetGlobalH (mparam.maxh);
		//mesh->CalcLocalH();
		return;
	      }

	    istream * infile;
	    NgArray<char> buf; // for distributing geometry!
	    int strs;

	    if( id == 0) {

	      if (filename.substr (filename.length()-3, 3) == ".gz")
		infile = new igzstream (filename.c_str());
	      else
		infile = new ifstream (filename.c_str());
	      mesh -> Load(*infile);

	      // make string from rest of file (for geometry info!)
	      // (this might be empty, in which case we take the global ng_geometry)
	      stringstream geom_part;
	      geom_part << infile->rdbuf();
	      string geom_part_string = geom_part.str();
	      strs = geom_part_string.size();
	      // buf = new char[strs];
	      buf.SetSize(strs);
	      memcpy(&buf[0], geom_part_string.c_str(), strs*sizeof(char));

	      delete infile;

	      if (ntasks > 1)
		{

		  char * weightsfilename = new char [filename.size()+1];
		  strcpy (weightsfilename, filename.c_str());            
		  weightsfilename[strlen (weightsfilename)-3] = 'w';
		  weightsfilename[strlen (weightsfilename)-2] = 'e';
		  weightsfilename[strlen (weightsfilename)-1] = 'i';

		  ifstream weightsfile(weightsfilename);      
		  delete [] weightsfilename;  
	  
		  if (!(weightsfile.good()))
		    {
		      // cout << "regular distribute" << endl;
		      mesh -> Distribute();
		    }
		  else
		    {
		      char str[20];   
		      bool endfile = false;
		      int n, dummy;
	      
		      NgArray<int> segment_weights;
		      NgArray<int> surface_weights;
		      NgArray<int> volume_weights;
	      
		      while (weightsfile.good() && !endfile)
			{
			  weightsfile >> str;
		  
			  if (strcmp (str, "edgeweights") == 0)
			    {
			      weightsfile >> n;
			      segment_weights.SetSize(n);
			      for (int i = 0; i < n; i++)
				weightsfile >> dummy >> segment_weights[i];
			    }
		  
			  if (strcmp (str, "surfaceweights") == 0)
			    {
			      weightsfile >> n;
			      surface_weights.SetSize(n);
			      for (int i=0; i<n; i++)
				weightsfile >> dummy >> surface_weights[i];
			    }
		  
			  if (strcmp (str, "volumeweights") == 0)
			    {
			      weightsfile >> n;
			      volume_weights.SetSize(n);
			      for (int i=0; i<n; i++)
				weightsfile >> dummy >> volume_weights[i];
			    }
		  
			  if (strcmp (str, "endfile") == 0)
			    endfile = true;  
			}     
	      
		      mesh -> Distribute(volume_weights, surface_weights, segment_weights);
		    }
		} // ntasks>1 end
	    } // id==0 end
	    else {
	      mesh->SendRecvMesh();
	    }

	    if(ntasks>1) {
#ifdef PARALLEL
	      /** Scatter the geometry-string (no dummy-implementation in mpi_interface) **/
	      int strs = buf.Size();
	      MyMPI_Bcast(strs, comm);
	      if(strs>0)
		MyMPI_Bcast(buf, comm);
#endif
	    }

	    shared_ptr<NetgenGeometry> geo;
	    if(buf.Size()) { // if we had geom-info in the file, take it
	      istringstream geom_infile(string((const char*)&buf[0], buf.Size()));
	      geo = geometryregister.LoadFromMeshFile(geom_infile);
	    }
	    if(geo!=nullptr) mesh->SetGeometry(geo);
	    else if(ng_geometry!=nullptr) mesh->SetGeometry(ng_geometry);
	  }),py::call_guard<py::gil_scoped_release>())
    // static_cast<void(Mesh::*)(const string & name)>(&Mesh::Load))
    .def("Save", static_cast<void(Mesh::*)(const string & name)const>(&Mesh::Save),py::call_guard<py::gil_scoped_release>())
    .def("Export",
         [] (Mesh & self, string filename, string format)
          {
            if (WriteUserFormat (format, self, /* *self.GetGeometry(), */ filename))
              {
                string err = string ("nothing known about format")+format;
                NgArray<const char*> names, extensions;
                RegisterUserFormats (names, extensions);
                err += "\navailable formats are:\n";
                for (auto name : names)
                  err += string("'") + name + "'\n";
                throw NgException (err);
              }
          },
         py::arg("filename"), py::arg("format"),py::call_guard<py::gil_scoped_release>())
    
    .def_property("dim", &Mesh::GetDimension, &Mesh::SetDimension)

    .def("Elements3D", 
         static_cast<Array<Element>&(Mesh::*)()> (&Mesh::VolumeElements),
         py::return_value_policy::reference)

    .def("Elements2D", 
         static_cast<Array<Element2d,SurfaceElementIndex>&(Mesh::*)()> (&Mesh::SurfaceElements),
         py::return_value_policy::reference)

    .def("Elements1D", 
         static_cast<Array<Segment>&(Mesh::*)()> (&Mesh::LineSegments),
         py::return_value_policy::reference)

    .def("Elements0D", FunctionPointer([] (Mesh & self) -> Array<Element0d>&
                                       {
                                         return self.pointelements;
                                       } ),
         py::return_value_policy::reference)

    .def("Points", 
         static_cast<Mesh::T_POINTS&(Mesh::*)()> (&Mesh::Points),
         py::return_value_policy::reference)

    .def("FaceDescriptor", static_cast<FaceDescriptor&(Mesh::*)(int)> (&Mesh::GetFaceDescriptor),
         py::return_value_policy::reference)
    .def("GetNFaceDescriptors", &Mesh::GetNFD)

    .def("GetVolumeNeighboursOfSurfaceElement", [](Mesh & self, size_t sel)
                                                {
                                                  int elnr1, elnr2;
                                                  self.GetTopology().GetSurface2VolumeElement(sel+1, elnr1, elnr2);
                                                  return py::make_tuple(elnr1, elnr2);
                                                }, "Returns element nrs of volume element connected to surface element, -1 if no volume element")

    .def("GetNCD2Names", &Mesh::GetNCD2Names)
    

    .def("__getitem__", [](const Mesh & self, PointIndex id) { return self[id]; })
    .def("__getitem__", [](const Mesh & self, ElementIndex id) { return self[id]; })
    .def("__getitem__", [](const Mesh & self, SurfaceElementIndex id) { return self[id]; })
    .def("__getitem__", [](const Mesh & self, SegmentIndex id) { return self[id]; })

    .def("__setitem__", [](Mesh & self, PointIndex id, const MeshPoint & mp) { return self[id] = mp; })
    
    .def ("Add", [](Mesh & self, MeshPoint p)
          {
            return self.AddPoint (Point3d(p));
          })
          
    .def ("Add", [](Mesh & self, const Element & el)
          {
            return self.AddVolumeElement (el);
          })
          
    .def ("Add", [](Mesh & self, const Element2d & el)
          {
            return self.AddSurfaceElement (el);
          })

    .def ("Add", [](Mesh & self, const Segment & el)
          {
            return self.AddSegment (el);
          })
          
    .def ("Add", [](Mesh & self, const Element0d & el)
          {
            return self.pointelements.Append (el);
          })

    .def ("Add", [](Mesh & self, const FaceDescriptor & fd)
          {
            return self.AddFaceDescriptor (fd);
          })
    
    .def ("DeleteSurfaceElement",
          [](Mesh & self, SurfaceElementIndex i)
          {
            return self.Delete(i);
          })
          
    .def ("Compress", [](Mesh & self)
          {
            return self.Compress ();
          } ,py::call_guard<py::gil_scoped_release>())
          
    .def ("AddRegion", [] (Mesh & self, string name, int dim) -> int
         {
           auto & regionnames = self.GetRegionNamesCD(self.GetDimension()-dim);
           regionnames.Append (new string(name));
           int idx = regionnames.Size();
           if (dim == 2)
             {
               FaceDescriptor fd;
               fd.SetBCName(regionnames.Last());
               fd.SetBCProperty(idx);
               self.AddFaceDescriptor(fd);
             }
           return idx;
         }, py::arg("name"), py::arg("dim"))
    
    .def ("SetBCName", &Mesh::SetBCName)
    .def ("GetBCName", FunctionPointer([](Mesh & self, int bc)->string 
                                       { return self.GetBCName(bc); }))
    .def ("SetMaterial", &Mesh::SetMaterial)
    .def ("GetMaterial", FunctionPointer([](Mesh & self, int domnr)
                                         { return string(self.GetMaterial(domnr)); }))

    .def ("GetCD2Name", &Mesh::GetCD2Name)
    .def ("SetCD2Name", &Mesh::SetCD2Name)

    .def ("GetCD3Name", &Mesh::GetCD3Name)
    .def ("SetCD3Name", &Mesh::SetCD3Name)

    .def ("AddPointIdentification", [](Mesh & self, py::object pindex1, py::object pindex2, int identnr, int type)
                           {
			     if(py::extract<PointIndex>(pindex1).check() && py::extract<PointIndex>(pindex2).check())
			       {
				 self.GetIdentifications().Add (py::extract<PointIndex>(pindex1)(), py::extract<PointIndex>(pindex2)(), identnr);
				 self.GetIdentifications().SetType(identnr, Identifications::ID_TYPE(type)); // type = 2 ... periodic
			       }
                           },
          //py::default_call_policies(),
          py::arg("pid1"),
           py::arg("pid2"),
           py::arg("identnr"),
           py::arg("type"))
    .def ("CalcLocalH", &Mesh::CalcLocalH)
    .def ("SetMaxHDomain", [] (Mesh& self, py::list maxhlist)
          {
            NgArray<double> maxh;
            for(auto el : maxhlist)
              maxh.Append(py::cast<double>(el));
            self.SetMaxHDomain(maxh);
          })
    .def ("GenerateVolumeMesh", 
          [](Mesh & self, MeshingParameters* pars,
             py::kwargs kwargs)
           {
             MeshingParameters mp;
             if(pars) mp = *pars;
             {
               py::gil_scoped_acquire acquire;
               CreateMPfromKwargs(mp, kwargs);
             }
             MeshVolume (mp, self);
             OptimizeVolume (mp, self);
           }, py::arg("mp")=nullptr,
          meshingparameter_description.c_str(),
          py::call_guard<py::gil_scoped_release>())

    .def ("OptimizeVolumeMesh", [](Mesh & self)
          {
            MeshingParameters mp;
            mp.optsteps3d = 5;
            OptimizeVolume (mp, self);
          },py::call_guard<py::gil_scoped_release>())

    .def ("OptimizeMesh2d", [](Mesh & self)
          {
            self.CalcLocalH(0.5);
            MeshingParameters mp;
            mp.optsteps2d = 5;
            Optimize2d (self, mp);
          },py::call_guard<py::gil_scoped_release>())
    
    .def ("Refine", FunctionPointer
          ([](Mesh & self)
           {
             self.GetGeometry()->GetRefinement().Refine(self);
             self.UpdateTopology();
           }),py::call_guard<py::gil_scoped_release>())

    .def ("SecondOrder", FunctionPointer
          ([](Mesh & self)
           {
             self.GetGeometry()->GetRefinement().MakeSecondOrder(self);
           }))

    .def ("GetGeometry", [] (Mesh& self) { return self.GetGeometry(); })
    .def ("SetGeometry", [](Mesh & self, shared_ptr<NetgenGeometry> geo)
           {
             self.SetGeometry(geo);
           })

    /*
    .def ("SetGeometry", FunctionPointer
          ([](Mesh & self, shared_ptr<CSGeometry> geo)
           {
             self.SetGeometry(geo);
           }))
    */
    
    .def ("BuildSearchTree", &Mesh::BuildElementSearchTree,py::call_guard<py::gil_scoped_release>())

    .def ("BoundaryLayer", FunctionPointer 
          ([](Mesh & self, int bc, py::list thicknesses, int volnr, py::list materials)
           {
             int n = py::len(thicknesses);
             BoundaryLayerParameters blp;

             for (int i = 1; i <= self.GetNFD(); i++)
               if (self.GetFaceDescriptor(i).BCProperty() == bc)
                   blp.surfid.Append (i);

             cout << "add layer at surfaces: " << blp.surfid << endl;

             blp.prismlayers = n;
             blp.growthfactor = 1.0;

             // find max domain nr
             int maxind = 0;
             for (ElementIndex ei = 0; ei < self.GetNE(); ei++)
               maxind = max (maxind, self[ei].GetIndex());
             cout << "maxind = " << maxind << endl;
             for ( int i=0; i<n; i++ )
               {
                 blp.heights.Append( py::extract<double>(thicknesses[i])()) ;
                 blp.new_matnrs.Append( maxind+1+i );
                 self.SetMaterial (maxind+1+i, py::extract<string>(materials[i])().c_str());
               }
             blp.bulk_matnr = volnr;
             GenerateBoundaryLayer (self, blp);
           }
           ))

    .def ("BoundaryLayer", FunctionPointer
          ([](Mesh & self, int bc, double thickness, int volnr, string material)
           {
             BoundaryLayerParameters blp;

             for (int i = 1; i <= self.GetNFD(); i++)
               if (self.GetFaceDescriptor(i).BCProperty() == bc)
                   blp.surfid.Append (i);

             cout << "add layer at surfaces: " << blp.surfid << endl;

             blp.prismlayers = 1;
             blp.hfirst = thickness;
             blp.growthfactor = 1.0;

             // find max domain nr
             int maxind = 0;
             for (ElementIndex ei = 0; ei < self.GetNE(); ei++)
               maxind = max (maxind, self[ei].GetIndex());
             cout << "maxind = " << maxind << endl;
             self.SetMaterial (maxind+1, material.c_str());
             blp.new_matnr = maxind+1;
             blp.bulk_matnr = volnr;
             GenerateBoundaryLayer (self, blp);
           }
           ))

    .def ("EnableTable", [] (Mesh & self, string name, bool set)
          {
            if (name == "edges")
              const_cast<MeshTopology&>(self.GetTopology()).SetBuildEdges(set);
            if (name == "faces")
              const_cast<MeshTopology&>(self.GetTopology()).SetBuildFaces(set);
          },
          py::arg("name"), py::arg("set")=true)
    
    .def ("Scale", [](Mesh & self, double factor)
          {
            for(auto i = 0; i<self.GetNP();i++)
              self.Point(i).Scale(factor);
          })
    .def ("Copy", [](Mesh & self)
          {
            auto m2 = make_shared<Mesh> ();
            *m2 = self;
            return m2;
          })
    .def ("CalcMinMaxAngle", [](Mesh & self, double badel_limit)
          {
            double values[4];
            self.CalcMinMaxAngle (badel_limit, values);
            py::dict res;
            res["trig"] = py::make_tuple( values[0], values[1] );
            res["tet"] = py::make_tuple( values[2], values[3] );
            return res;
          }, py::arg("badelement_limit")=175.0)
    .def ("CalcTotalBadness", &Mesh::CalcTotalBad)
    .def ("GetQualityHistogram", &Mesh::GetQualityHistogram)
    ;

  m.def("ImportMesh", [](const string& filename)
                      {
                        auto mesh = make_shared<Mesh>();
                        ReadFile(*mesh, filename);
                        return mesh;
                      }, py::arg("filename"),
    R"delimiter(Import mesh from other file format, supported file formats are:
 Neutral format (*.mesh, *.emt)
 Surface file (*.surf)
 Universal format (*.unv)
 Olaf format (*.emt)
 Tet format (*.tet)
 Pro/ENGINEER format (*.fnf)
)delimiter");
  py::enum_<MESHING_STEP>(m,"MeshingStep")
    .value("ANALYSE", MESHCONST_ANALYSE)
    .value("MESHEDGES", MESHCONST_MESHEDGES)
    .value("MESHSURFACE", MESHCONST_OPTSURFACE)
    .value("MESHVOLUME", MESHCONST_OPTVOLUME)
    ;
         
  typedef MeshingParameters MP;
  auto mp = py::class_<MP> (m, "MeshingParameters")
    .def(py::init<>())
            .def(py::init([](MeshingParameters* other, py::kwargs kwargs)
                  {
                    MeshingParameters mp;
                    if(other) mp = *other;
                    CreateMPfromKwargs(mp, kwargs, false);
                    return mp;
                  }), py::arg("mp")=nullptr, meshingparameter_description.c_str())
    .def("__str__", &ToString<MP>)
    .def("RestrictH", FunctionPointer
         ([](MP & mp, double x, double y, double z, double h)
          {
            mp.meshsize_points.Append ( MeshingParameters::MeshSizePoint (Point<3> (x,y,z), h));
          }),
         py::arg("x"), py::arg("y"), py::arg("z"), py::arg("h")
         )
    ;

  m.def("SetTestoutFile", FunctionPointer ([] (const string & filename)
                                             {
                                               delete testout;
                                               testout = new ofstream (filename);
                                             }));

  m.def("SetMessageImportance", FunctionPointer ([] (int importance)
                                                   {
                                                     int old = printmessage_importance;
                                                     printmessage_importance = importance;
                                                     return old;
                                                   }));


}

PYBIND11_MODULE(libmesh, m) {
  ExportNetgenMeshing(m);
}
#endif




