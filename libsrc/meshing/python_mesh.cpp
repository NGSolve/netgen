#ifdef NG_PYTHON

#include <../general/ngpython.hpp>

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
}


template <typename T, int BASE = 0, typename TIND = int>
void ExportArray (py::module &m)
{
  using TA = Array<T,BASE,TIND>;
  string name = string("Array_") + typeid(T).name();
  py::class_<Array<T,BASE,TIND>>(m, name.c_str())
    .def ("__len__", [] ( Array<T,BASE,TIND> &self ) { return self.Size(); } )
    .def ("__getitem__", 
          FunctionPointer ([](Array<T,BASE,TIND> & self, TIND i) -> T&
                           {
                             if (i < BASE || i >= BASE+self.Size())
                               throw py::index_error();
                             return self[i];
                           }),
          py::return_value_policy::reference)
    .def("__iter__", [] ( TA & self) {
	return py::make_iterator (self.begin(),self.end());
      }, py::keep_alive<0,1>()) // keep array alive while iterator is used

    ;
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

  py::class_<NGDummyArgument>(m, "NGDummyArgument")
    .def("__bool__", []( NGDummyArgument &self ) { return false; } )
    ;
  
  py::class_<Point<2>> (m, "Point2d")
    .def(py::init<double,double>())
    .def ("__str__", &ToString<Point<2>>)
    .def(py::self-py::self)
    .def(py::self+Vec<2>())
    .def(py::self-Vec<2>())
    ;

  py::class_<Point<3>> (m, "Point3d")
    .def(py::init<double,double,double>())
    .def ("__str__", &ToString<Point<3>>)
    .def(py::self-py::self)
    .def(py::self+Vec<3>())
    .def(py::self-Vec<3>())
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
    ;

  py::class_<Vec<3>> (m, "Vec3d")
    .def(py::init<double,double,double>())
    .def ("__str__", &ToString<Vec<3>>)
    .def(py::self+py::self)
    .def(py::self-py::self)
    .def(-py::self)
    .def(double()*py::self)
    .def("Norm", &Vec<3>::Length)
    ;

  m.def ("Vec", FunctionPointer
           ([] (double x, double y, double z) { return global_trafo(Vec<3>(x,y,z)); }));
  m.def ("Vec", FunctionPointer
           ([] (double x, double y) { return Vec<2>(x,y); }));

  py::class_<Transformation<3>> (m, "Trafo")
    .def(py::init<Vec<3>>())
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
    .def_property_readonly("p", FunctionPointer([](const MeshPoint & self)
                                       {
                                         py::list l;
                                         l.append ( py::cast(self[0]) );
                                         l.append ( py::cast(self[1]) );
                                         l.append ( py::cast(self[2]) );
                                         return py::tuple(l);
                                       }))
    .def("__getitem__", FunctionPointer([](const MeshPoint & self, int index) {
	  if(index<0 || index>2)
              throw py::index_error();
	  return self[index];
	}))
    .def("__setitem__", FunctionPointer([](MeshPoint & self, int index, double val) {
	  if(index<0 || index>2)
              throw py::index_error();
	  self(index) = val;
	}))
    ;
  
  py::class_<Element>(m, "Element3D")
    .def(py::init([](int index, py::list vertices)
                  {
                    Element * newel = nullptr;
                    if (py::len(vertices) == 4)
                      {
                        newel = new Element(TET);
                        for (int i = 0; i < 4; i++)
                          (*newel)[i] = py::extract<PointIndex>(vertices[i])();
                        newel->SetIndex(index);
                      }
                    else if (py::len(vertices) == 5)
                      {
                        newel = new Element(PYRAMID);
                        for (int i = 0; i < 5; i++)
                          (*newel)[i] = py::extract<PointIndex>(vertices[i])();
                        newel->SetIndex(index);
                      }
                    else if (py::len(vertices) == 6)
                      {
                        newel = new Element(PRISM);
                        for (int i = 0; i < 6; i++)
                          (*newel)[i] = py::extract<PointIndex>(vertices[i])();
                        newel->SetIndex(index);
                      }
                    else if (py::len(vertices) == 8)
                      {
                        newel = new Element(HEX);
                        for (int i = 0; i < 8; i++)
                          (*newel)[i] = py::extract<PointIndex>(vertices[i])();
                        newel->SetIndex(index);
                      }
                    else
                      throw NgException ("cannot create element");
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
    .def(py::init ([](int index, py::list vertices)
                   {
                     Element2d * newel = nullptr;
                     if (py::len(vertices) == 3)
                       {
                         newel = new Element2d(TRIG);
                         for (int i = 0; i < 3; i++)
                           (*newel)[i] = py::extract<PointIndex>(vertices[i])();
                         newel->SetIndex(index);
                       }
                     else if (py::len(vertices) == 4)
                       {
                         newel = new Element2d(QUAD);
                         for (int i = 0; i < 4; i++)
                           (*newel)[i] = py::extract<PointIndex>(vertices[i])();
                         newel->SetIndex(index);
                       }
                     else if (py::len(vertices) == 6)
                       {
                         newel = new Element2d(TRIG6);
                         for(int i = 0; i<6; i++)
                           (*newel)[i] = py::extract<PointIndex>(vertices[i])();
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

  

  ExportArray<Element,0,size_t>(m);
  ExportArray<Element2d,0,size_t>(m);
  ExportArray<Segment,0,size_t>(m);
  ExportArray<Element0d>(m);
  ExportArray<MeshPoint,PointIndex::BASE,PointIndex>(m);
  ExportArray<FaceDescriptor>(m);

  py::implicitly_convertible< int, PointIndex>();

  py::class_<NetgenGeometry, shared_ptr<NetgenGeometry>> (m, "NetgenGeometry", py::dynamic_attr())
             ;
  
  py::class_<Mesh,shared_ptr<Mesh>>(m, "Mesh")
    // .def(py::init<>("create empty mesh"))

    .def(py::init( [] (int dim)
                   {
                     auto mesh = make_shared<Mesh>();
                     mesh -> SetDimension(dim);
                     SetGlobalMesh(mesh);  // for visualization
                     mesh -> SetGeometry (nullptr);
                     return mesh;
                   } ),
         py::arg("dim")=3         
         )
    .def(NGSPickle<Mesh>())

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
    .def("Load",  FunctionPointer 
	 ([](Mesh & self, const string & filename)
	  {
	    istream * infile;

#ifdef PARALLEL
	    MPI_Comm_rank(MPI_COMM_WORLD, &id);
	    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
	    
	    char* buf = nullptr;
	    int strs = 0;
	    if(id==0) {
#endif
            if (filename.find(".vol.gz") != string::npos)
              infile = new igzstream (filename.c_str());
            else
              infile = new ifstream (filename.c_str());
	    // ifstream input(filename);
#ifdef PARALLEL
	    //still inside id==0-bracket...
	        self.Load(*infile);	      
		self.Distribute();

		/** Copy the rest of the file into a string (for geometry) **/
		stringstream geom_part;
		geom_part << infile->rdbuf();
		string geom_part_string = geom_part.str();
		strs = geom_part_string.size();
		buf = new char[strs];
		memcpy(buf, geom_part_string.c_str(), strs*sizeof(char));
	      }
	    else {
	      self.SendRecvMesh();
	    }

	    /** Scatter the geometry-string **/
	    MPI_Bcast(&strs, 1, MPI_INT, 0, MPI_COMM_WORLD); 
	    if(id!=0)
	      buf = new char[strs];
	    MPI_Bcast(buf, strs, MPI_CHAR, 0, MPI_COMM_WORLD);
	    if(id==0)
	      delete infile;
	    infile = new istringstream(string((const char*)buf, (size_t)strs));
	    delete[] buf;
	    
#else
	    self.Load(*infile);
#endif
	    for (int i = 0; i < geometryregister.Size(); i++)
	      {
		NetgenGeometry * hgeom = geometryregister[i]->LoadFromMeshFile (*infile);
		if (hgeom)
		  {
		    ng_geometry.reset (hgeom);
                    self.SetGeometry(ng_geometry);
		    break;
		  }
	      }
	    self.SetGeometry(ng_geometry);
	    delete infile;
	  }),py::call_guard<py::gil_scoped_release>())
    // static_cast<void(Mesh::*)(const string & name)>(&Mesh::Load))
    .def("Save", static_cast<void(Mesh::*)(const string & name)const>(&Mesh::Save),py::call_guard<py::gil_scoped_release>())
    .def("Export",
         [] (Mesh & self, string filename, string format)
          {
            if (WriteUserFormat (format, self, /* *self.GetGeometry(), */ filename))
              {
                string err = string ("nothing known about format")+format;
                Array<const char*> names, extensions;
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
         static_cast<Array<Element,0,size_t>&(Mesh::*)()> (&Mesh::VolumeElements),
         py::return_value_policy::reference)

    .def("Elements2D", 
         static_cast<Array<Element2d,0,size_t>&(Mesh::*)()> (&Mesh::SurfaceElements),
         py::return_value_policy::reference)

    .def("Elements1D", 
         static_cast<Array<Segment,0,size_t>&(Mesh::*)()> (&Mesh::LineSegments),
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

    .def("GetNCD2Names", &Mesh::GetNCD2Names)
    

    .def("__getitem__", FunctionPointer ([](const Mesh & self, PointIndex pi)
                                         {
                                           return self[pi];
                                         }))

    .def ("Add", FunctionPointer ([](Mesh & self, MeshPoint p)
                                  {
                                    return self.AddPoint (Point3d(p));
                                  }))

    .def ("Add", FunctionPointer ([](Mesh & self, const Element & el)
                                  {
                                    return self.AddVolumeElement (el);
                                  }))

    .def ("Add", FunctionPointer ([](Mesh & self, const Element2d & el)
                                  {
                                    return self.AddSurfaceElement (el);
                                  }))

    .def ("Add", FunctionPointer ([](Mesh & self, const Segment & el)
                                  {
                                    return self.AddSegment (el);
                                  }))
    
    .def ("Add", FunctionPointer ([](Mesh & self, const Element0d & el)
                                  {
                                    return self.pointelements.Append (el);
                                  }))

    .def ("Add", FunctionPointer ([](Mesh & self, const FaceDescriptor & fd)
                                  {
                                    return self.AddFaceDescriptor (fd);
                                  }))
    
    .def ("DeleteSurfaceElement",
          FunctionPointer ([](Mesh & self, SurfaceElementIndex i)
                           {
                             return self.DeleteSurfaceElement (i);
                           }))
    
    .def ("Compress", FunctionPointer ([](Mesh & self)
                                       {
                                         return self.Compress ();
                                       }),py::call_guard<py::gil_scoped_release>())
    
    
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
            Array<double> maxh;
            for(auto el : maxhlist)
              maxh.Append(py::cast<double>(el));
            self.SetMaxHDomain(maxh);
          })
    .def ("GenerateVolumeMesh", 
          [](Mesh & self, py::object pymp)
           {
             cout << "generate vol mesh" << endl;

             MeshingParameters mp;
             {
               py::gil_scoped_acquire acquire;
             if (py::extract<MeshingParameters>(pymp).check())
               mp = py::extract<MeshingParameters>(pymp)();
             else
               {
                 mp.optsteps3d = 5;
               }
             }
             MeshVolume (mp, self);
             OptimizeVolume (mp, self);
           },
          py::arg("mp")=NGDummyArgument(),py::call_guard<py::gil_scoped_release>())

   .def ("OptimizeVolumeMesh", FunctionPointer
         ([](Mesh & self)
          {
            MeshingParameters mp;
            mp.optsteps3d = 5;
            OptimizeVolume (mp, self);
          }),py::call_guard<py::gil_scoped_release>())

    .def ("Refine", FunctionPointer
          ([](Mesh & self)
           {
             if (self.GetGeometry())
               self.GetGeometry()->GetRefinement().Refine(self);
             else
               Refinement().Refine(self);
             self.UpdateTopology();
           }),py::call_guard<py::gil_scoped_release>())

    .def ("SecondOrder", FunctionPointer
          ([](Mesh & self)
           {
             if (self.GetGeometry())
               self.GetGeometry()->GetRefinement().MakeSecondOrder(self);
             else
               Refinement().MakeSecondOrder(self);
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
    
    .def ("Scale", FunctionPointer([](Mesh & self, double factor)
				   {
				     for(auto i = 0; i<self.GetNP();i++)
				       self.Point(i).Scale(factor);
				   }))
                                            
    ;
  

  py::enum_<MESHING_STEP>(m,"MeshingStep")
    .value("MESHEDGES",MESHCONST_MESHEDGES)
    .value("MESHSURFACE",MESHCONST_OPTSURFACE)
    .value("MESHVOLUME",MESHCONST_OPTVOLUME)
    ;
         
  typedef MeshingParameters MP;
  py::class_<MP> (m, "MeshingParameters")
    .def(py::init<>())
    .def(py::init([](double maxh, double minh, bool quad_dominated, int optsteps2d, int optsteps3d,
                     MESHING_STEP perfstepsend, int only3D_domain, const string & meshsizefilename,
                     double grading, double curvaturesafety, double segmentsperedge)
                  {
                    MP * instance = new MeshingParameters;
                    instance->maxh = maxh;
                    instance->minh = minh;
                    instance->quad = int(quad_dominated);
                    instance->optsteps2d = optsteps2d;
                    instance->optsteps3d = optsteps3d;			     
                    instance->only3D_domain_nr = only3D_domain;
                    instance->perfstepsend = perfstepsend;
                    instance->meshsizefilename = meshsizefilename;
                    
                    instance->grading = grading;
                    instance->curvaturesafety = curvaturesafety;
                    instance->segmentsperedge = segmentsperedge;
                    return instance;
                  }),
         py::arg("maxh")=1000,
         py::arg("minh")=0,
         py::arg("quad_dominated")=false,
         py::arg("optsteps2d") = 3,
	 py::arg("optsteps3d") = 3,
	 py::arg("perfstepsend") = MESHCONST_OPTVOLUME,
	 py::arg("only3D_domain") = 0,
         py::arg("meshsizefilename") = "",
         py::arg("grading")=0.3,
         py::arg("curvaturesafety")=2,
         py::arg("segmentsperedge")=1,
         "create meshing parameters"
         )
    .def("__str__", &ToString<MP>)
    .def_property("maxh", 
                  FunctionPointer ([](const MP & mp ) { return mp.maxh; }),
                  FunctionPointer ([](MP & mp, double maxh) { return mp.maxh = maxh; }))
    .def_property("minh",
                  FunctionPointer ([](const MP & mp ) { return mp.minh; }),
                  FunctionPointer ([](MP & mp, double minh) { return mp.minh = minh; }))
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




