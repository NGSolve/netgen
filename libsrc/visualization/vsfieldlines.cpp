#ifndef NOTCL

#include <mystdlib.h>
#include <incopengl.hpp>


#include <myadt.hpp>
#include <meshing.hpp>
#include <csg.hpp>
#include <stlgeom.hpp>

#include <visual.hpp>
#include <meshing/fieldlines.hpp>


namespace netgen
{

  // extern shared_ptr<Mesh> mesh;

  void VisualSceneSolution :: BuildFieldLinesFromBox(Array<Point<3>> & startpoints)
  {
    shared_ptr<Mesh> mesh = GetMesh();
    if (!mesh) return;

    if(fieldlines_startarea_parameter[0] > fieldlines_startarea_parameter[3] ||
       fieldlines_startarea_parameter[1] > fieldlines_startarea_parameter[4] ||
       fieldlines_startarea_parameter[2] > fieldlines_startarea_parameter[5])
      {
	Point3d pmin, pmax;
	mesh->GetBox (pmin, pmax);
	
	fieldlines_startarea_parameter[0] = pmin.X();	
	fieldlines_startarea_parameter[1] = pmin.Y();
	fieldlines_startarea_parameter[2] = pmin.Z();
	fieldlines_startarea_parameter[3] = pmax.X();	
	fieldlines_startarea_parameter[4] = pmax.Y();
	fieldlines_startarea_parameter[5] = pmax.Z();
      }
    
    for (int i = 1; i <= startpoints.Size(); i++)
      {
	Point<3> p (fieldlines_startarea_parameter[0] + double (rand()) / RAND_MAX * (fieldlines_startarea_parameter[3]-fieldlines_startarea_parameter[0]),
		   fieldlines_startarea_parameter[1] + double (rand()) / RAND_MAX * (fieldlines_startarea_parameter[4]-fieldlines_startarea_parameter[1]),
		   fieldlines_startarea_parameter[2] + double (rand()) / RAND_MAX * (fieldlines_startarea_parameter[5]-fieldlines_startarea_parameter[2]));
	
	startpoints[i-1] = p;
      }
  }

  void VisualSceneSolution :: BuildFieldLinesFromLine(Array<Point<3>> & startpoints)
  {
    shared_ptr<Mesh> mesh = GetMesh();
    if (!mesh) return;


    for (int i = 1; i <= startpoints.Size(); i++)
      {
	double s = double (rand()) / RAND_MAX;

	Point<3> p (fieldlines_startarea_parameter[0] + s * (fieldlines_startarea_parameter[3]-fieldlines_startarea_parameter[0]),
		   fieldlines_startarea_parameter[1] + s * (fieldlines_startarea_parameter[4]-fieldlines_startarea_parameter[1]),
		   fieldlines_startarea_parameter[2] + s * (fieldlines_startarea_parameter[5]-fieldlines_startarea_parameter[2]));
	
	startpoints[i-1] = p;
      }
  }


  void VisualSceneSolution :: BuildFieldLinesFromFile(Array<Point<3>> & startpoints)
  {
    shared_ptr<Mesh> mesh = GetMesh();
    if (!mesh) return;

    ifstream * infile;

    infile = new ifstream(fieldlines_filename.c_str());

    //cout << "reading from file " << fieldlines_filename << endl;

    int numpoints = 0;

    string keyword;

    
    double dparam;
    int iparam;

    while(infile->good())
      {
	(*infile) >> keyword;

	if(keyword == "point") numpoints++;
	else if(keyword == "line" || keyword == "box")
	  {
	    for(int i=0; i<6; i++) (*infile) >> dparam;
	    (*infile) >> iparam;
	    numpoints += iparam;
	  }
      }

    delete infile;


    //cout << numpoints << " startpoints" << endl;

    startpoints.SetSize(numpoints);
    
    infile = new ifstream(fieldlines_filename.c_str());

    numpoints = 0;

    while(infile->good())
      {
	(*infile) >> keyword;

	if (keyword == "point")
	  {
	    (*infile) >> startpoints[numpoints][0];
            (*infile) >> startpoints[numpoints][1];
            (*infile) >> startpoints[numpoints][2];
	    numpoints++;
	  }
	else if (keyword == "line" || keyword == "box")
	  {
	    for(int i=0; i<6; i++) (*infile) >> fieldlines_startarea_parameter[i];
	    (*infile) >> iparam;

	    Array<Point<3>> auxpoints(iparam);
	    
	    if (keyword == "box")
	      BuildFieldLinesFromBox(auxpoints);
	    else if (keyword == "line")
	      BuildFieldLinesFromLine(auxpoints);
	    
	    for(int i=0; i<iparam; i++)
	      {
		startpoints[numpoints] = auxpoints[i];
		numpoints++;
	      }
	  }

	//cout << "startpoints " << startpoints << endl;
      }

    delete infile;
    
    

    
  }

  
  void VisualSceneSolution :: BuildFieldLinesFromFace(Array<Point<3>> & startpoints)
  {
    shared_ptr<Mesh> mesh = GetMesh();
    if (!mesh) return;

    Array<SurfaceElementIndex> elements_2d;
    
    //cout << "fieldlines_startface " << fieldlines_startface << endl;
    mesh->GetSurfaceElementsOfFace(fieldlines_startface,elements_2d);
    if(elements_2d.Size() == 0)
      {
	cerr << "No Elements on selected face (?)" << endl;
	return;
      }
    Vec3d v1,v2,cross;
    
    double area = 0;

	int i;
    for(i=0; i<elements_2d.Size(); i++)
      {
	const Element2d & elem = (*mesh)[elements_2d[i]];
	
	v1 = mesh->Point(elem[1]) - mesh->Point(elem[0]);
	v2 = mesh->Point(elem[2]) - mesh->Point(elem[0]);
	cross = Cross(v1,v2);
	area += cross.Length();
	
	if(elem.GetNV() == 4)
	  {
	    v1 = mesh->Point(elem[2]) - mesh->Point(elem[0]);
	    v2 = mesh->Point(elem[3]) - mesh->Point(elem[0]);
	    cross = Cross(v1,v2);
	    area += cross.Length();
	  }
      }
    
    int startpointsp = 0;
    i = 0;
    
    while(startpointsp < startpoints.Size())
      {
	const Element2d & elem = (*mesh)[elements_2d[i]];
	
	int numtri = (elem.GetNV() == 3) ? 1 : 2;
	
	for(int tri = 0; startpointsp < startpoints.Size() && tri<numtri; tri++)
	  {
	    
	    if(tri == 0)
	      {
		v1 = mesh->Point(elem[1]) - mesh->Point(elem[0]);
		v2 = mesh->Point(elem[2]) - mesh->Point(elem[0]);
		cross = Cross(v1,v2);
	      }
	    else if(tri == 1)
	      {
		v1 = mesh->Point(elem[2]) - mesh->Point(elem[0]);
		v2 = mesh->Point(elem[3]) - mesh->Point(elem[0]);
		cross = Cross(v1,v2);
	      }
	    
	    double thisarea = cross.Length();
	    
	    int numloc = int(startpoints.Size()*thisarea/area);
	    if(double (rand()) / RAND_MAX < startpoints.Size()*thisarea/area - numloc)
	      numloc++;
	    
	    for(int j=0; startpointsp < startpoints.Size() && j<numloc; j++)
	      {
		double s = double (rand()) / RAND_MAX;
		double t = double (rand()) / RAND_MAX;
		if(s+t > 1)
		  {
		    s = 1.-s; t = 1.-t;
		  }
		startpoints[startpointsp] = mesh->Point(elem[0]) + s*v1 +t*v2;
		startpointsp++;
	      }
	  }
	i++;
	if(i == elements_2d.Size()) i = 0;
      } 
    
  }


  void VisualSceneSolution :: BuildFieldLinesPlot ()
  {
    shared_ptr<Mesh> mesh = GetMesh();
    if (!mesh) return;

    if (fieldlinestimestamp >= solutiontimestamp) 
      return;
    fieldlinestimestamp = solutiontimestamp;
    

    if (fieldlineslist)
      glDeleteLists (fieldlineslist, num_fieldlineslists);

    if (vecfunction == -1)
      return;

    const SolData * vsol = soldata[fieldlines_vecfunction];

    num_fieldlineslists = (vsol -> iscomplex && !fieldlines_fixedphase) ? 100 : 1;
   
    double phaser=1.0;
    double phasei=0.0;
    std::function<bool(int, const double *, Vec<3> &)> eval_func = [&](int elnr, const double * lami, Vec<3> & vec)
    {
        double values[6] = {0., 0., 0., 0., 0., 0.};
        bool drawelem;
        auto mesh = GetMesh();
        if (mesh->GetDimension()==3)
            drawelem = GetValues (vsol, elnr, lami[0], lami[1], lami[2], values);
        else
            drawelem = GetSurfValues (vsol, elnr, -1, lami[0], lami[1], values);

        Vec3d v;
        RealVec3d (values, v, vsol->iscomplex, phaser, phasei);
        vec = v;
        return drawelem;
    };

    FieldLineCalc linecalc(*mesh, eval_func,
			   fieldlines_rellength,fieldlines_maxpoints,fieldlines_relthickness,fieldlines_reltolerance,fieldlines_rktype);

    if(fieldlines_randomstart) 
      linecalc.Randomized();

    fieldlineslist = glGenLists (num_fieldlineslists);

    int num_startpoints = num_fieldlines / num_fieldlineslists;
    if (num_fieldlines % num_fieldlineslists != 0) num_startpoints++;

    if(fieldlines_randomstart)
      num_startpoints *= 10;

    
    Array<Point<3>> startpoints(num_startpoints);
    

    for (int ln = 0; ln < num_fieldlineslists; ln++)
      {
	if(fieldlines_startarea == 0)
	  BuildFieldLinesFromBox(startpoints);
	else if(fieldlines_startarea == 1)
	  BuildFieldLinesFromFile(startpoints);
	else if(fieldlines_startarea == 2)
	  BuildFieldLinesFromFace(startpoints);


	    
	double phi;
	
	if(vsol -> iscomplex)
	  {
	    if(fieldlines_fixedphase)
	      phi = fieldlines_phase;
	    else
	      phi = 2*M_PI*ln / num_fieldlineslists;
	  }
	else
	  phi = 0;

	cout << "phi = " << phi << endl;

	phaser = cos(phi);
        phasei = sin(phi);
	

	linecalc.GenerateFieldLines(startpoints,num_fieldlines / num_fieldlineslists+1);

        auto & pstart = linecalc.GetPStart();
        auto & pend = linecalc.GetPEnd();
        auto & values = linecalc.GetValues();
        auto nlines = values.Size();

        glNewList(fieldlineslist+ln, GL_COMPILE);
        SetTextureMode (usetexture);

        for(auto i : Range(nlines))
          {
            SetOpenGlColor  (values[i]);
            DrawCylinder (pstart[i], pend[i], fieldlines_relthickness);
          }

        glEndList ();
      }
  }



  
}


#endif // NOTCL
