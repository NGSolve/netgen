#include <mystdlib.h>
#include <myadt.hpp>
#include <meshing.hpp>

#include <visual.hpp>
// #include <parallel.hpp>



#ifndef WIN32
#define GLX_GLXEXT_LEGACY

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xatom.h>  /* for XA_RGB_DEFAULT_MAP atom */
// #include <GL/glx.h>    // for parallel GL ???
#endif





namespace netgen
{
  NGGUI_API Point3d VisualScene :: center;
  NGGUI_API double VisualScene :: rad;
  NGGUI_API GLdouble VisualScene :: backcolor;
  NGGUI_API VisualScene visual_scene_cross;
  NGGUI_API VisualScene *visual_scene = &visual_scene_cross;

  /*
#if TOGL_MAJOR_VERSION!=2
  GLuint VisualScene :: fontbase = 0;
#else
  Tcl_Obj * VisualScene :: fontbase = NULL;
  Togl * VisualScene :: globtogl;
#endif
  */

  void (*opengl_text_function)(const char * text) = NULL;
  int opengl_text_width = 0;
  void Set_OpenGLText_Callback ( void (*fun) (const char * text), int width )
  {
    opengl_text_function = fun;
    opengl_text_width = width;
  }

  void MyOpenGLText (const char * text)
  {
    if (opengl_text_function)
      (*opengl_text_function) (text);
    // cout << "MyOpenGLText: " << text << endl;
  }

  int MyOpenGLTextWidth ()
  {
      return opengl_text_width;
  }


  // texture for color decoding
  // GLubyte * VisualScene :: colortexture = NULL;
  GLuint VisualScene :: coltexname = 1;
  int VisualScene :: ntexcols = -1;


  double VisualScene :: lookatmat[16];
  double VisualScene :: transmat[16];
  double VisualScene :: rotmat[16];
  double VisualScene :: centermat[16];
  double VisualScene :: transformationmat[16];

  int VisualScene :: selface;
  int VisualScene :: selelement;
  PointIndex VisualScene :: selpoint;
  PointIndex VisualScene :: selpoint2;
  int VisualScene :: locpi;
  int VisualScene :: seledge;

  optional<Point<3>> VisualScene :: marker = nullopt;

  int VisualScene :: subdivision_timestamp = -1;
  int VisualScene :: subdivisions = 2;

  int VisualScene :: viewport[4];

  VisualizationParameters :: VisualizationParameters()
  {
    lightamb = 0.3;
    lightdiff = 0.7;
    lightspec = 1;
    shininess = 50;
    transp = 0.3;
    locviewer = 0;
    showstltrias = 0;
    centerpoint = PointIndex::INVALID;
    usedispllists = 1;
    strcpy (selectvisual, "cross");

    use_center_coords = false;
  };
  VisualizationParameters vispar;



  double dist = 0;
  // double dist = 6;
  // vorher: pnear = 2;
  // double pnear = 0.1;
  // double pfar = 10;



  VisualScene :: VisualScene ()
  {
    changeval = -1;
    backcolor = 0;
  }


  VisualScene :: ~VisualScene()
  {
    ;
  }



  void VisualScene :: BuildScene (int zoomall)
  {
    center = Point3d (0,0,0);
    rad = 1;

    CalcTransformationMatrices();

    glEnable(GL_DEPTH_TEST);
    glDisable (GL_DITHER);
  
    GLfloat ambvals[] = { 0.4f, 0.4f, 0.4f, 1.0f };
    GLfloat diffvals[] = { 0.5f, 0.5f, 0.5f, 1.0f };
    GLfloat specvals[] =  { 0.7f, 0.7f, 0.7f, 1.0f };
    glLightfv(GL_LIGHT0, GL_AMBIENT, ambvals);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, diffvals);
    glLightfv(GL_LIGHT0, GL_SPECULAR, specvals);
  
    GLfloat light_position[] = { 1, 3, 3, 0 };
    glLightfv(GL_LIGHT0, GL_POSITION, light_position);
  
    glLightModeli (GL_LIGHT_MODEL_TWO_SIDE, 0);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
  }


  void VisualScene :: DrawScene ()
  {
    if (changeval == -1)
      BuildScene();
    changeval = 0;

    glClearColor(backcolor, backcolor, backcolor, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glEnable (GL_COLOR_MATERIAL);
    glColor3f (1.0f, 1.0f, 1.0f);
    glLineWidth (1.0f);

    DrawCoordinateCross ();
    DrawNetgenLogo ();
    glFinish();  
  }


  void VisualScene :: CalcTransformationMatrices()
  {
    // prepare model view matrix
  
    glPushMatrix();

    glLoadIdentity();
    gluLookAt (0, 0, 6, 0, 0, 0, 0, 1, 0);
    glGetDoublev (GL_MODELVIEW_MATRIX, lookatmat);

    glLoadIdentity();
    glTranslatef(0.0f, 0.0f, -dist);
    glGetDoublev (GL_MODELVIEW_MATRIX, transmat);
  
    glLoadIdentity();
    glGetDoublev (GL_MODELVIEW_MATRIX, rotmat);

    glScaled (1/rad, 1/rad, 1/rad);
    glTranslated (-center.X(), -center.Y(), -center.Z());
    glGetDoublev (GL_MODELVIEW_MATRIX, centermat);

    glLoadIdentity();
    glMultMatrixd (lookatmat);
    glMultMatrixd (transmat);
    glMultMatrixd (rotmat);
    glMultMatrixd (centermat);
    glGetDoublev (GL_MODELVIEW_MATRIX, transformationmat);

    glPopMatrix();
  }


  void VisualScene :: ArbitraryRotation (const NgArray<double> & alpha, const NgArray<Vec3d> & vec)
  {
    glPushMatrix();

    glLoadIdentity();

    for(int i=0; i<alpha.Size() && i<vec.Size(); i++)
      {
	glRotatef(alpha[i], vec[i].X(), vec[i].Y(), vec[i].Z());
      }

    glGetDoublev (GL_MODELVIEW_MATRIX, rotmat);

    glLoadIdentity();
    glMultMatrixd (lookatmat);
    glMultMatrixd (transmat);
    glMultMatrixd (rotmat);
    glMultMatrixd (centermat);
    glGetDoublev (GL_MODELVIEW_MATRIX, transformationmat);
  
    glPopMatrix();
  } 



  void VisualScene :: ArbitraryRotation (const double alpha, const Vec3d & vec)
  {
    NgArray<double> a(1); a[0] = alpha;
    NgArray<Vec3d> v(1); v[0] = vec;

    ArbitraryRotation(a,v);
  } 

  void VisualScene :: StandardRotation (const char * dir)
  {
    glPushMatrix();

    glLoadIdentity();
  
    if (strcmp (dir, "xy") == 0)
      ;
    else if (strcmp (dir, "yx") == 0)
      glRotatef(180.0, 1.0f, 1.0f, 0.0f);    
    else if (strcmp (dir, "xz") == 0)
      glRotatef(-90.0, 1.0f, 0.0f, 0.0f);    
    else if (strcmp (dir, "zx") == 0)
      {
	glRotatef(180.0, 1.0f, 1.0f, 0.0f);    
	glRotatef(-90.0, 1.0f, 0.0f, 0.0f);    
      }
    else if (strcmp (dir, "yz") == 0)
      {
	glRotatef(-90.0, 0.0f, 0.0f, 1.0f);    
	glRotatef(-90.0, 0.0f, 1.0f, 0.0f);    
      }
    else if (strcmp (dir, "zy") == 0)
      glRotatef(90.0, 0.0f, 1.0f, 0.0f);    


    glGetDoublev (GL_MODELVIEW_MATRIX, rotmat);

    glLoadIdentity();
    glMultMatrixd (lookatmat);
    glMultMatrixd (transmat);
    glMultMatrixd (rotmat);
    glMultMatrixd (centermat);
    glGetDoublev (GL_MODELVIEW_MATRIX, transformationmat);
  
    glPopMatrix();
  }

  void VisualScene :: MouseMove(int oldx, int oldy,
				int newx, int newy,
				char mode)
  {
    int deltax = newx - oldx;
    int deltay = newy - oldy;
  
    glPushMatrix();
    glLoadIdentity ();
  
    switch (mode)
      {
      case 'r':
	{	
	  glRotatef(float(deltax)/2, 0.0f, 1.0f, 0.0f);
	  glRotatef(float(deltay)/2, 1.0f, 0.0f, 0.0f);
	  glMultMatrixd (rotmat);
	  glGetDoublev (GL_MODELVIEW_MATRIX, rotmat);
	  break;
	}
      case 'm':
	{
	  GLdouble projmat[16], modelviewmat[16];
	  GLint viewport[4];
	  glGetDoublev (GL_PROJECTION_MATRIX, projmat);
	  glGetDoublev (GL_MODELVIEW_MATRIX, modelviewmat);
	  glGetIntegerv (GL_VIEWPORT, viewport);
	
	  // vorher pvz1/2 = 0
	  GLdouble pvx1 = 0, pvy1 = 0, pvz1 = 0.99; //  0.95;
	  GLdouble pvx2 = deltax, pvy2 = -deltay, pvz2 = 0.99; // 0.95;

	  GLdouble px1, py1, pz1;
	  GLdouble px2, py2, pz2;
	
	  gluUnProject (pvx1, pvy1, pvz1, 
			modelviewmat, projmat, viewport,
			&px1, &py1, &pz1);
	  gluUnProject (pvx2, pvy2, pvz2, 
			modelviewmat, projmat, viewport,
			&px2, &py2, &pz2);
	  /*
	    gluUnProject (oldx, oldy, 1, 
	    modelviewmat, projmat, viewport,
	    &px1, &py1, &pz1);
	    gluUnProject (newx, newy, 1, 
	    modelviewmat, projmat, viewport,
	    &px2, &py2, &pz2);
	  */

	  /*	
	    cout << "pv1 = " << pvx1 << ", " << pvy1 << ", " << pvz1 << endl;
	    cout << "p1 = " << px1 << ", " << py1 << ", " << pz1 << endl;
	  */

	  glTranslated (px2-px1, py2-py1, pz2-pz1);
	
	  glMultMatrixd (transmat);
	  glGetDoublev (GL_MODELVIEW_MATRIX, transmat);
	  break;
	}
      case 'z':
	{
	  // glTranslatef(0.0f, 0.0f, -dist);

	  // cout << "deltay = " << deltay << endl;
	  // cout << "float_bug = " << (float(deltay)/100) << endl;   gives wrong result with icc 9.0.021
	  glScaled (exp (double (-deltay)/100), 
		    exp (double (-deltay)/100), 
		    exp (double (-deltay)/100));
	  // glTranslatef(0.0f, 0.0f, dist);
	  glMultMatrixd (transmat);
	  glGetDoublev (GL_MODELVIEW_MATRIX, transmat);
	  break;
	}
      }

    glLoadIdentity();
    glMultMatrixd (lookatmat);
    glMultMatrixd (transmat);
    glMultMatrixd (rotmat);
    glMultMatrixd (centermat);
    glGetDoublev (GL_MODELVIEW_MATRIX, transformationmat);
  
    glPopMatrix();
  }


  void VisualScene :: LookAt (const Point<3> & cam, const Point<3> & obj,
			      const Point<3> & camup)
  {
    glPushMatrix();
    glLoadIdentity ();
    gluLookAt (cam(0), cam(1), cam(2), 
	       obj(0), obj(1), obj(2),
	       camup(0), camup(1), camup(2));
    glMultMatrixd (centermat);
    glGetDoublev (GL_MODELVIEW_MATRIX, transformationmat);
    glPopMatrix();
  }

  
  void VisualScene :: SetClippingPlane ()
  {
    if (vispar.clipping.enable)
      {
	Vec3d n = vispar.clipping.normal;
	n /= (n.Length()+1e-10);
	clipplane[0] = n.X();
	clipplane[1] = n.Y();
	clipplane[2] = n.Z();
	clipplane[3] = -(Vec3d(center) * n) + rad * vispar.clipping.dist;

	double clipplane2[4];
	clipplane2[0] = n.X();
	clipplane2[1] = n.Y();
	clipplane2[2] = n.Z();
	clipplane2[3] = -(Vec3d(center) * n) + 
	  rad * (vispar.clipping.dist + vispar.clipping.dist2);

	glClipPlane(GL_CLIP_PLANE0, clipplane2);
	glEnable(GL_CLIP_PLANE0);
      }
    else
      glDisable (GL_CLIP_PLANE0);
  }




  void VisualScene :: MouseDblClick (int /* px */, int /* py */)
  {
    ;
  }



  void VisualScene :: SetLight()
  {
    GLfloat vals[3];
    double lightamb = vispar.lightamb;
    vals[0] = vals[1] = vals[2] = lightamb;
    glLightfv(GL_LIGHT0, GL_AMBIENT, vals);

    double lightdiff = vispar.lightdiff;
    vals[0] = vals[1] = vals[2] = lightdiff;
    glLightfv(GL_LIGHT0, GL_DIFFUSE, vals);

    double lightspec = vispar.lightspec;
    vals[0] = vals[1] = vals[2] = lightspec;
    glLightfv(GL_LIGHT0, GL_SPECULAR, vals);

    glMaterialf (GL_FRONT_AND_BACK, GL_SHININESS, vispar.shininess);
    glLightModeli (GL_LIGHT_MODEL_LOCAL_VIEWER, vispar.locviewer);

    float mat_spec_col[] = { 1, 1, 1, 1 };
    glMaterialfv (GL_FRONT_AND_BACK, GL_SPECULAR, mat_spec_col);

    glEnable (GL_LIGHTING);
    glEnable (GL_LIGHT0);
  }




  void VisualScene :: SetOpenGlColor(double val, double valmin, double valmax,
				     int logscale)
  {
    double value;

    if (!logscale)
      value = (val - valmin) / (valmax - valmin);
    else
      {
	if (valmax <= 0) valmax = 1;
	if (valmin <= 0) valmin = 1e-4 * valmax;
	value = (log(fabs(val)) - log(valmin)) / (log(valmax) - log(valmin));
      }

    if (!invcolor)
      value = 1 - value;

    glTexCoord1f ( 0.998 * value + 0.001);
    // glTexCoord1f ( val ); 

    glTexCoord2f ( 0.998 * value + 0.001, 1.5);
    // glTexCoord1f ( value ); 

    if (value > 1) value = 1;
    if (value < 0) value = 0;

    value *= 4;

    static const double colp[][3] =
      {
	{ 1, 0, 0 },
	{ 1, 1, 0 },
	{ 0, 1, 0 },
	{ 0, 1, 1 },
	{ 0, 0, 1 },
	//	{ 1, 0, 1 },
	//	{ 1, 0, 0 },
      };
  
    int i = int(value);
    double r = value - i;

    GLdouble col[3];
    for (int j = 0; j < 3; j++)
      col[j] = (1-r) * colp[i][j] + r * colp[i+1][j];
  
    glColor3d (col[0], col[1], col[2]);
  }



  void VisualScene :: CreateTexture (int ncols, int linear, double alpha, int typ)
  {
    if (linear) ncols = 32;

    if (ntexcols != ncols) 
      {
	ntexcols = ncols;
      
	GLubyte colortexture[4*32];

	const double colp[][3] =
	  {
	    { 1, 0, 0 },
	    { 1, 1, 0 },
	    { 0, 1, 0 },
	    { 0, 1, 1 },
	    { 0, 0, 1 },
	  };
  
	for (int i = 0; i < ncols; i++)
	  {
	    double value = 4.0 * i / (ncols-1);

	    int iv = int(value);
	    double r = value - iv;

	    GLdouble col[3];

	    if(r > 1e-3)
	      for (int j = 0; j < 3; j++)
		col[j] = (1.-r) * colp[iv][j] + r * colp[iv+1][j];
	    else
	      for (int j = 0; j < 3; j++)
		col[j] = colp[iv][j];

	    colortexture[4*i] = GLubyte (255 * col[0]);
	    colortexture[4*i+1] = GLubyte (255 * col[1]);
	    colortexture[4*i+2] = GLubyte (255 * col[2]);
	    colortexture[4*i+3] = GLubyte(255*alpha);
	  }

	// glPixelStorei (GL_UNPACK_ALIGNMENT, 1);

     	glTexImage1D (GL_TEXTURE_1D, 0, 4, ncols, 0, GL_RGBA, GL_UNSIGNED_BYTE, colortexture);
	glTexImage2D (GL_TEXTURE_2D, 0, 4, ncols, 1, 0, GL_RGBA, GL_UNSIGNED_BYTE, colortexture);

	glTexEnvi (GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, typ);  // DECAL or MODULATE
	
	GLfloat bcol[] = { 1, 1, 1, 1.0 };
	glTexParameterfv (GL_TEXTURE_1D, GL_TEXTURE_BORDER_COLOR, bcol);
	glTexParameteri (GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);

	glTexParameterfv (GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, bcol);
	glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

	
	if (linear)
	  {
	    glTexParameteri (GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	    glTexParameteri (GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

	    glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	    glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	  }
	else
	  {
	    glTexParameteri (GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	    glTexParameteri (GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

	    glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	    glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	  }
      }
  }
  



  void VisualScene :: DrawColorBar (double minval, double maxval, int logscale, bool linear, string format, string unit)
  {
    if (!vispar.drawcolorbar) return;

    CreateTexture (GetVSSolution().numtexturecols, linear, 1, GL_DECAL);

    if (logscale && maxval <= 0) maxval = 1;
    if (logscale && minval <= 0) minval = 1e-4 * maxval;

    double minx = -1;
    double maxx = 1;
    double miny = 0.75;
    double maxy = 0.8;

    glDisable (GL_LIGHTING);
    glEnable (GL_COLOR_MATERIAL);
    glEnable (GL_TEXTURE_1D);
    glNormal3d (0, 0, 1);
    glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
    
    glDisable (GL_DEPTH_TEST);
    glBegin (GL_QUAD_STRIP);

    for (auto i : Range(50))
      {
        double x = minx + i*1.0/49*(maxx-minx);
	SetOpenGlColor (x, minx, maxx);
	glVertex3d (x, miny, -5);
	glVertex3d (x, maxy, -5);
      }
    glEnd();

    glDisable (GL_TEXTURE_1D);
    
    glEnable (GL_COLOR_MATERIAL);
    GLfloat textcol[3] = { GLfloat(1 - backcolor), 
                           GLfloat(1 - backcolor), 
                           GLfloat(1 - backcolor) };
    glColor3fv (textcol);
    
    glPushAttrib (GL_LIST_BIT);
    // glListBase (fontbase);

    constexpr size_t buf_size = 20;
    char buf[buf_size];
    GLint viewport[4];
    glGetIntegerv (GL_VIEWPORT, viewport);
    double char_width = 2.0*MyOpenGLTextWidth()/(viewport[3]);
    for (int i = 0; i <= 4; i++)
      {
	double val;
	if (logscale)
	  val = minval * pow (maxval / minval, i / 4.0);
	else
	  val = minval + i * (maxval-minval) / 4;

	snprintf (buf, buf_size, format.c_str(), val);
        auto n = strlen(buf);
	double x = minx + i * (maxx-minx) / 4;
        x -= 0.5*char_width * n; // center text
	glRasterPos3d (x, 0.7,-5);

	MyOpenGLText (buf);
      }

    if(unit != "")
        MyOpenGLText (unit.c_str());

    glPopAttrib ();
    glEnable (GL_DEPTH_TEST);
  }

  void VisualScene :: DrawTitle (string title)
  {
    if(title=="")
      return;
    glDisable (GL_LIGHTING);
    glDisable (GL_DEPTH_TEST);

    glEnable (GL_COLOR_MATERIAL);
    GLfloat textcol[3] = { GLfloat(1 - backcolor),
                           GLfloat(1 - backcolor),
                           GLfloat(1 - backcolor) };
    glColor3fv (textcol);

    glPushAttrib (GL_LIST_BIT);

    GLint viewport[4];
    glGetIntegerv (GL_VIEWPORT, viewport);
    double char_width = 2.0*MyOpenGLTextWidth()/(viewport[3]);
    double x = -0.5*char_width * title.size(); // center text
    glRasterPos3d (x, 0.82,-5);
    MyOpenGLText (title.c_str());
    glPopAttrib ();
    glEnable (GL_DEPTH_TEST);
  }


  void VisualScene :: DrawCoordinateCross ()
  {
    if (!vispar.drawcoordinatecross) return;

    glDisable (GL_DEPTH_TEST);
    glMatrixMode (GL_PROJECTION); 
    glPushMatrix();
    glLoadIdentity();

    glMatrixMode (GL_MODELVIEW); 
    glPushMatrix();
    glLoadIdentity();

    GLint viewport[4];
    glGetIntegerv (GL_VIEWPORT, viewport);

    glTranslatef (-1, -1, 0.0);
    glScalef (40.0 / viewport[2], 40.0 / viewport[3], 1);
    glTranslatef (2.0, 2.0, 0.0);
    glMultMatrixd (rotmat);

    glEnable (GL_COLOR_MATERIAL);
    glDisable (GL_LIGHTING);

    glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);

    GLfloat textcol[3] = { GLfloat(1 - backcolor),
			   GLfloat(1 - backcolor),
			   GLfloat(1 - backcolor) };
    glColor3fv (textcol);

    glLineWidth (1.0f);

    double len = 1;

    glBegin(GL_LINES);
    glVertex3d (0, 0, 0);
    glVertex3d (len, 0, 0);
    glVertex3d (0.0f, 0.0f, 0.0f);
    glVertex3d (0.0f, len, 0.0f);
    glVertex3d (0.0f, 0.0f, 0.0f);
    glVertex3d (0.0f, 0.0f, len);
    glEnd ();

    glPushAttrib (GL_LIST_BIT);
    // glListBase (fontbase);

    char buf[20];

    glRasterPos3d (len, 0.0f, 0.0f);
    snprintf (buf, size(buf), "x");
    // glCallLists (GLsizei(strlen (buf)), GL_UNSIGNED_BYTE, buf);
    MyOpenGLText (buf);
    glRasterPos3d (0.0f, len, 0.0f);
    snprintf (buf, size(buf), "y");
    // glCallLists (GLsizei(strlen (buf)), GL_UNSIGNED_BYTE, buf);
    MyOpenGLText (buf);
    glRasterPos3d (0.0f, 0.0f, len);
    snprintf (buf, size(buf), "z");
    // glCallLists (GLsizei(strlen (buf)), GL_UNSIGNED_BYTE, buf);
    MyOpenGLText (buf);

    glPopAttrib ();

    glEnable (GL_LIGHTING);

    glMatrixMode (GL_PROJECTION); 
    glPopMatrix();
    glMatrixMode (GL_MODELVIEW); 
    glPopMatrix();
    glEnable (GL_DEPTH_TEST);
  }


  void VisualScene :: DrawMarker()
  {
    static constexpr GLubyte cross[] = { 0xc6, 0xee, 0x7c, 0x38, 0x7c, 0xee, 0xc6 };

    if(!marker)
        return;

    glColor3d (0, 0, 1);

    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

    glDisable (GL_COLOR_MATERIAL);
    glDisable (GL_LIGHTING);
    glDisable (GL_CLIP_PLANE0);

    auto & p = *marker;
    glRasterPos3d (p[0], p[1], p[2]);
    glBitmap (7, 7, 3, 3, 0, 0, &cross[0]);
  }


  void VisualScene :: DrawNetgenLogo ()
  {
    if (!vispar.drawnetgenlogo) return;

    glDisable (GL_DEPTH_TEST);
    glMatrixMode (GL_PROJECTION); 
    glPushMatrix();
    glLoadIdentity();

    glMatrixMode (GL_MODELVIEW); 
    glPushMatrix();
    glLoadIdentity();

    GLint viewport[4];
    glGetIntegerv (GL_VIEWPORT, viewport);

    glTranslatef (1, -1, 0.0);
    glScalef (40.0 / viewport[2], 40.0 / viewport[3], 1);
    glTranslatef (-7.0, 2.0, 0.0);

    glDisable (GL_CLIP_PLANE0);
    glDisable (GL_LIGHTING);

    glEnable (GL_COLOR_MATERIAL);
    GLfloat textcol[3] = { GLfloat(1 - backcolor),
			   GLfloat(1 - backcolor),
			   GLfloat(1 - backcolor) };
    glColor3fv (textcol);
    glLineWidth (1.0f);

    glPushAttrib (GL_LIST_BIT);
    // glListBase (fontbase);

    char buf[] = "Netgen " PACKAGE_VERSION;

    glRasterPos3d (0.0f, 0.0f, 0.0f);
    // glCallLists (GLsizei(strlen (buf)), GL_UNSIGNED_BYTE, buf);
    MyOpenGLText (buf);

    glPopAttrib ();

    glEnable (GL_LIGHTING);
    glMatrixMode (GL_PROJECTION); 
    glPopMatrix();
    glMatrixMode (GL_MODELVIEW); 
    glPopMatrix();
    glEnable (GL_DEPTH_TEST);
  }


  void VisualSceneSurfaceMeshing::MouseMove(int oldx, int oldy,
                                            int newx, int newy,
                                            char mode)
  {
    double fac = 0.001;
    if(mode == 'M')
      {
        shiftx += fac * (newx - oldx);
        shifty += fac * (oldy - newy);
        return;
      }
    else if(mode == 'Z')
      {
        scalex *= (1 - fac * (newy - oldy));
        scaley *= (1 - fac * (newy - oldy));
        return;
      }
    
    VisualScene::MouseMove(oldx, oldy, newx, newy, mode);
  }

  std::vector<unsigned char> Snapshot( int w, int h )
  {
    // save current settings
    GLint viewport[4];
    glGetIntegerv (GL_VIEWPORT, viewport);

    glMatrixMode (GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();

    double pnear = 0.1;
    double pfar = 10;

    gluPerspective(20.0f, double(w) / h, pnear, pfar);

    glMatrixMode (GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glViewport(0,0,w,h);

    GLuint fb = 0;
    glGenFramebuffers(1, &fb);
    glBindFramebuffer(GL_FRAMEBUFFER, fb);

    // create, reserve and attach color and depth renderbuffer
    GLuint rbs[2];
    glGenRenderbuffers(2, rbs);
    glBindRenderbuffer(GL_RENDERBUFFER, rbs[0]);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_RGBA8, w, h);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, rbs[0]);

    glBindRenderbuffer(GL_RENDERBUFFER, rbs[1]);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, w, h);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, rbs[1]);

    // check if framebuffer status is complete
    if(int fbstatus; (fbstatus = glCheckFramebufferStatus(GL_FRAMEBUFFER)) != GL_FRAMEBUFFER_COMPLETE)
        cerr << "no frame buffer " << fbstatus << endl;

    visual_scene->DrawScene();
    glFinish();

    std::vector<unsigned char> buffer(w*h*3);
    glPixelStorei(GL_UNPACK_ALIGNMENT,1);
    glPixelStorei(GL_PACK_ALIGNMENT,1);
    glReadPixels (0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, &buffer[0]);

    glDeleteRenderbuffers(2, rbs);
    glDeleteFramebuffers(1, &fb);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);

    // restore previous settings
    glViewport(viewport[0], viewport[1], viewport[2], viewport[3]);
    glMatrixMode (GL_PROJECTION);
    glPopMatrix();
    glMatrixMode (GL_MODELVIEW);
    glPopMatrix();
    return buffer;
  }

  VisualSceneSurfaceMeshing :: VisualSceneSurfaceMeshing ()
    : VisualScene()
  {
    ;
  }

  VisualSceneSurfaceMeshing :: ~VisualSceneSurfaceMeshing ()
  {
    ;
  }

  void VisualSceneSurfaceMeshing :: DrawScene ()
  {
    // int i, j, k;
    if(!locpointsptr)
      return;
    auto& locpoints = *locpointsptr;
    auto& loclines = *loclinesptr;
    auto& plainpoints = *plainpointsptr;

    if (loclines.Size() != changeval)
      {
	center = Point<3>(0,0,-5);
	rad = 0.1;

	// CalcTransformationMatrices();
	changeval = loclines.Size();
      }

  glClearColor(backcolor, backcolor, backcolor, 1.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  SetLight();

  //  glEnable (GL_COLOR_MATERIAL);

  //  glDisable (GL_SHADING);
  //  glColor3f (0.0f, 1.0f, 1.0f);
  //  glLineWidth (1.0f);
  //  glShadeModel (GL_SMOOTH);

  //  glCallList (linelists.Get(1));

  //  SetLight();

  glPushMatrix();
  glMultMatrixd (transformationmat);

  glShadeModel (GL_SMOOTH);
  // glDisable (GL_COLOR_MATERIAL);
  glEnable (GL_COLOR_MATERIAL);
  glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);

  glEnable (GL_BLEND);
  glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  //  glEnable (GL_LIGHTING);

  double shine = vispar.shininess;
  double transp = vispar.transp;

  glMaterialf (GL_FRONT_AND_BACK, GL_SHININESS, shine);
  glLogicOp (GL_COPY);





  float mat_col[] = { 0.2, 0.2, 0.8, 1 };
  glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_col);

  glPolygonOffset (1, 1);
  glEnable (GL_POLYGON_OFFSET_FILL);

    float mat_colbl[] = { 0.8, 0.2, 0.2, 1 };
    float mat_cololdl[] = { 0.2, 0.8, 0.2, 1 };
    float mat_colnewl[] = { 0.8, 0.8, 0.2, 1 };


    glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
    glPolygonOffset (1, -1);
    glLineWidth (3);

    for (int i = 1; i <= loclines.Size(); i++)
      {
	if (i == 1)
	  {
	    glEnable (GL_POLYGON_OFFSET_FILL);
	    glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_colbl);
	  }
	else if (i <= oldnl)
	  glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_cololdl);
	else
	  glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_colnewl);

	int pi1 = loclines.Get(i).I1();
	int pi2 = loclines.Get(i).I2();

	if (pi1 >= 1 && pi2 >= 1)
	  {
	    Point3d p1 = locpoints.Get(pi1);
	    Point3d p2 = locpoints.Get(pi2);

	    glBegin (GL_LINES);
	    glVertex3f (p1.X(), p1.Y(), p1.Z());
	    glVertex3f (p2.X(), p2.Y(), p2.Z());
	    glEnd();
	  }

	glDisable (GL_POLYGON_OFFSET_FILL);
      }


    glLineWidth (1);


    glPointSize (5);
    float mat_colp[] = { 1, 0, 0, 1 };
    glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_colp);
    glBegin (GL_POINTS);
    for (int i = 1; i <= locpoints.Size(); i++)
      {
	Point3d p = locpoints.Get(i);
	glVertex3f (p.X(), p.Y(), p.Z());
      }
    glEnd();


    glPopMatrix();


    // float mat_colp[] = { 1, 0, 0, 1 };

    float mat_col2d1[] = { 1, 0.5, 0.5, 1 };
    float mat_col2d[] = { 1, 1, 1, 1 };
    glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_col2d);

    glBegin (GL_LINES);
    for (int i = 1; i <= loclines.Size(); i++)
      {
	glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_col2d);
	if (i == 1)
	  glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_col2d1);

	int pi1 = loclines.Get(i).I1();
	int pi2 = loclines.Get(i).I2();

	if (pi1 >= 1 && pi2 >= 1)
	  {
	    const auto& p1 = plainpoints.Get(pi1);
	    const auto& p2 = plainpoints.Get(pi2);

	    glBegin (GL_LINES);
	    glVertex3f (scalex * p1[0] + shiftx, scaley * p1[1] + shifty, -5);
	    glVertex3f (scalex * p2[0] + shiftx, scaley * p2[1] + shifty, -5);
	    glEnd();
	  }
      }
    glEnd ();


    glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_colp);
    glBegin (GL_POINTS);
    for (int i = 1; i <= plainpoints.Size(); i++)
      {
	const auto& p = plainpoints.Get(i);
	glVertex3f (scalex * p[0] + shiftx, scaley * p[1] + shifty, -5);
      }
    glEnd();






  glDisable (GL_POLYGON_OFFSET_FILL);

  glPopMatrix();
  DrawCoordinateCross ();
  DrawNetgenLogo ();
  glFinish();

  }


  void VisualSceneSurfaceMeshing :: BuildScene (int zoomall)
  {
  }

  VisualSceneSurfaceMeshing vssurfacemeshing;

  void Impl_Render (bool blocking)
  {
    if (blocking && multithread.running)
      {
        multithread.redraw = 2;
        while (multithread.redraw == 2) ;
      }
    else
      multithread.redraw = 1;
  }

  void Impl_UpdateVisSurfaceMeshData(int oldnl,
            shared_ptr<NgArray<Point<3>>> locpointsptr,
            shared_ptr<NgArray<INDEX_2>> loclinesptr,
            shared_ptr<NgArray<Point<2>>> plainpointsptr)
  {
      vssurfacemeshing.oldnl = oldnl;
      if(locpointsptr) vssurfacemeshing.locpointsptr = locpointsptr;
      if(loclinesptr) vssurfacemeshing.loclinesptr = loclinesptr;
      if(plainpointsptr) vssurfacemeshing.plainpointsptr = plainpointsptr;
  }

  static bool set_function_pointers = []()
  {
      Ptr_Render = Impl_Render;
      Ptr_UpdateVisSurfaceMeshData = Impl_UpdateVisSurfaceMeshData;
      return true;
  }();



#ifdef PARALLELGL
  void VisualScene :: InitParallelGL ()
  {
    static int init = 0;

    if (!init)
      {
	init = 1;

	if (id == 0)
	  {
	    string displname;
	    
	    Display * dpy = glXGetCurrentDisplay();
	    GLXDrawable drawable = glXGetCurrentDrawable();
	    GLXContext ctx = glXGetCurrentContext();
	    GLXContextID xid = glXGetContextIDEXT (ctx);
	    
	    displname = XDisplayName (0);

	    if( glXIsDirect ( dpy, ctx ) )
	      cout << "WARNING: direct rendering enabled; this might break mpi-parallel netgen (especially if X-forwarding is used! (to disable, change -indirect to true in ng/drawing.tcl)" << endl;
	      
	    /*
	    cout << "Init Parallel GL" << endl;
	    cout << "DisplayName = " << displname << endl;
	    cout << "current display = " << dpy << endl;
	    cout << "current drawable = " << drawable << endl;                  
	    cout << "current context = " << ctx << endl;                  
	    
	    cout << "contextid = " << xid << endl;
	    cout << "isdirect = " << glXIsDirect ( dpy, ctx ) << endl;                  
	    cout << "extensionstring = " << glXQueryExtensionsString( dpy, 0 ) << endl;
	    */

	    MyMPI_SendCmd ("redraw");
	    MyMPI_SendCmd ("init");
		
	    for (int dest = 1; dest < ntasks; dest++)
	      {
		MyMPI_Send (displname, dest, MPI_TAG_VIS);
		MyMPI_Send (int (drawable), dest, MPI_TAG_VIS);
		MyMPI_Send (int (xid), dest, MPI_TAG_VIS);
	      } 
	  }
      }
  }


  void VisualScene :: Broadcast ()
  {
    if (ntasks == 1) return;

    if (id == 0)
      {
	/*
	for (int dest = 1; dest < ntasks; dest++)
	  {
	    MyMPI_Send ("redraw", dest, MPI_TAG_CMD);
	    MyMPI_Send ("broadcast", dest, MPI_TAG_VIS);
	  }
	*/

	MyMPI_SendCmd ("redraw");
	MyMPI_SendCmd ("broadcast");
      }

    MyMPI_Bcast (selface);

    netgen::GetVSSolution().Broadcast ();
  }
#endif 

}
