#ifndef FILE_VISPAR
#define FILE_VISPAR

namespace netgen
{

class VisualizationParameters
{
public:
  double lightamb;
  double lightdiff;
  double lightspec;
  double shininess;
  double transp;
  int locviewer;
  char selectvisual[20];
  int showstltrias;
  
  /*
  Vec3d clipnormal;
  double clipdist;
  int clipenable;
  int clipplanetimestamp;
  */
  class Clipping
  {
  public:
    Vec3d normal;
    double dist;
    double dist2;
    int enable;
    int timestamp;
    bool operator== (Clipping & clip2)
    {
      return 
	(normal == clip2.normal) && 
	(dist == clip2.dist) && 
	// (dist2 == clip2.dist2) && 
	(enable == clip2.enable);
    }
  };
  Clipping clipping;

  int colormeshsize;

  int drawfilledtrigs;
  int drawbadels;
  int drawoutline;
  int drawedges;
  int subdivisions;

  int drawprisms;
  int drawpyramids;
  int drawhexes;
  double shrink;
  int drawidentified;
  int drawpointnumbers;
  int drawedgenumbers;
  int drawfacenumbers;
  int drawelementnumbers;
  int drawsurfaceelementnumbers;
  int drawsegmentnumbers;
  int drawdomainsurf;
  int drawtets;
  int drawtetsdomain;

  int clipdomain;
  int donotclipdomain;

  int drawededges;
  int drawedpoints;
  int drawedpointnrs;
  int drawedtangents;
  int drawededgenrs;
  int drawmetispartition;

  int drawcurveproj;
  int drawcurveprojedge;
  

  PointIndex centerpoint;
  int drawelement;

  // stl:
  int stlshowtrias;
  int stlshowfilledtrias;
  int stlshowedges;
  int stlshowmarktrias;
  int stlshowactivechart;
  int stlchartnumber;
  int stlchartnumberoffset;

  // occ:
  int occshowvolumenr;
  bool occshowsurfaces;
  bool occshowedges;
  bool occvisproblemfaces;
  bool occzoomtohighlightedentity;
  double occdeflection;

  // ACIS

  bool ACISshowfaces;
  bool ACISshowedges;
  int ACISshowsolidnr;
  int ACISshowsolidnr2;

  bool whitebackground;
  int stereo;
  bool usedispllists;
  bool drawcoordinatecross;
  bool drawcolorbar;
  bool drawnetgenlogo;

  bool use_center_coords;
  double centerx,centery,centerz;

  bool drawspecpoint;
  double specpointx,specpointy,specpointz;

  
public:
  VisualizationParameters();
};
NGGUI_API extern VisualizationParameters vispar;
}

#endif
