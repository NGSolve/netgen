
#include <pybind11/pybind11.h>
#include "meshing.hpp"

namespace netgen
{
  // TODO: Clarify a lot of these parameters
  static string meshingparameter_description = R"delimiter(
Meshing Parameters
-------------------

maxh: float = 1e10
  Global upper bound for mesh size.

grading: float = 0.3
  Mesh grading how fast the local mesh size can change.

meshsizefilename: str = None
  Load meshsize from file. Can set local mesh size for points
  and along edges. File must have the format:

    nr_points
    x1, y1, z1, meshsize
    x2, y2, z2, meshsize
    ...
    xn, yn, zn, meshsize

    nr_edges
    x11, y11, z11, x12, y12, z12, meshsize
    ...
    xn1, yn1, zn1, xn2, yn2, zn2, meshsize

segmentsperedge: float = 1.
  Minimal number of segments per edge.

quad: bool = False
  Quad-dominated surface meshing.

blockfill: bool = True
  Do fast blockfilling.

filldist: float = 0.1
  Block fill up to distance

delaunay: bool = True
  Use delaunay meshing.

Optimization Parameters
-----------------------

optimize3d: str = "cmdmustm"
  3d optimization strategy:
    m .. move nodes
    M .. move nodes, cheap functional
    s .. swap faces
    c .. combine elements
    d .. divide elements
    p .. plot, no pause
    P .. plot, Pause
    h .. Histogramm, no pause
    H .. Histogramm, pause

optsteps3d: int = 3
  Number of 3d optimization steps.

optimize2d: str = "smsmsmSmSmSm"
  2d optimization strategy:
    s .. swap, opt 6 lines/node
    S .. swap, optimal elements
    m .. move nodes
    p .. plot, no pause
    P .. plot, pause
    c .. combine

optsteps2d: int = 3
  Number of 2d optimization steps.

elsizeweight: float = 0.2
  Weight of element size w.r.t. element shape in optimization.

)delimiter";

  inline void CreateMPfromKwargs(MeshingParameters& mp, py::kwargs kwargs)
  {
    if(kwargs.contains("optimize3d"))
      mp.optimize3d = py::cast<string>(kwargs["optimize3d"]);
    if(kwargs.contains("optsteps3d"))
      mp.optsteps3d = py::cast<int>(kwargs["optsteps3d"]);
    if(kwargs.contains("optimize2d"))
      mp.optimize2d = py::cast<string>(kwargs["optimize2d"]);
    if(kwargs.contains("optsteps2d"))
      mp.optsteps2d = py::cast<int>(kwargs["optsteps2d"]);
    if(kwargs.contains("opterrpow"))
      mp.opterrpow = py::cast<double>(kwargs["opterrpow"]);
    if(kwargs.contains("blockfill"))
      mp.blockfill = py::cast<bool>(kwargs["blockfill"]);
    if(kwargs.contains("filldist"))
      mp.filldist = py::cast<double>(kwargs["filldist"]);
    if(kwargs.contains("safety"))
      mp.safety = py::cast<double>(kwargs["safety"]);
    if(kwargs.contains("relinnersafety"))
      mp.relinnersafety = py::cast<double>(kwargs["relinnersafety"]);
    if(kwargs.contains("uselocalh"))
      mp.uselocalh = py::cast<bool>(kwargs["uselocalh"]);
    if(kwargs.contains("grading"))
      mp.grading = py::cast<double>(kwargs["grading"]);
    if(kwargs.contains("delaunay"))
      mp.delaunay = py::cast<bool>(kwargs["delaunay"]);
    if(kwargs.contains("maxh"))
      mp.maxh = py::cast<double>(kwargs["maxh"]);
    if(kwargs.contains("minh"))
      mp.minh = py::cast<double>(kwargs["minh"]);
    if(kwargs.contains("meshsizefilename"))
      mp.meshsizefilename = py::cast<string>(kwargs["meshsizefilename"]);
    if(kwargs.contains("startinsurface"))
      mp.startinsurface = py::cast<bool>(kwargs["startinsurface"]);
    if(kwargs.contains("checkoverlap"))
      mp.checkoverlap = py::cast<bool>(kwargs["checkoverlap"]);
    if(kwargs.contains("checkoverlappingboundary"))
      mp.checkoverlappingboundary = py::cast<bool>(kwargs["checkoverlappingboundary"]);
    if(kwargs.contains("checkchartboundary"))
      mp.checkchartboundary = py::cast<bool>(kwargs["checkchartboundary"]);
    if(kwargs.contains("curvaturesafety"))
      mp.curvaturesafety = py::cast<double>(kwargs["curvaturesafety"]);
    if(kwargs.contains("segmentsperedge"))
      mp.segmentsperedge = py::cast<double>(kwargs["segmentsperedge"]);
    if(kwargs.contains("parthread"))
      mp.parthread = py::cast<bool>(kwargs["parthread"]);
    if(kwargs.contains("elsizeweight"))
      mp.elsizeweight = py::cast<double>(kwargs["elsizeweight"]);
    if(kwargs.contains("perfstepsstart"))
      mp.perfstepsstart = py::cast<int>(kwargs["perfstepsstart"]);
    if(kwargs.contains("perfstepsend"))
      mp.perfstepsend = py::cast<int>(kwargs["perfstepsend"]);
    if(kwargs.contains("giveuptol2d"))
      mp.giveuptol2d = py::cast<int>(kwargs["giveuptol2d"]);
    if(kwargs.contains("giveuptol"))
      mp.giveuptol = py::cast<int>(kwargs["giveuptol"]);
    if(kwargs.contains("maxoutersteps"))
      mp.maxoutersteps = py::cast<int>(kwargs["maxoutersteps"]);
    if(kwargs.contains("starshapeclass"))
      mp.starshapeclass = py::cast<int>(kwargs["starshapeclass"]);
    if(kwargs.contains("baseelnp"))
      mp.baseelnp = py::cast<int>(kwargs["baseelnp"]);
    if(kwargs.contains("sloppy"))
      mp.sloppy = py::cast<int>(kwargs["sloppy"]);
    if(kwargs.contains("badellimit"))
      mp.badellimit = py::cast<double>(kwargs["badellimit"]);
    if(kwargs.contains("check_impossible"))
      mp.check_impossible = py::cast<bool>(kwargs["check_impossible"]);
    if(kwargs.contains("only3D_domain_nr"))
      mp.only3D_domain_nr = py::cast<int>(kwargs["only3D_domain_nr"]);
    if(kwargs.contains("secondorder"))
      mp.secondorder = py::cast<bool>(kwargs["secondorder"]);
    if(kwargs.contains("elementorder"))
      mp.elementorder = py::cast<int>(kwargs["elementorder"]);
    if(kwargs.contains("quad"))
      mp.quad = py::cast<bool>(kwargs["quad"]);
    if(kwargs.contains("try_hexes"))
      mp.try_hexes = py::cast<bool>(kwargs["try_hexes"]);
    if(kwargs.contains("inverttets"))
      mp.inverttets = py::cast<bool>(kwargs["inverttets"]);
    if(kwargs.contains("inverttrigs"))
      mp.inverttrigs = py::cast<bool>(kwargs["inverttrigs"]);
    if(kwargs.contains("autozrefine"))
      mp.autozrefine = py::cast<bool>(kwargs["autozrefine"]);
  }
} // namespace netgen
