#ifndef NETGEN_MESHING_PYTHON_MESH_HPP
#define NETGEN_MESHING_PYTHON_MESH_HPP

#include <core/python_ngcore.hpp>

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

quad_dominated: bool = False
  Quad-dominated surface meshing.

blockfill: bool = True
  Do fast blockfilling.

filldist: float = 0.1
  Block fill up to distance

delaunay: bool = True
  Use delaunay meshing.

delaunay2d : bool = True
  Use delaunay meshing for 2d geometries.

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

optimize2d: str = "smcmSmcmSmcm"
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

inline void CreateMPfromKwargs(MeshingParameters& mp, py::kwargs kwargs, bool throw_if_not_all_parsed=true)
  {
    if(kwargs.contains("optimize3d"))
      mp.optimize3d = py::cast<string>(kwargs.attr("pop")("optimize3d"));
    if(kwargs.contains("optsteps3d"))
      mp.optsteps3d = py::cast<int>(kwargs.attr("pop")("optsteps3d"));
    if(kwargs.contains("optimize2d"))
      mp.optimize2d = py::cast<string>(kwargs.attr("pop")("optimize2d"));
    if(kwargs.contains("optsteps2d"))
      mp.optsteps2d = py::cast<int>(kwargs.attr("pop")("optsteps2d"));
    if(kwargs.contains("opterrpow"))
      mp.opterrpow = py::cast<double>(kwargs.attr("pop")("opterrpow"));
    if(kwargs.contains("blockfill"))
      mp.blockfill = py::cast<bool>(kwargs.attr("pop")("blockfill"));
    if(kwargs.contains("filldist"))
      mp.filldist = py::cast<double>(kwargs.attr("pop")("filldist"));
    if(kwargs.contains("safety"))
      mp.safety = py::cast<double>(kwargs.attr("pop")("safety"));
    if(kwargs.contains("relinnersafety"))
      mp.relinnersafety = py::cast<double>(kwargs.attr("pop")("relinnersafety"));
    if(kwargs.contains("uselocalh"))
      mp.uselocalh = py::cast<bool>(kwargs.attr("pop")("uselocalh"));
    if(kwargs.contains("grading"))
      mp.grading = py::cast<double>(kwargs.attr("pop")("grading"));
    if(kwargs.contains("delaunay"))
      mp.delaunay = py::cast<bool>(kwargs.attr("pop")("delaunay"));
    if(kwargs.contains("delaunay2d"))
      mp.delaunay2d = py::cast<bool>(kwargs.attr("pop")("delaunay2d"));
    if(kwargs.contains("maxh"))
      mp.maxh = py::cast<double>(kwargs.attr("pop")("maxh"));
    if(kwargs.contains("minh"))
      mp.minh = py::cast<double>(kwargs.attr("pop")("minh"));
    if(kwargs.contains("meshsizefilename"))
      mp.meshsizefilename = py::cast<string>(kwargs.attr("pop")("meshsizefilename"));
    if(kwargs.contains("startinsurface"))
      mp.startinsurface = py::cast<bool>(kwargs.attr("pop")("startinsurface"));
    if(kwargs.contains("checkoverlap"))
      mp.checkoverlap = py::cast<bool>(kwargs.attr("pop")("checkoverlap"));
    if(kwargs.contains("checkoverlappingboundary"))
      mp.checkoverlappingboundary = py::cast<bool>(kwargs.attr("pop")("checkoverlappingboundary"));
    if(kwargs.contains("checkchartboundary"))
      mp.checkchartboundary = py::cast<bool>(kwargs.attr("pop")("checkchartboundary"));
    if(kwargs.contains("curvaturesafety"))
      mp.curvaturesafety = py::cast<double>(kwargs.attr("pop")("curvaturesafety"));
    if(kwargs.contains("segmentsperedge"))
      mp.segmentsperedge = py::cast<double>(kwargs.attr("pop")("segmentsperedge"));
    if(kwargs.contains("parthread"))
      mp.parthread = py::cast<bool>(kwargs.attr("pop")("parthread"));
    if(kwargs.contains("elsizeweight"))
      mp.elsizeweight = py::cast<double>(kwargs.attr("pop")("elsizeweight"));
    if(kwargs.contains("perfstepsstart"))
      mp.perfstepsstart = py::cast<int>(kwargs.attr("pop")("perfstepsstart"));
    if(kwargs.contains("perfstepsend"))
      mp.perfstepsend = py::cast<int>(kwargs.attr("pop")("perfstepsend"));
    if(kwargs.contains("giveuptol2d"))
      mp.giveuptol2d = py::cast<int>(kwargs.attr("pop")("giveuptol2d"));
    if(kwargs.contains("giveuptol"))
      mp.giveuptol = py::cast<int>(kwargs.attr("pop")("giveuptol"));
    if(kwargs.contains("giveuptolopenquads"))
      mp.giveuptolopenquads = py::cast<int>(kwargs.attr("pop")("giveuptolopenquads"));
    if(kwargs.contains("maxoutersteps"))
      mp.maxoutersteps = py::cast<int>(kwargs.attr("pop")("maxoutersteps"));
    if(kwargs.contains("starshapeclass"))
      mp.starshapeclass = py::cast<int>(kwargs.attr("pop")("starshapeclass"));
    if(kwargs.contains("baseelnp"))
      mp.baseelnp = py::cast<int>(kwargs.attr("pop")("baseelnp"));
    if(kwargs.contains("sloppy"))
      mp.sloppy = py::cast<int>(kwargs.attr("pop")("sloppy"));
    if(kwargs.contains("badellimit"))
      mp.badellimit = py::cast<double>(kwargs.attr("pop")("badellimit"));
    if(kwargs.contains("check_impossible"))
      mp.check_impossible = py::cast<bool>(kwargs.attr("pop")("check_impossible"));
    if(kwargs.contains("only3D_domain_nr"))
      mp.only3D_domain_nr = py::cast<int>(kwargs.attr("pop")("only3D_domain_nr"));
    if(kwargs.contains("secondorder"))
      mp.secondorder = py::cast<bool>(kwargs.attr("pop")("secondorder"));
    if(kwargs.contains("elementorder"))
      mp.elementorder = py::cast<int>(kwargs.attr("pop")("elementorder"));
    if(kwargs.contains("quad"))
      {
        cout << "WARNING: Meshing parameter 'quad' is deprecated, use 'quad_dominated' instead!" << endl;
        mp.quad = py::cast<bool>(kwargs.attr("pop")("quad"));
      }
    if(kwargs.contains("quad_dominated"))
      mp.quad = py::cast<bool>(kwargs.attr("pop")("quad_dominated"));
    if(kwargs.contains("try_hexes"))
      mp.try_hexes = py::cast<bool>(kwargs.attr("pop")("try_hexes"));
    if(kwargs.contains("inverttets"))
      mp.inverttets = py::cast<bool>(kwargs.attr("pop")("inverttets"));
    if(kwargs.contains("inverttrigs"))
      mp.inverttrigs = py::cast<bool>(kwargs.attr("pop")("inverttrigs"));
    if(kwargs.contains("autozrefine"))
      mp.autozrefine = py::cast<bool>(kwargs.attr("pop")("autozrefine"));
    if(kwargs.contains("parallel_meshing"))
      mp.parallel_meshing = py::cast<bool>(kwargs.attr("pop")("parallel_meshing"));
    if(kwargs.contains("nthreads"))
      mp.nthreads = py::cast<int>(kwargs.attr("pop")("nthreads"));
    if(kwargs.contains("closeedgefac"))
      mp.closeedgefac = py::cast<optional<double>>(kwargs.attr("pop")("closeedgefac"));

    if(kwargs.contains("boundary_layers"))
    {
      auto layers = py::list(kwargs.attr("pop")("boundary_layers"));
      for(auto layer : layers)
        mp.boundary_layers.Append(py::cast<BoundaryLayerParameters>(layer));
    }

    if(kwargs.size())
    {
      if(throw_if_not_all_parsed)
        throw Exception(string("Not all kwargs given to GenerateMesh could be parsed:") + string(py::str(kwargs)));
      mp.geometrySpecificParameters = CreateFlagsFromKwArgs(kwargs);
    }
  }
} // namespace netgen

#endif // NETGEN_MESHING_PYTHON_MESH_HPP
 
