
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

  MeshingParameters CreateMPfromKwargs(py::kwargs kwargs);
} // namespace netgen
