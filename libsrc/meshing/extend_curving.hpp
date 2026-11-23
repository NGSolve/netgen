#ifndef NETGEN_EXTEND_CURVING_HPP
#define NETGEN_EXTEND_CURVING_HPP

/**************************************************************************/
/* File:   extend_curving.hpp                                             */
/* Author: Christopher Lackner                                            */
/* Date:   Oct 2026                                                       */
/**************************************************************************/

#include "meshclass.hpp"

namespace netgen
{

  struct ExtendCurvingParameters
  {
    double quality = 0.2;        // target min det J_c / det J_s, bending and straightening
    double eps = 0.05;           // hard bound, below this the boundary curving is reduced
    int max_depth = 4;           // element layers the bending may go into the interior
    int max_sweeps = 20;         // Gauss-Seidel sweeps per depth / continuation step
    bool reduce = true;          // allow reducing the boundary curving
    int max_release_passes = 4;
    double straighten_tol = 0.;  // boundary entities with |c| < tol*h are straightened
    int bezier_depth = 3;        // subdivision depth of the validity check
    double beta = 5.;            // stiffness growth per layer
    bool move_vertices = false;
    int vertex_layers = 1;       // vertex layers around the bad elements that move
    int vertex_sweeps = 10;
  };

  struct ExtendCurvingInfo
  {
    size_t nbad = 0;             // elements below quality after Curve
    size_t ncurved_before = 0;
    size_t ncurved_after = 0;
    size_t nreleased = 0;        // boundary entities with reduced curving
    size_t nremaining = 0;       // elements still below eps
    size_t nundecided = 0;       // elements the Bezier check could not decide
    double geometry_error = 0;   // max distance of the released boundary from the geometry
    Array<size_t> depth_histogram; // bad elements handled at each depth
    size_t nmoved = 0;           // inner vertices moved
    double max_vertex_move = 0;  // largest vertex move relative to the local element size
    Array<size_t> causes;
  };

  DLL_HEADER ExtendCurvingInfo
  ExtendCurving (Mesh & mesh, const ExtendCurvingParameters & params);

  DLL_HEADER void CurvedElementQuality (const Mesh & mesh,
                                        Array<double, ElementIndex> & lower,
                                        Array<double, ElementIndex> & upper,
                                        int bezier_depth = 6);
}

#endif
