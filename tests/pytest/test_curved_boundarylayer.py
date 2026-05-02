"""
Tests for curved boundary layers (TDD).

These tests are expected to FAIL until the curved BL implementation is done.
Each test corresponds to a phase in doc/curved_boundarylayer_plan.md.

Test execution order matches implementation phases:
  T1 - BL point map storage and serialization
  T2 - Volume convergence on sphere (single layer)
  T3 - Volume convergence on cylinder
  T4 - Multi-layer volume convergence on sphere
  T5 - Positive Jacobian on thick BL
  T6 - Press-through prevention on extreme geometry
"""

import math
import pytest

occ = pytest.importorskip("netgen.occ")
from netgen.meshing import BoundaryLayerParameters, Mesh


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------

def make_sphere_with_bl(r=1.0, bl_thickness=[0.05], maxh=0.4, outside=False):
    """Sphere with inward BL. Returns (geo, bl_params, maxh)."""
    sphere = occ.Sphere(occ.Pnt(0, 0, 0), r)
    geo = occ.OCCGeometry(sphere)
    bl = [BoundaryLayerParameters(".*", bl_thickness, "layer",
                                  outside=outside, disable_curving=False)]
    return geo, bl, maxh


def make_cylinder_with_bl(r=0.5, h=2.0, bl_thickness=[0.02], maxh=0.4):
    """Cylinder with BL on all faces. Returns (geo, bl_params, maxh)."""
    cyl = occ.Cylinder(occ.Pnt(0, 0, 0), occ.Z, r=r, h=h)
    cyl.faces.Max(occ.Z).name = "top"
    cyl.faces.Min(occ.Z).name = "bottom"
    geo = occ.OCCGeometry(cyl)
    bl = [BoundaryLayerParameters(".*", bl_thickness, "layer",
                                  outside=False, disable_curving=False)]
    return geo, bl, maxh


def sphere_volume(r):
    return 4.0 / 3.0 * math.pi * r**3


def sphere_shell_volume(r_outer, r_inner):
    return sphere_volume(r_outer) - sphere_volume(r_inner)


def cylinder_volume(r, h):
    return math.pi * r**2 * h


# ---------------------------------------------------------------------------
# T1: BL point map storage and serialization
# ---------------------------------------------------------------------------

class TestBLPointMap:
    """Phase 1: boundary_layer_point_map is populated and survives save/load."""

    def test_bl_point_map_exists(self):
        """After BL generation with disable_curving=False, the mesh must have
        a non-empty BL point map."""
        geo, bl, maxh = make_sphere_with_bl()
        mesh = geo.GenerateMesh(maxh=maxh, boundary_layers=bl)
        bl_map = mesh.GetBoundaryLayerPointMap()
        valid_count = sum(1 for info in bl_map if info.IsValid())
        assert valid_count > 0, "BL point map has no valid entries after BL generation"

    def test_bl_point_map_save_load(self, tmp_path):
        """BL point map must survive save/load cycle."""
        geo, bl, maxh = make_sphere_with_bl()
        mesh = geo.GenerateMesh(maxh=maxh, boundary_layers=bl)
        bl_map_before = mesh.GetBoundaryLayerPointMap()
        entries_before = {i: (info.base_pi, info.height)
                          for i, info in enumerate(bl_map_before, start=1) if info.IsValid()}

        path = str(tmp_path / "sphere_bl.vol")
        mesh.Save(path)
        mesh2 = Mesh()
        mesh2.Load(path)
        bl_map_after = mesh2.GetBoundaryLayerPointMap()
        entries_after = {i: (info.base_pi, info.height)
                         for i, info in enumerate(bl_map_after, start=1) if info.IsValid()}

        assert len(entries_after) == len(entries_before), \
            "BL point map size changed after save/load"
        for key in entries_before:
            assert key in entries_after, f"Point {key} missing after load"
            assert entries_after[key][0] == entries_before[key][0], \
                f"base_pi mismatch at {key}"
            assert entries_after[key][1] == pytest.approx(entries_before[key][1]), \
                f"height mismatch at {key}"


# ---------------------------------------------------------------------------
# T2: Volume convergence on sphere (single layer)
#
# The key test is the INNER CORE volume: 4/3 * pi * (r-h)^3.
# Without BL curving the inner surface is flat → inner core volume is wrong.
# The total sphere volume passes trivially because the outer boundary is
# curved to CAD geometry by standard curving.
# ---------------------------------------------------------------------------

class TestSphereVolume:
    """Phase 2+3: curved BL on sphere must give correct inner core volume."""

    @pytest.mark.parametrize("order", [2, 3, 4])
    def test_sphere_bl_core_volume(self, order):
        """Inner core (non-BL domain) of sphere must have correct volume."""
        ngs = pytest.importorskip("ngsolve")
        r = 1.0
        h = 0.05
        geo, bl, maxh = make_sphere_with_bl(r=r, bl_thickness=[h], maxh=0.4)
        ngmesh = geo.GenerateMesh(maxh=maxh, boundary_layers=bl)
        mesh = ngs.Mesh(ngmesh)
        mesh.Curve(order)
        core_vol = ngs.Integrate(1.0, mesh, definedon=mesh.Materials("default"))
        expected = sphere_volume(r - h)
        assert core_vol == pytest.approx(expected, rel=0.01), \
            f"Inner core volume {core_vol} != expected {expected} at order {order}"

    def test_sphere_bl_shell_volume(self):
        """The BL shell must have correct volume."""
        ngs = pytest.importorskip("ngsolve")
        r = 1.0
        h = 0.05
        geo, bl, maxh = make_sphere_with_bl(r=r, bl_thickness=[h], maxh=0.4)
        ngmesh = geo.GenerateMesh(maxh=maxh, boundary_layers=bl)
        mesh = ngs.Mesh(ngmesh)
        mesh.Curve(4)
        layer_vol = ngs.Integrate(1.0, mesh, definedon=mesh.Materials("layer"))
        expected_shell = sphere_shell_volume(r, r - h)
        assert layer_vol == pytest.approx(expected_shell, rel=0.05), \
            f"Shell volume {layer_vol} != expected {expected_shell}"


# ---------------------------------------------------------------------------
# T3: Volume convergence on cylinder
# ---------------------------------------------------------------------------

class TestCylinderVolume:
    """Phase 2+3: curved BL on cylinder must give correct inner core volume."""

    @pytest.mark.parametrize("order", [2, 3])
    def test_cylinder_bl_core_volume(self, order):
        """Cylinder with BL: inner core volume must match analytic value."""
        ngs = pytest.importorskip("ngsolve")
        r = 0.5
        h_cyl = 2.0
        bl_h = 0.02
        geo, bl, maxh = make_cylinder_with_bl(r=r, h=h_cyl, bl_thickness=[bl_h],
                                               maxh=0.4)
        ngmesh = geo.GenerateMesh(maxh=maxh, boundary_layers=bl)
        mesh = ngs.Mesh(ngmesh)
        mesh.Curve(order)
        core_vol = ngs.Integrate(1.0, mesh, definedon=mesh.Materials("default"))
        # Inner core: cylinder with r-bl_h radius and h-2*bl_h height (BL on all faces)
        expected = cylinder_volume(r - bl_h, h_cyl - 2*bl_h)
        assert core_vol == pytest.approx(expected, rel=0.01), \
            f"Cylinder core volume {core_vol} != expected {expected} at order {order}"


# ---------------------------------------------------------------------------
# T4: Multi-layer volume convergence
# ---------------------------------------------------------------------------

class TestMultiLayer:
    """Phase 4: multiple BL layers must all be curved correctly."""

    @pytest.mark.parametrize("nlayers", [2, 3])
    def test_sphere_multilayer_core_volume(self, nlayers):
        """Sphere with multiple BL layers: inner core volume must match."""
        ngs = pytest.importorskip("ngsolve")
        r = 1.0
        layer_h = 0.02
        thicknesses = [layer_h] * nlayers
        total_h = sum(thicknesses)
        geo, bl, maxh = make_sphere_with_bl(r=r, bl_thickness=thicknesses, maxh=0.4)
        ngmesh = geo.GenerateMesh(maxh=maxh, boundary_layers=bl)
        mesh = ngs.Mesh(ngmesh)
        mesh.Curve(3)
        core_vol = ngs.Integrate(1.0, mesh, definedon=mesh.Materials("default"))
        expected = sphere_volume(r - total_h)
        assert core_vol == pytest.approx(expected, rel=0.01), \
            f"Multi-layer core volume {core_vol} != expected {expected}"

    def test_sphere_multilayer_shell_volume(self):
        """Total BL shell (all layers combined) must have correct volume."""
        ngs = pytest.importorskip("ngsolve")
        r = 1.0
        thicknesses = [0.02, 0.03]
        total_h = sum(thicknesses)
        geo, bl, maxh = make_sphere_with_bl(r=r, bl_thickness=thicknesses, maxh=0.4)
        ngmesh = geo.GenerateMesh(maxh=maxh, boundary_layers=bl)
        mesh = ngs.Mesh(ngmesh)
        mesh.Curve(4)
        layer_vol = ngs.Integrate(1.0, mesh, definedon=mesh.Materials("layer"))
        expected_shell = sphere_shell_volume(r, r - total_h)
        assert layer_vol == pytest.approx(expected_shell, rel=0.05), \
            f"Multi-layer shell volume {layer_vol} != expected {expected_shell}"
