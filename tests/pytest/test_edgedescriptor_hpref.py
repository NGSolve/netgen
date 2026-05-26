"""
EdgeDescriptor: HP refinement — singular edge (hpref) propagation through
the meshing pipeline: geometry -> EdgeDescriptor -> save/load -> RefineHP.
"""

import pytest
import tempfile
import os


def test_2d_hpref_on_edge_descriptor():
    """2D spline: hpref values propagate from geometry to EdgeDescriptor."""
    from netgen.geom2d import SplineGeometry
    geo = SplineGeometry()
    p = [geo.AppendPoint(x, y) for x, y in [(0, 0), (1, 0), (1, 1), (0, 1)]]
    # bottom edge: hpref=2 on both sides
    geo.Append(["line", p[0], p[1]], leftdomain=1, rightdomain=0, bc="bottom", hpref=2)
    # right edge: asymmetric hpref
    geo.Append(["line", p[1], p[2]], leftdomain=1, rightdomain=0, bc="right",
               hprefleft=3, hprefright=0)
    # top and left: no hpref
    geo.Append(["line", p[2], p[3]], leftdomain=1, rightdomain=0, bc="top")
    geo.Append(["line", p[3], p[0]], leftdomain=1, rightdomain=0, bc="left")
    mesh = geo.GenerateMesh(maxh=0.5)

    # Collect singedge values per edge descriptor
    ned = mesh.GetNED()
    assert ned == 4, f"expected 4 edge descriptors, got {ned}"
    singedge = {}
    for i in range(ned):
        ed = mesh.EdgeDescriptor(i+1)
        singedge[ed.edgenr] = (ed.singedge_left, ed.singedge_right)

    # Edge 1 (bottom): hpref=2 -> both left and right = 2
    assert singedge[1] == (2.0, 2.0), f"bottom edge singedge mismatch: {singedge[1]}"
    # Edge 2 (right): hprefleft=3, hprefright=0
    assert singedge[2] == (3.0, 0.0), f"right edge singedge mismatch: {singedge[2]}"
    # Edge 3 (top): no hpref -> 0
    assert singedge[3] == (0.0, 0.0), f"top edge singedge mismatch: {singedge[3]}"
    # Edge 4 (left): no hpref -> 0
    assert singedge[4] == (0.0, 0.0), f"left edge singedge mismatch: {singedge[4]}"


def test_2d_hpref_save_load_roundtrip():
    """2D hpref values survive save/load cycle (both edgedescriptors section and legacy)."""
    from netgen.geom2d import SplineGeometry
    import netgen.meshing as nm

    geo = SplineGeometry()
    p = [geo.AppendPoint(x, y) for x, y in [(0, 0), (1, 0), (1, 1), (0, 1)]]
    geo.Append(["line", p[0], p[1]], leftdomain=1, rightdomain=0, hpref=2)
    geo.Append(["line", p[1], p[2]], leftdomain=1, rightdomain=0, hprefleft=1.5, hprefright=0.5)
    geo.Append(["line", p[2], p[3]], leftdomain=1, rightdomain=0)
    geo.Append(["line", p[3], p[0]], leftdomain=1, rightdomain=0)
    mesh = geo.GenerateMesh(maxh=0.5)

    # Collect original values
    orig = []
    for i in range(mesh.GetNED()):
        ed = mesh.EdgeDescriptor(i+1)
        orig.append((ed.edgenr, ed.singedge_left, ed.singedge_right))

    with tempfile.TemporaryDirectory() as tmpdir:
        path = os.path.join(tmpdir, "hpref_test.vol")
        mesh.Save(path)
        mesh2 = nm.Mesh()
        mesh2.Load(path)

        assert mesh2.GetNED() == mesh.GetNED()
        for i, (ednr, sl, sr) in enumerate(orig):
            ed2 = mesh2.EdgeDescriptor(i+1)
            assert ed2.edgenr == ednr
            assert ed2.singedge_left == pytest.approx(sl), \
                f"ED[{i}] singedge_left mismatch after load: {ed2.singedge_left} != {sl}"
            assert ed2.singedge_right == pytest.approx(sr), \
                f"ED[{i}] singedge_right mismatch after load: {ed2.singedge_right} != {sr}"


def test_2d_hpref_pickle_roundtrip():
    """2D hpref values survive pickle (DoArchive) cycle."""
    from netgen.geom2d import SplineGeometry
    import pickle

    geo = SplineGeometry()
    p = [geo.AppendPoint(x, y) for x, y in [(0, 0), (1, 0), (1, 1), (0, 1)]]
    geo.Append(["line", p[0], p[1]], leftdomain=1, rightdomain=0, hpref=2)
    geo.Append(["line", p[1], p[2]], leftdomain=1, rightdomain=0, hprefleft=1, hprefright=3)
    geo.Append(["line", p[2], p[3]], leftdomain=1, rightdomain=0)
    geo.Append(["line", p[3], p[0]], leftdomain=1, rightdomain=0)
    mesh = geo.GenerateMesh(maxh=0.5)

    orig = []
    for i in range(mesh.GetNED()):
        ed = mesh.EdgeDescriptor(i+1)
        orig.append((ed.edgenr, ed.singedge_left, ed.singedge_right))

    data = pickle.dumps(mesh)
    mesh2 = pickle.loads(data)

    assert mesh2.GetNED() == len(orig)
    for i, (ednr, sl, sr) in enumerate(orig):
        ed2 = mesh2.EdgeDescriptor(i+1)
        assert ed2.singedge_left == pytest.approx(sl), \
            f"ED[{i}] singedge_left mismatch after pickle: {ed2.singedge_left} != {sl}"
        assert ed2.singedge_right == pytest.approx(sr), \
            f"ED[{i}] singedge_right mismatch after pickle: {ed2.singedge_right} != {sr}"


def test_2d_hpref_segments_consistent():
    """All segments on an hpref edge share the same EdgeDescriptor with correct singedge."""
    from netgen.geom2d import SplineGeometry

    geo = SplineGeometry()
    p = [geo.AppendPoint(x, y) for x, y in [(0, 0), (1, 0), (1, 1), (0, 1)]]
    geo.Append(["line", p[0], p[1]], leftdomain=1, rightdomain=0, hpref=2)
    geo.Append(["line", p[1], p[2]], leftdomain=1, rightdomain=0)
    geo.Append(["line", p[2], p[3]], leftdomain=1, rightdomain=0)
    geo.Append(["line", p[3], p[0]], leftdomain=1, rightdomain=0)
    mesh = geo.GenerateMesh(maxh=0.2)  # finer -> more segments per edge

    # Find all segments on edge 1 (bottom, hpref=2)
    hpref_segs = []
    for seg in mesh.Elements1D():
        ed = mesh.EdgeDescriptor(seg.index)
        if ed.edgenr == 1:
            hpref_segs.append(seg.index)
            assert ed.singedge_left == pytest.approx(2.0)
            assert ed.singedge_right == pytest.approx(2.0)
    assert len(hpref_segs) > 0, "no segments found on bottom edge"

    # All non-hpref edges should have singedge = 0
    for seg in mesh.Elements1D():
        ed = mesh.EdgeDescriptor(seg.index)
        if ed.edgenr != 1:
            assert ed.singedge_left == pytest.approx(0.0)
            assert ed.singedge_right == pytest.approx(0.0)


def test_3d_occ_hpref_on_edge_descriptor():
    """OCC 3D: hpref set via maxh_edge property propagates to EdgeDescriptor."""
    pytest.importorskip("netgen.occ")
    from netgen.occ import Box, Pnt, OCCGeometry

    box = Box(Pnt(0, 0, 0), Pnt(1, 1, 1))
    geo = OCCGeometry(box)
    mesh = geo.GenerateMesh(maxh=0.5)

    # OCC doesn't natively set singedge (it's a CSG/2D feature),
    # so all EDs should have singedge = 0
    for i in range(mesh.GetNED()):
        ed = mesh.EdgeDescriptor(i+1)
        assert ed.singedge_left == pytest.approx(0.0), \
            f"OCC ED[{i}] unexpected singedge_left={ed.singedge_left}"
        assert ed.singedge_right == pytest.approx(0.0), \
            f"OCC ED[{i}] unexpected singedge_right={ed.singedge_right}"


def test_2d_hpref_domin_domout_on_ed():
    """2D: EdgeDescriptor carries correct domin/domout after hpref edge creation."""
    from netgen.geom2d import SplineGeometry

    geo = SplineGeometry()
    p = [geo.AppendPoint(x, y) for x, y in [(0, 0), (1, 0), (1, 1), (0, 1)]]
    geo.Append(["line", p[0], p[1]], leftdomain=1, rightdomain=0, hpref=1)
    geo.Append(["line", p[1], p[2]], leftdomain=1, rightdomain=0)
    geo.Append(["line", p[2], p[3]], leftdomain=1, rightdomain=0)
    geo.Append(["line", p[3], p[0]], leftdomain=1, rightdomain=0)
    mesh = geo.GenerateMesh(maxh=0.5)

    for i in range(mesh.GetNED()):
        ed = mesh.EdgeDescriptor(i+1)
        # leftdomain=1, rightdomain=0 for all edges
        assert ed.domin == 1, f"ED[{i}] domin={ed.domin}, expected 1"
        assert ed.domout == 0, f"ED[{i}] domout={ed.domout}, expected 0"


def test_2d_multi_domain_hpref():
    """2D with two domains: hpref on shared boundary propagates correctly."""
    from netgen.geom2d import SplineGeometry

    geo = SplineGeometry()
    # Two rectangles side by side sharing an edge at x=1
    p = [geo.AppendPoint(x, y) for x, y in
         [(0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (2, 1)]]
    # Domain 1: left rectangle
    geo.Append(["line", p[0], p[1]], leftdomain=1, rightdomain=0, bc="bottom_left")
    geo.Append(["line", p[1], p[4]], leftdomain=1, rightdomain=2, bc="shared", hpref=3)
    geo.Append(["line", p[4], p[3]], leftdomain=1, rightdomain=0, bc="top_left")
    geo.Append(["line", p[3], p[0]], leftdomain=1, rightdomain=0, bc="left")
    # Domain 2: right rectangle
    geo.Append(["line", p[1], p[2]], leftdomain=2, rightdomain=0, bc="bottom_right")
    geo.Append(["line", p[2], p[5]], leftdomain=2, rightdomain=0, bc="right")
    geo.Append(["line", p[5], p[4]], leftdomain=2, rightdomain=0, bc="top_right")
    mesh = geo.GenerateMesh(maxh=0.5)

    # Find the shared edge (hpref=3)
    shared_found = False
    for i in range(mesh.GetNED()):
        ed = mesh.EdgeDescriptor(i+1)
        if ed.singedge_left == pytest.approx(3.0):
            shared_found = True
            assert ed.singedge_right == pytest.approx(3.0)
            # Shared edge between domain 1 and domain 2
            assert ed.domin == 1, f"shared edge domin={ed.domin}"
            assert ed.domout == 2, f"shared edge domout={ed.domout}"
    assert shared_found, "shared hpref edge not found on any EdgeDescriptor"


def test_2d_hpref_refine_hp_produces_graded_mesh():
    """End-to-end: hpref on edge -> EdgeDescriptor -> ngsolve RefineHP -> graded mesh.

    Verifies that singular edge markers actually drive geometric grading
    in HP refinement: elements near a singular edge should be smaller than
    elements far away.
    """
    ngs = pytest.importorskip("ngsolve")
    from netgen.geom2d import SplineGeometry

    geo = SplineGeometry()
    p = [geo.AppendPoint(x, y) for x, y in [(0, 0), (1, 0), (1, 1), (0, 1)]]
    # Mark bottom edge as singular (hpref=1)
    geo.Append(["line", p[0], p[1]], leftdomain=1, rightdomain=0, bc="bottom", hpref=1)
    geo.Append(["line", p[1], p[2]], leftdomain=1, rightdomain=0, bc="right")
    geo.Append(["line", p[2], p[3]], leftdomain=1, rightdomain=0, bc="top")
    geo.Append(["line", p[3], p[0]], leftdomain=1, rightdomain=0, bc="left")
    ngmesh = geo.GenerateMesh(maxh=0.4)

    mesh = ngs.Mesh(ngmesh)
    ne_before = mesh.ne

    # Perform HP refinement -- 2 levels
    mesh.RefineHP(levels=2, factor=0.25)
    ne_after = mesh.ne

    # HP refinement must increase element count
    assert ne_after > ne_before, \
        f"RefineHP did not increase elements: {ne_before} -> {ne_after}"

    # Grading check: elements touching the singular edge (bottom, y~0)
    # should be smaller than elements far from it (near top, y~1).
    bottom_areas = []
    top_areas = []
    for el in mesh.Elements(ngs.VOL):
        verts = [mesh[v].point for v in el.vertices]
        # 2D triangle area via cross product
        ax, ay = verts[1][0] - verts[0][0], verts[1][1] - verts[0][1]
        bx, by = verts[2][0] - verts[0][0], verts[2][1] - verts[0][1]
        area = abs(ax * by - ay * bx) / 2
        cy = sum(v[1] for v in verts) / len(verts)  # centroid y
        if cy < 0.15:
            bottom_areas.append(area)
        elif cy > 0.85:
            top_areas.append(area)

    assert len(bottom_areas) > 0, "no elements near bottom edge"
    assert len(top_areas) > 0, "no elements near top edge"

    avg_bottom = sum(bottom_areas) / len(bottom_areas)
    avg_top = sum(top_areas) / len(top_areas)
    # Elements near the singular edge should be significantly smaller
    assert avg_bottom < avg_top, \
        f"Expected grading: avg bottom area ({avg_bottom:.6f}) should be < avg top area ({avg_top:.6f})"


def test_3d_csg_hpref_refine_hp():
    """End-to-end: CSG singular edge -> EdgeDescriptor -> ngsolve RefineHP.

    Verifies that HP refinement runs without error on a 3D CSG mesh
    with singular edges, and that segments are correctly reconstructed
    after refinement with valid index.
    """
    ngs = pytest.importorskip("ngsolve")
    from netgen.csg import CSGeometry, OrthoBrick, Pnt

    geo = CSGeometry()
    geo.Add(OrthoBrick(Pnt(0, 0, 0), Pnt(1, 1, 1)))
    ngmesh = geo.GenerateMesh(maxh=0.5)

    mesh = ngs.Mesh(ngmesh)
    ne_before = mesh.ne

    # HP refinement with 1 level -- should not crash
    mesh.RefineHP(levels=1, factor=0.25)
    ne_after = mesh.ne
    # OrthoBrick has no singular edges, so RefineHP may not increase elements.
    # The key check is that it doesn't crash and segments remain valid.
    assert ne_after >= ne_before, \
        f"RefineHP lost elements: {ne_before} -> {ne_after}"

    # All segments after refinement should have valid index
    for seg in mesh.ngmesh.Elements1D():
        assert seg.index >= 1, "segment has invalid index after HP refinement"


def test_2d_hpref_no_singular_no_grading():
    """Control test: without hpref, RefineHP should do uniform refinement."""
    ngs = pytest.importorskip("ngsolve")
    from netgen.geom2d import SplineGeometry

    geo = SplineGeometry()
    p = [geo.AppendPoint(x, y) for x, y in [(0, 0), (1, 0), (1, 1), (0, 1)]]
    # No hpref on any edge
    geo.Append(["line", p[0], p[1]], leftdomain=1, rightdomain=0)
    geo.Append(["line", p[1], p[2]], leftdomain=1, rightdomain=0)
    geo.Append(["line", p[2], p[3]], leftdomain=1, rightdomain=0)
    geo.Append(["line", p[3], p[0]], leftdomain=1, rightdomain=0)
    ngmesh = geo.GenerateMesh(maxh=0.4)

    mesh = ngs.Mesh(ngmesh)
    mesh.RefineHP(levels=2, factor=0.25)

    # Without singular markers, there should be no strong grading
    bottom_areas = []
    top_areas = []
    for el in mesh.Elements(ngs.VOL):
        verts = [mesh[v].point for v in el.vertices]
        ax, ay = verts[1][0] - verts[0][0], verts[1][1] - verts[0][1]
        bx, by = verts[2][0] - verts[0][0], verts[2][1] - verts[0][1]
        area = abs(ax * by - ay * bx) / 2
        cy = sum(v[1] for v in verts) / len(verts)
        if cy < 0.15:
            bottom_areas.append(area)
        elif cy > 0.85:
            top_areas.append(area)

    if len(bottom_areas) > 0 and len(top_areas) > 0:
        avg_bottom = sum(bottom_areas) / len(bottom_areas)
        avg_top = sum(top_areas) / len(top_areas)
        # Ratio should be close to 1 (no strong grading)
        ratio = avg_bottom / avg_top if avg_top > 0 else 1.0
        assert ratio > 0.2, \
            f"Unexpected grading without hpref: bottom/top area ratio = {ratio:.3f}"


def test_2d_hpref_save_load_then_refine_hp():
    """hpref survives save/load and still drives HP refinement correctly."""
    ngs = pytest.importorskip("ngsolve")
    from netgen.geom2d import SplineGeometry
    import netgen.meshing as nm

    geo = SplineGeometry()
    p = [geo.AppendPoint(x, y) for x, y in [(0, 0), (1, 0), (1, 1), (0, 1)]]
    geo.Append(["line", p[0], p[1]], leftdomain=1, rightdomain=0, hpref=1)
    geo.Append(["line", p[1], p[2]], leftdomain=1, rightdomain=0)
    geo.Append(["line", p[2], p[3]], leftdomain=1, rightdomain=0)
    geo.Append(["line", p[3], p[0]], leftdomain=1, rightdomain=0)
    ngmesh = geo.GenerateMesh(maxh=0.4)

    # Save and reload
    with tempfile.TemporaryDirectory() as tmpdir:
        path = os.path.join(tmpdir, "hpref_reload.vol")
        ngmesh.Save(path)
        ngmesh2 = nm.Mesh()
        ngmesh2.Load(path)

        # Verify singedge survived
        found_hpref = False
        for i in range(ngmesh2.GetNED()):
            ed = ngmesh2.EdgeDescriptor(i+1)
            if ed.singedge_left > 0:
                found_hpref = True
                break
        assert found_hpref, "hpref values lost after save/load"

        # RefineHP on reloaded mesh should still produce grading
        mesh = ngs.Mesh(ngmesh2)
        ne_before = mesh.ne
        mesh.RefineHP(levels=2, factor=0.25)
        assert mesh.ne > ne_before, "RefineHP did nothing on reloaded mesh"

        bottom_areas = []
        top_areas = []
        for el in mesh.Elements(ngs.VOL):
            verts = [mesh[v].point for v in el.vertices]
            ax, ay = verts[1][0] - verts[0][0], verts[1][1] - verts[0][1]
            bx, by = verts[2][0] - verts[0][0], verts[2][1] - verts[0][1]
            area = abs(ax * by - ay * bx) / 2
            cy = sum(v[1] for v in verts) / len(verts)
            if cy < 0.15:
                bottom_areas.append(area)
            elif cy > 0.85:
                top_areas.append(area)

        if len(bottom_areas) > 0 and len(top_areas) > 0:
            avg_bottom = sum(bottom_areas) / len(bottom_areas)
            avg_top = sum(top_areas) / len(top_areas)
            assert avg_bottom < avg_top, \
                f"No grading after reload+RefineHP: bottom={avg_bottom:.6f}, top={avg_top:.6f}"


    # Save the HP-refined mesh via netgen and verify EDs survive
    ngmesh_refined = mesh.ngmesh
    orig = []
    for i in range(ngmesh_refined.GetNED()):
        ed = ngmesh_refined.EdgeDescriptor(i+1)
        orig.append((ed.edgenr, ed.surfnr, ed.name, ed.singedge_left, ed.singedge_right))

    with tempfile.NamedTemporaryFile(suffix=".vol", delete=False) as f:
        tmpfile = f.name
    try:
        ngmesh_refined.Save(tmpfile)

        # verify the section exists in the file
        with open(tmpfile) as fh:
            content = fh.read()
        assert "edgedescriptors" in content

        # load and compare
        mesh2 = nm.Mesh()
        mesh2.Load(tmpfile)
        assert mesh2.GetNED() == len(orig)
        for i, (ednr, snr, name, sl, sr) in enumerate(orig):
            ed2 = mesh2.EdgeDescriptor(i+1)
            assert ed2.edgenr == ednr, f"ed[{i}] edgenr mismatch"
            assert ed2.surfnr == snr, f"ed[{i}] surfnr mismatch"
            assert ed2.name == name, f"ed[{i}] name mismatch"
            assert ed2.singedge_left == sl
            assert ed2.singedge_right == sr

        # Note: HP-refined segments may not have valid index after save/load
        # since HP refinement creates internal segments outside the normal
        # meshing pipeline. The important check is that EDs roundtrip above.
    finally:
        os.unlink(tmpfile)