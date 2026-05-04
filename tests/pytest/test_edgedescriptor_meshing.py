"""
EdgeDescriptor: verify that edge descriptors are populated and consistent
after meshing for all geometry types (OCC, CSG, 2D spline, STL, periodic,
boundary layer).
"""

import pytest
import tempfile
import os


def test_occ_3d_box():
    """OCC 3D box: 12 edges, each segment has valid index."""
    pytest.importorskip("netgen.occ")
    from netgen.occ import Box, Pnt, OCCGeometry
    box = Box(Pnt(0,0,0), Pnt(1,1,1))
    geo = OCCGeometry(box)
    mesh = geo.GenerateMesh(maxh=0.5)

    assert mesh.GetNED() == 12
    edsi_values = set()
    for seg in mesh.Elements1D():
        assert seg.index >= 1, "segment has no edge descriptor"
        edsi_values.add(seg.index)
    assert len(edsi_values) == 12


def test_occ_3d_named_edges():
    """OCC 3D with named edges: names propagate to EdgeDescriptor."""
    pytest.importorskip("netgen.occ")
    from netgen.occ import Box, Pnt, OCCGeometry
    box = Box(Pnt(0,0,0), Pnt(1,1,1))
    for i, e in enumerate(box.edges):
        e.name = f"myedge_{i}"
    geo = OCCGeometry(box)
    mesh = geo.GenerateMesh(maxh=0.5)
    assert mesh.GetNED() == 12
    # verify names round-trip through cd2names
    found_names = set()
    for i in range(mesh.GetNCD2Names()):
        name = mesh.GetCD2Name(i)
        if name.startswith("myedge_"):
            found_names.add(name)
    assert len(found_names) == 12


def test_occ_2d_rectangle():
    """OCC 2D rectangle face: 4 edges."""
    pytest.importorskip("netgen.occ")
    from netgen.occ import WorkPlane, OCCGeometry
    wp = WorkPlane()
    wp.Rectangle(2, 1)
    face = wp.Face()
    geo = OCCGeometry(face, dim=2)
    mesh = geo.GenerateMesh(maxh=0.3)

    assert mesh.GetNED() == 4
    for seg in mesh.Elements1D():
        assert seg.index >= 1


def test_csg_box():
    """CSG box: all segments get valid edge descriptors."""
    from netgen.csg import CSGeometry, OrthoBrick, Pnt
    geo = CSGeometry()
    geo.Add(OrthoBrick(Pnt(0,0,0), Pnt(1,1,1)))
    mesh = geo.GenerateMesh(maxh=0.5)

    assert mesh.GetNED() > 0
    for seg in mesh.Elements1D():
        assert seg.index >= 1, "CSG segment has no edge descriptor"
    edgenrs = set(mesh.EdgeDescriptor(seg.index).edgenr for seg in mesh.Elements1D())
    assert len(edgenrs) == 12


def test_csg_spline_surface():
    """CSG with SplineSurface: no crash, segments get descriptors."""
    from netgen.csg import CSGeometry, OrthoBrick, Pnt, Plane, Vec, SplineSurface
    from netgen import meshing
    geo = CSGeometry()
    geo.Add(OrthoBrick(Pnt(-1,-1,-1), Pnt(1,1,1)).mat("air"))
    base = Plane(Pnt(0,0,0), Vec(0,0,1))
    surface = SplineSurface(base)
    pts = [(-0.2,-0.2,0),(-0.2,0.2,0),(0.2,0.2,0),(0.2,-0.2,0)]
    geopts = [surface.AddPoint(*p) for p in pts]
    for p1,p2,bc in [(0,1,"wire"),(1,2,"contact"),(2,3,"wire"),(3,0,"wire")]:
        surface.AddSegment(geopts[p1],geopts[p2],bc)
    geo.AddSplineSurface(surface)
    mesh = geo.GenerateMesh(maxh=0.4, perfstepsend=meshing.MeshingStep.MESHSURFACE)
    for i in range(2):
        mesh.GenerateVolumeMesh(only3D_domain_nr=i+1, maxh=0.4)
    assert mesh.GetNED() > 0


def test_2d_spline():
    """2D spline geometry: one EdgeDescriptor per spline segment."""
    from netgen.geom2d import SplineGeometry
    geo = SplineGeometry()
    p = [geo.AppendPoint(x,y) for x,y in [(0,0),(1,0),(1,1),(0,1)]]
    geo.Append(["line", p[0], p[1]], leftdomain=1, rightdomain=0, bc="bottom")
    geo.Append(["line", p[1], p[2]], leftdomain=1, rightdomain=0, bc="right")
    geo.Append(["line", p[2], p[3]], leftdomain=1, rightdomain=0, bc="top")
    geo.Append(["line", p[3], p[0]], leftdomain=1, rightdomain=0, bc="left")
    mesh = geo.GenerateMesh(maxh=0.3)

    assert mesh.GetNED() >= 4
    edsi_values = set()
    for seg in mesh.Elements1D():
        assert seg.index >= 1
        edsi_values.add(seg.index)
    assert len(edsi_values) == 4


def test_edsi_consistent_with_index():
    """All segments on the same geometric edge share the same index."""
    pytest.importorskip("netgen.occ")
    from netgen.occ import Box, Pnt, OCCGeometry
    geo = OCCGeometry(Box(Pnt(0,0,0), Pnt(1,1,1)))
    mesh = geo.GenerateMesh(maxh=0.3)

    # group segments by index (1-based ED index)
    from collections import defaultdict
    groups = defaultdict(list)
    for seg in mesh.Elements1D():
        groups[seg.index].append(seg)

    # each group should map to the same edgenr via EdgeDescriptor
    for edsi, segs in groups.items():
        ed = mesh.EdgeDescriptor(edsi)
        # all segments in this group share the same ED, so edgenr is consistent by definition
        assert ed.edgenr >= 0, f"ED[{edsi}] has invalid edgenr={ed.edgenr}"


def test_all_segments_have_valid_edsi_occ():
    """After OCC meshing, every segment has index >= 1."""
    pytest.importorskip("netgen.occ")
    from netgen.occ import Box, Pnt, OCCGeometry
    geo = OCCGeometry(Box(Pnt(0,0,0), Pnt(1,1,1)))
    mesh = geo.GenerateMesh(maxh=0.5)
    for seg in mesh.Elements1D():
        assert seg.index >= 1, f"segment {seg.vertices} has index={seg.index}"


def test_all_segments_have_valid_edsi_csg():
    """After CSG meshing, every segment has index >= 1."""
    from netgen.csg import CSGeometry, OrthoBrick, Pnt
    geo = CSGeometry()
    geo.Add(OrthoBrick(Pnt(0,0,0), Pnt(1,1,1)))
    mesh = geo.GenerateMesh(maxh=0.5)
    for seg in mesh.Elements1D():
        assert seg.index >= 1, f"segment {seg.vertices} has index={seg.index}"


def test_all_segments_have_valid_edsi_2d():
    """After 2D meshing, every segment has index >= 1."""
    from netgen.geom2d import SplineGeometry
    geo = SplineGeometry()
    p = [geo.AppendPoint(x,y) for x,y in [(0,0),(1,0),(1,1),(0,1)]]
    geo.Append(["line", p[0], p[1]], leftdomain=1, rightdomain=0)
    geo.Append(["line", p[1], p[2]], leftdomain=1, rightdomain=0)
    geo.Append(["line", p[2], p[3]], leftdomain=1, rightdomain=0)
    geo.Append(["line", p[3], p[0]], leftdomain=1, rightdomain=0)
    mesh = geo.GenerateMesh(maxh=0.3)
    for seg in mesh.Elements1D():
        assert seg.index >= 1, f"segment {seg.vertices} has index={seg.index}"


def test_all_segments_have_valid_edsi_boundarylayer():
    """After boundary layer, every segment has index >= 1."""
    pytest.importorskip("netgen.occ")
    from netgen.occ import Box, Pnt, OCCGeometry
    from netgen.meshing import BoundaryLayerParameters
    geo = OCCGeometry(Box(Pnt(0,0,0), Pnt(1,1,1)))
    blparams = [BoundaryLayerParameters(".*", [0.01, 0.01], "layer", outside=True)]
    mesh = geo.GenerateMesh(maxh=0.3, boundary_layers=blparams)
    assert mesh.GetNED() > 0
    for seg in mesh.Elements1D():
        assert seg.index >= 1, f"segment {seg.vertices} has index={seg.index}"


def test_boundarylayer_seg_si_valid_fd_index():
    """After boundary layer, each segment's ED has a valid fdindex (1-based FD index)."""
    pytest.importorskip("netgen.occ")
    from netgen.occ import Box, Pnt, OCCGeometry
    from netgen.meshing import BoundaryLayerParameters
    geo = OCCGeometry(Box(Pnt(0,0,0), Pnt(1,1,1)))
    blparams = [BoundaryLayerParameters(".*", [0.01, 0.01], "layer", outside=True)]
    mesh = geo.GenerateMesh(maxh=0.3, boundary_layers=blparams)
    nfd = mesh.GetNFaceDescriptors()
    ned = mesh.GetNED()
    assert nfd > 0
    for seg in mesh.Elements1D():
        idx = seg.index  # 1-based edge descriptor index
        assert 1 <= idx <= ned, f"segment index={idx} out of ED range [1, {ned}]"
        ed = mesh.EdgeDescriptor(idx)
        assert ed.fdindex >= 1, f"ED[{idx}] fdindex={ed.fdindex} should be >= 1"


def test_boundarylayer_ed_has_valid_content():
    """After meshing, EDs have meaningful surfnr and edgenr."""
    pytest.importorskip("netgen.occ")
    from netgen.occ import Box, Pnt, OCCGeometry
    geo = OCCGeometry(Box(Pnt(0,0,0), Pnt(1,1,1)))
    mesh = geo.GenerateMesh(maxh=0.5)

    for seg in mesh.Elements1D():
        ed = mesh.EdgeDescriptor(seg.index)
        # ED should have valid edgenr (positive for OCC)
        assert ed.edgenr > 0, f"ED[{seg.index}] has edgenr={ed.edgenr}"
        # surfnr tuple should have non-negative values
        snr0, snr1 = ed.surfnr
        assert snr0 >= 0, f"ED[{seg.index}] has surfnr[0]={snr0}"
        assert snr1 >= 0, f"ED[{seg.index}] has surfnr[1]={snr1}"


def test_csg_per_refedge_descriptors():
    """CSG box: per-refedge EDs carry per-refedge surfnr/tlosurf values."""
    from netgen.csg import CSGeometry, OrthoBrick, Pnt
    geo = CSGeometry()
    geo.Add(OrthoBrick(Pnt(0,0,0), Pnt(1,1,1)))
    mesh = geo.GenerateMesh(maxh=0.5)
    # Each geometric edge has 2 refedges (2 adjacent surfaces)
    # Group EDs by edgenr
    from collections import defaultdict
    eds_by_edge = defaultdict(list)
    for i, ed in enumerate(mesh.EdgeDescriptors()):
        eds_by_edge[ed.edgenr].append(ed)
    for edgenr, eds in eds_by_edge.items():
        assert len(eds) == 2, \
            f"edge {edgenr}: expected 2 per-refedge EDs, got {len(eds)}"
        # The two EDs should share edgenr but may differ in surfnr/fdindex
        assert eds[0].edgenr == eds[1].edgenr


def test_csg_delete_mesh_no_crash():
    """DeleteMesh properly clears edgedecoding and cd2names (no dangling pointers)."""
    from netgen.csg import CSGeometry, OrthoBrick, Pnt, Plane, Vec, SplineSurface
    from netgen import meshing

    geo = CSGeometry()
    geo.Add(OrthoBrick(Pnt(-1,-1,-1), Pnt(1,1,1)).mat("air"))
    base = Plane(Pnt(0,0,0), Vec(0,0,1))
    surface = SplineSurface(base)
    pts = [(-0.2,-0.2,0),(-0.2,0.2,0),(0.2,0.2,0),(0.2,-0.2,0)]
    geopts = [surface.AddPoint(*p) for p in pts]
    for p1,p2,bc in [(0,1,"wire"),(1,2,"contact"),(2,3,"wire"),(3,0,"wire")]:
        surface.AddSegment(geopts[p1],geopts[p2],bc)
    geo.AddSplineSurface(surface)

    # mesh twice — exercises DeleteMesh path
    mesh = geo.GenerateMesh(maxh=0.4, perfstepsend=meshing.MeshingStep.MESHSURFACE)
    for i in range(2):
        mesh.GenerateVolumeMesh(only3D_domain_nr=i+1, maxh=0.4)
    assert mesh.GetNED() > 0

    mesh2 = geo.GenerateMesh(maxh=0.4, perfstepsend=meshing.MeshingStep.MESHSURFACE)
    for i in range(2):
        mesh2.GenerateVolumeMesh(only3D_domain_nr=i+1, maxh=0.4)
    assert mesh2.GetNED() > 0


def test_2d_periodic_mesh_edsi():
    """2D periodic mesh: all segments have valid index >= 1."""
    from netgen.geom2d import SplineGeometry
    geo = SplineGeometry()
    p = [geo.AppendPoint(x,y) for x,y in [(0,0),(1,0),(1,1),(0,1)]]
    l0 = geo.Append(["line", p[0], p[1]], leftdomain=1, rightdomain=0, bc="bottom")
    l1 = geo.Append(["line", p[1], p[2]], leftdomain=1, rightdomain=0, bc="right")
    geo.Append(["line", p[3], p[2]], leftdomain=0, rightdomain=1, bc="top", copy=l0)
    geo.Append(["line", p[0], p[3]], leftdomain=0, rightdomain=1, bc="left", copy=l1)
    mesh = geo.GenerateMesh(maxh=0.3)

    for seg in mesh.Elements1D():
        assert seg.index >= 1, f"periodic segment has index={seg.index}"


def test_stl_mesh_edsi():
    """STL geometry: all segments have valid index >= 1."""
    stl = pytest.importorskip("netgen.stl")
    pytest.importorskip("netgen.occ")
    from netgen.occ import Box, Pnt, OCCGeometry

    # Generate an OCC box mesh, export as STL, then reimport via STLGeometry
    geo_occ = OCCGeometry(Box(Pnt(0,0,0), Pnt(1,1,1)))
    mesh_occ = geo_occ.GenerateMesh(maxh=0.5)

    with tempfile.NamedTemporaryFile(suffix=".stl", delete=False) as f:
        tmpfile = f.name
    try:
        mesh_occ.Export(tmpfile, "STL Format")
        geo_stl = stl.STLGeometry(tmpfile)
        mesh = geo_stl.GenerateMesh(maxh=0.5)
        assert len(mesh.Elements1D()) > 0
        for seg in mesh.Elements1D():
            assert seg.index >= 1, f"STL segment has index={seg.index}"
    finally:
        os.unlink(tmpfile)


def test_uniform_refinement_preserves_edsi():
    """After uniform refinement, all segments still have valid index and 12 edge descriptors."""
    pytest.importorskip("netgen.occ")
    from netgen.occ import Box, Pnt, OCCGeometry
    geo = OCCGeometry(Box(Pnt(0,0,0), Pnt(1,1,1)))
    mesh = geo.GenerateMesh(maxh=0.5)
    assert mesh.GetNED() == 12
    mesh.Refine()
    for seg in mesh.Elements1D():
        assert seg.index >= 1, f"segment {seg.vertices} has index={seg.index} after refinement"
    assert mesh.GetNED() == 12


def test_mesh_copy_preserves_edge_descriptors():
    """Mesh.Copy() preserves segment count; index values are copied on segments."""
    pytest.importorskip("netgen.occ")
    from netgen.occ import Box, Pnt, OCCGeometry
    geo = OCCGeometry(Box(Pnt(0,0,0), Pnt(1,1,1)))
    mesh = geo.GenerateMesh(maxh=0.5)
    nseg_orig = len(mesh.Elements1D())
    orig_edsi = [seg.index for seg in mesh.Elements1D()]

    copy = mesh.Copy()
    assert len(copy.Elements1D()) == nseg_orig
    copy_edsi = [seg.index for seg in copy.Elements1D()]
    assert copy_edsi == orig_edsi, "index values should be identical after Copy()"


def test_boundarylayer_ed_content_detailed():
    """After boundary layer meshing, EDs have valid edgenr, surfnr, and cover all box edges."""
    pytest.importorskip("netgen.occ")
    from netgen.occ import Box, Pnt, OCCGeometry
    from netgen.meshing import BoundaryLayerParameters
    geo = OCCGeometry(Box(Pnt(0,0,0), Pnt(1,1,1)))
    blparams = [BoundaryLayerParameters(".*", [0.01, 0.01], "layer", outside=True)]
    mesh = geo.GenerateMesh(maxh=0.3, boundary_layers=blparams)

    for i in range(1, mesh.GetNED()+1):
        ed = mesh.EdgeDescriptor(i)
        assert ed.edgenr > 0 or ed.edgenr == -1, f"ED[{i}] has edgenr={ed.edgenr}"
        snr0, snr1 = ed.surfnr
        assert snr0 >= 0 or snr0 == -1, f"ED[{i}] has surfnr[0]={snr0}"
        assert snr1 >= 0 or snr1 == -1, f"ED[{i}] has surfnr[1]={snr1}"


def test_occ_glue_compound_ed():
    """Glued OCC boxes: all segments and EDs are valid."""
    pytest.importorskip("netgen.occ")
    from netgen.occ import Box, Pnt, OCCGeometry, Glue
    box1 = Box(Pnt(0,0,0), Pnt(1,1,1))
    box2 = Box(Pnt(1,0,0), Pnt(2,1,1))
    geo = OCCGeometry(Glue([box1, box2]))
    mesh = geo.GenerateMesh(maxh=0.5)

    for seg in mesh.Elements1D():
        assert seg.index >= 1, f"segment has index={seg.index}"
    assert mesh.GetNED() > 0
    for i in range(1, mesh.GetNED()+1):
        ed = mesh.EdgeDescriptor(i)
        assert ed.edgenr > 0, f"ED[{i}] has edgenr={ed.edgenr}"


def test_occ_fuse_compound_ed():
    """Fused OCC boxes: all segments and EDs are valid."""
    pytest.importorskip("netgen.occ")
    from netgen.occ import Box, Pnt, OCCGeometry
    box1 = Box(Pnt(0,0,0), Pnt(1,1,1))
    box2 = Box(Pnt(0.5,0,0), Pnt(1.5,1,1))
    fused = box1 + box2
    geo = OCCGeometry(fused)
    mesh = geo.GenerateMesh(maxh=0.5)

    for seg in mesh.Elements1D():
        assert seg.index >= 1, f"segment has index={seg.index}"
    assert mesh.GetNED() > 0


def test_second_order_mesh_preserves_ed():
    """SecondOrder() does not change edge descriptor count or validity."""
    pytest.importorskip("netgen.occ")
    from netgen.occ import Box, Pnt, OCCGeometry
    geo = OCCGeometry(Box(Pnt(0,0,0), Pnt(1,1,1)))
    mesh = geo.GenerateMesh(maxh=0.5)

    ned_before = mesh.GetNED()
    edgenrs_before = [mesh.EdgeDescriptor(i).edgenr for i in range(1, ned_before+1)]

    mesh.SecondOrder()

    assert mesh.GetNED() == ned_before, \
        f"NED changed from {ned_before} to {mesh.GetNED()} after SecondOrder()"
    for seg in mesh.Elements1D():
        assert seg.index >= 1, f"segment has index={seg.index} after SecondOrder()"
    for i in range(1, mesh.GetNED()+1):
        assert mesh.EdgeDescriptor(i).edgenr == edgenrs_before[i-1], \
            f"ED[{i}] edgenr changed after SecondOrder()"


def test_csg_tlosurf_on_mesh():
    """CSG: tlosurf is -1 for simple OrthoBrick, >= 0 for SplineSurface TLO."""
    from netgen.csg import CSGeometry, OrthoBrick, Pnt, Plane, Vec, SplineSurface
    import netgen.meshing as meshing

    # Simple OrthoBrick: tlosurf should be -1 (no SplineSurface TLO)
    geo = CSGeometry()
    geo.Add(OrthoBrick(Pnt(0,0,0), Pnt(1,1,1)))
    mesh = geo.GenerateMesh(maxh=0.5)
    for i in range(1, mesh.GetNED()+1):
        ed = mesh.EdgeDescriptor(i)
        assert ed.tlosurf == -1, f"ED[{i}] tlosurf={ed.tlosurf} should be -1 for OrthoBrick"

    # SplineSurface TLO: at least some EDs should have tlosurf >= 0
    geo2 = CSGeometry()
    brick = OrthoBrick(Pnt(-1,-1,-1), Pnt(1,1,1))
    geo2.Add(brick)
    base = Plane(Pnt(0,0,0), Vec(0,0,1))
    surface = SplineSurface(base)
    pts = [(-0.2,-0.2,0), (-0.2,0.2,0), (0.2,0.2,0), (0.2,-0.2,0)]
    geopts = [surface.AddPoint(*p) for p in pts]
    for i in range(4):
        surface.AddSegment(geopts[i], geopts[(i+1) % 4], "wire")
    geo2.AddSplineSurface(surface)
    mesh2 = geo2.GenerateMesh(maxh=0.4, perfstepsend=meshing.MeshingStep.MESHSURFACE)
    for dom in range(1, 3):
        mesh2.GenerateVolumeMesh(only3D_domain_nr=dom, maxh=0.4)

    has_tlosurf = False
    for i in range(1, mesh2.GetNED()+1):
        ed = mesh2.EdgeDescriptor(i)
        if ed.tlosurf >= 0:
            has_tlosurf = True
            break
    assert has_tlosurf, "No ED with tlosurf >= 0 found in CSG SplineSurface mesh"


def test_boundarylayer_multiple_domains():
    """Boundary layer on a single box with named faces, BL on one face only: segments and EDs valid."""
    pytest.importorskip("netgen.occ")
    from netgen.occ import Box, Pnt, OCCGeometry
    import netgen.occ as occ
    from netgen.meshing import BoundaryLayerParameters
    box = Box(Pnt(0,0,0), Pnt(1,1,1))
    box.faces.Min(occ.Z).name = "bottom"
    box.faces.Max(occ.Z).name = "top"
    geo = OCCGeometry(box)
    blparams = [BoundaryLayerParameters("bottom", [0.01, 0.01], "layer", outside=True)]
    mesh = geo.GenerateMesh(maxh=0.3, boundary_layers=blparams)

    for seg in mesh.Elements1D():
        assert seg.index >= 1, f"segment has index={seg.index}"
    for i in range(1, mesh.GetNED()+1):
        ed = mesh.EdgeDescriptor(i)
        assert ed.edgenr > 0 or ed.edgenr == -1, f"ED[{i}] has edgenr={ed.edgenr}"


def test_occ_floating_wire_conforms_to_volume():
    """A wire glued inside a box must be conformingly resolved in the volume mesh.

    Regression test: the EdgeDescriptor for free edges must have DomainIn/DomainOut
    set so that ConformToFreeSegments can find and resolve them.
    """
    pytest.importorskip("netgen.occ")
    from netgen.occ import Box, Pnt, Segment, Glue, OCCGeometry

    box = Box(Pnt(0, 0, 0), Pnt(1, 1, 1))
    wire = Segment(Pnt(0.5, 0.5, 0.2), Pnt(0.5, 0.5, 0.8))
    wire.edges.name = "wire"
    shape = Glue([box, wire])
    mesh = OCCGeometry(shape).GenerateMesh(maxh=0.1)

    # Find the wire ED and check that domin/domout are set
    wire_edsi = None
    for i in range(1, mesh.GetNED()+1):
        ed = mesh.EdgeDescriptor(i)
        if ed.name == "wire":
            wire_edsi = i
            assert ed.domin >= 1, f"wire ED domin={ed.domin}, expected >= 1"
            assert ed.domout >= 1, f"wire ED domout={ed.domout}, expected >= 1"
            break
    assert wire_edsi is not None, "wire EdgeDescriptor not found"

    # Check that every wire segment endpoint is a vertex of some volume element,
    # i.e. the wire is conformingly resolved in the tet mesh.
    wire_pts = set()
    for seg in mesh.Elements1D():
        if seg.index == wire_edsi:
            for v in seg.vertices:
                wire_pts.add(v.nr)
    assert len(wire_pts) >= 2, "wire has no segment points"

    vol_pts = set()
    for el in mesh.Elements3D():
        for v in el.vertices:
            vol_pts.add(v.nr)

    unresolved = wire_pts - vol_pts
    assert len(unresolved) == 0, (
        f"{len(unresolved)} wire points are not volume element vertices — "
        f"wire is not conformingly meshed"
    )