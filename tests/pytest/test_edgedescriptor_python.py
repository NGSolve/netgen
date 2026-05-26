"""
EdgeDescriptor: Python bindings — EdgeDescriptor class, mesh access,
names, deprecated properties, numpy dtype, manual construction.
"""

import pytest


def test_edge_descriptor_python_class():
    """EdgeDescriptor Python class has all expected properties."""
    from netgen.meshing import EdgeDescriptor
    ed = EdgeDescriptor()
    assert ed.edgenr == -1
    assert ed.surfnr == (-1, -1)
    assert ed.name == "default"
    assert ed.singedge_left == 0.0
    assert ed.singedge_right == 0.0
    assert ed.tlosurf == -1

    ed.edgenr = 5
    ed.surfnr = (2, 3)
    ed.name = "myedge"
    ed.singedge_left = 1.5
    ed.singedge_right = 2.5
    ed.tlosurf = 7

    assert ed.edgenr == 5
    assert ed.surfnr == (2, 3)
    assert ed.name == "myedge"
    assert ed.singedge_left == 1.5
    assert ed.singedge_right == 2.5
    assert ed.tlosurf == 7
    assert "myedge" in repr(ed)


def test_mesh_edge_descriptors_access():
    """Mesh.EdgeDescriptors() returns iterable, Mesh.EdgeDescriptor(i) returns ref."""
    pytest.importorskip("netgen.occ")
    from netgen.occ import Box, Pnt, OCCGeometry
    geo = OCCGeometry(Box(Pnt(0,0,0), Pnt(1,1,1)))
    mesh = geo.GenerateMesh(maxh=0.5)

    eds = mesh.EdgeDescriptors()
    assert len(eds) == 12

    ed0 = mesh.EdgeDescriptor(1)
    assert ed0.edgenr >= 0
    assert ed0.surfnr[0] >= 0 or ed0.surfnr[1] >= 0

    # modifying via reference should stick
    old_name = ed0.name
    ed0.name = "test_rename"
    assert mesh.EdgeDescriptor(1).name == "test_rename"
    ed0.name = old_name  # restore


def test_mesh_add_edge_descriptor():
    """Mesh.Add(EdgeDescriptor) works."""
    from netgen.meshing import Mesh, EdgeDescriptor
    mesh = Mesh(dim=3)
    ed = EdgeDescriptor()
    ed.edgenr = 42
    ed.name = "added_edge"
    idx = mesh.Add(ed)
    assert idx == 1
    assert mesh.GetNED() == 1
    assert mesh.EdgeDescriptor(1).edgenr == 42
    assert mesh.EdgeDescriptor(1).name == "added_edge"


def test_edge_descriptor_names_from_mesh():
    """EdgeDescriptor names match cd2names for named edges."""
    pytest.importorskip("netgen.occ")
    from netgen.occ import Box, Pnt, OCCGeometry
    box = Box(Pnt(0,0,0), Pnt(1,1,1))
    for i, e in enumerate(box.edges):
        e.name = f"edge_{i}"
    geo = OCCGeometry(box)
    mesh = geo.GenerateMesh(maxh=0.5)

    for i in range(mesh.GetNED()):
        ed = mesh.EdgeDescriptor(i+1)
        cd2name = mesh.GetCD2Name(i)
        assert ed.name == cd2name, f"ed[{i}].name={ed.name} != cd2name={cd2name}"


def test_cd2names_syncs_to_edgedescriptor():
    """SetCD2Name updates EdgeDescriptor, GetCD2Name reads from EdgeDescriptor."""
    pytest.importorskip("netgen.occ")
    from netgen.occ import Box, Pnt, OCCGeometry
    geo = OCCGeometry(Box(Pnt(0,0,0), Pnt(1,1,1)))
    mesh = geo.GenerateMesh(maxh=0.5)

    # SetCD2Name is 1-based, GetCD2Name is 0-based
    mesh.SetCD2Name(1, "renamed_edge")
    assert mesh.EdgeDescriptor(1).name == "renamed_edge"
    assert mesh.GetCD2Name(0) == "renamed_edge"


def test_2d_spline_edgedescriptor_names():
    """2D spline BC names propagated to EdgeDescriptor after meshing."""
    from netgen.geom2d import SplineGeometry
    geo = SplineGeometry()
    p = [geo.AppendPoint(x,y) for x,y in [(0,0),(1,0),(1,1),(0,1)]]
    geo.Append(["line", p[0], p[1]], leftdomain=1, rightdomain=0, bc="bottom")
    geo.Append(["line", p[1], p[2]], leftdomain=1, rightdomain=0, bc="right")
    geo.Append(["line", p[2], p[3]], leftdomain=1, rightdomain=0, bc="top")
    geo.Append(["line", p[3], p[0]], leftdomain=1, rightdomain=0, bc="left")
    mesh = geo.GenerateMesh(maxh=0.3)

    expected_names = {"bottom", "right", "top", "left"}
    found_names = set()
    for i in range(1, mesh.GetNED()+1):
        ed = mesh.EdgeDescriptor(i)
        assert ed.edgenr >= 0
        # surfnr[0] is leftdom, surfnr[1] is rightdom
        assert ed.surfnr[0] == 1 or ed.surfnr[1] == 1  # at least one side is domain 1
        # EdgeDescriptor name should be the actual BC name, not "default"
        assert ed.name != "default", f"ED {i} has name 'default' instead of actual BC name"
        found_names.add(ed.name)
    assert found_names == expected_names, f"ED names {found_names} != expected {expected_names}"


def test_backward_compat_fields():
    """Segment index property is 0-based edge descriptor index; edgenr throws."""
    pytest.importorskip("netgen.occ")
    from netgen.occ import Box, Pnt, OCCGeometry
    geo = OCCGeometry(Box(Pnt(0,0,0), Pnt(1,1,1)))
    mesh = geo.GenerateMesh(maxh=0.5)

    seg = list(mesh.Elements1D())[0]
    assert seg.index >= 1, "1-based index should be >= 1"
    with pytest.raises(AttributeError):
        _ = seg.edgenr


def test_surfaces_property_throws():
    """Accessing seg.surfaces raises AttributeError directing to EdgeDescriptor."""
    pytest.importorskip("netgen.occ")
    from netgen.occ import Box, Pnt, OCCGeometry
    geo = OCCGeometry(Box(Pnt(0,0,0), Pnt(1,1,1)))
    mesh = geo.GenerateMesh(maxh=0.5)
    seg = list(mesh.Elements1D())[0]
    with pytest.raises(AttributeError, match="surfaces.*removed.*EdgeDescriptor"):
        _ = seg.surfaces


def test_singular_property_throws():
    """Accessing seg.singular raises AttributeError directing to EdgeDescriptor."""
    pytest.importorskip("netgen.occ")
    from netgen.occ import Box, Pnt, OCCGeometry
    geo = OCCGeometry(Box(Pnt(0,0,0), Pnt(1,1,1)))
    mesh = geo.GenerateMesh(maxh=0.5)
    seg = list(mesh.Elements1D())[0]
    with pytest.raises(AttributeError, match="singular.*EdgeDescriptor"):
        _ = seg.singular


def test_surfaces_constructor_param_ignored():
    """Segment(surfaces=[...]) is accepted but surfaces live on EdgeDescriptor."""
    from netgen.meshing import Mesh, Element1D, MeshPoint, Pnt
    mesh = Mesh(dim=2)
    p1 = mesh.Add(MeshPoint(Pnt(0, 0, 0)))
    p2 = mesh.Add(MeshPoint(Pnt(1, 0, 0)))
    # should not crash even though surfaces= is passed
    seg = Element1D([p1, p2], surfaces=[5, 6], index=1)
    mesh.Add(seg)
    # the segment is added; surfaces are not stored on segment
    assert len(list(mesh.Elements1D())) == 1


def test_numpy_element1d_edsi_field():
    """NumPy access to Elements1D has 'index' field, all >= 0."""
    pytest.importorskip("netgen.occ")
    from netgen.occ import Box, Pnt, OCCGeometry
    geo = OCCGeometry(Box(Pnt(0,0,0), Pnt(1,1,1)))
    mesh = geo.GenerateMesh(maxh=0.5)

    segs = mesh.Elements1D().NumPy()
    assert "index" in segs.dtype.names, f"fields: {segs.dtype.names}"
    indices = segs["index"]
    assert len(indices) > 0
    assert all(indices >= 1), f"found index values < 1 in numpy array"


def test_numpy_dtype_no_edgenr_field():
    """NumPy dtype for Elements1D should not contain 'edgenr' field."""
    pytest.importorskip("netgen.occ")
    from netgen.occ import Box, Pnt, OCCGeometry
    geo = OCCGeometry(Box(Pnt(0,0,0), Pnt(1,1,1)))
    mesh = geo.GenerateMesh(maxh=0.5)
    segs = mesh.Elements1D().NumPy()
    assert "edgenr" not in segs.dtype.names, \
        f"'edgenr' should be removed from numpy dtype, got fields: {segs.dtype.names}"


def test_numpy_index_field_matches_edsi():
    """NumPy 'index' field should match seg.index for every segment."""
    pytest.importorskip("netgen.occ")
    from netgen.occ import Box, Pnt, OCCGeometry
    geo = OCCGeometry(Box(Pnt(0,0,0), Pnt(1,1,1)))
    mesh = geo.GenerateMesh(maxh=0.5)
    segs_np = mesh.Elements1D().NumPy()
    segs_list = list(mesh.Elements1D())
    for i, seg in enumerate(segs_list):
        assert segs_np["index"][i] == seg.index, \
            f"seg[{i}]: numpy index={segs_np['index'][i]} != index={seg.index}"


def test_manual_mesh_construction_1d():
    """Manually constructed 1D mesh with Element1D has assigned index on segments."""
    from netgen.meshing import Mesh, MeshPoint, Element1D
    from netgen.libngpy._meshing import Pnt as MeshPnt

    mesh = Mesh(dim=1)
    pts = [mesh.Add(MeshPoint(MeshPnt(x, 0, 0))) for x in [0.0, 0.25, 0.5, 0.75]]

    for i in range(len(pts) - 1):
        mesh.Add(Element1D(index=1, vertices=[pts[i], pts[i+1]]))

    segs = list(mesh.Elements1D())
    assert len(segs) == 3
    for seg in segs:
        assert seg.index == 1


def test_fdindex_python_property():
    """EdgeDescriptor.fdindex is readable and writable from Python."""
    from netgen.meshing import EdgeDescriptor
    ed = EdgeDescriptor()
    assert ed.fdindex == -1, "default fdindex should be -1"
    ed.fdindex = 42
    assert ed.fdindex == 42