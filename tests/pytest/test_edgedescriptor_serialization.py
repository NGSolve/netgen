"""
EdgeDescriptor: serialization — save/load roundtrips, backward compat,
medit export/import, save-compare, boundary name survival.
"""

import pytest
import tempfile
import os


def test_save_load_roundtrip():
    """Save and load preserves segment count and reconstructs edge descriptors."""
    pytest.importorskip("netgen.occ")
    from netgen.occ import Box, Pnt, OCCGeometry
    from netgen.meshing import Mesh

    geo = OCCGeometry(Box(Pnt(0,0,0), Pnt(1,1,1)))
    mesh = geo.GenerateMesh(maxh=0.5)
    nseg_before = len(mesh.Elements1D())
    ned_before = mesh.GetNED()

    with tempfile.NamedTemporaryFile(suffix=".vol", delete=False) as f:
        tmpfile = f.name
    try:
        mesh.Save(tmpfile)
        mesh2 = Mesh()
        mesh2.Load(tmpfile)
        assert len(mesh2.Elements1D()) == nseg_before
        assert mesh2.GetNED() > 0
        for seg in mesh.Elements1D():
            assert seg.index >= 1
    finally:
        os.unlink(tmpfile)


def test_save_load_save_compare():
    """Save -> load -> save produces identical .vol output (excluding face_colours)."""
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
    mesh.SetGeometry(None)

    with tempfile.TemporaryDirectory() as tmpdir:
        f1 = os.path.join(tmpdir, "mesh1.vol")
        f2 = os.path.join(tmpdir, "mesh2.vol")
        mesh.Save(f1)
        mesh2 = meshing.Mesh()
        mesh2.Load(f1)
        mesh2.Save(f2)

        with open(f1) as a, open(f2) as b:
            lines_a = a.readlines()
            lines_b = b.readlines()
        # exclude face_colours section
        for lines in (lines_a, lines_b):
            for i, line in enumerate(lines):
                if line.strip().startswith("face_colours"):
                    del lines[i:]
                    break
        assert lines_a == lines_b, "save/load/save mismatch"


def test_save_load_edgedescriptors_section():
    """Save writes edgedescriptors section; load reads it back exactly."""
    pytest.importorskip("netgen.occ")
    from netgen.occ import Box, Pnt, OCCGeometry
    from netgen.meshing import Mesh
    import tempfile, os

    geo = OCCGeometry(Box(Pnt(0,0,0), Pnt(1,1,1)))
    mesh = geo.GenerateMesh(maxh=0.5)

    # collect original descriptor data
    orig = []
    for i in range(mesh.GetNED()):
        ed = mesh.EdgeDescriptor(i+1)
        orig.append((ed.edgenr, ed.surfnr, ed.name, ed.singedge_left, ed.singedge_right))

    with tempfile.NamedTemporaryFile(suffix=".vol", delete=False) as f:
        tmpfile = f.name
    try:
        mesh.Save(tmpfile)

        # verify the section exists in the file
        with open(tmpfile) as fh:
            content = fh.read()
        assert "edgedescriptors" in content

        # load and compare
        mesh2 = Mesh()
        mesh2.Load(tmpfile)
        assert mesh2.GetNED() == len(orig)
        for i, (ednr, snr, nm, sl, sr) in enumerate(orig):
            ed2 = mesh2.EdgeDescriptor(i+1)
            assert ed2.edgenr == ednr, f"ed[{i}] edgenr mismatch"
            assert ed2.surfnr == snr, f"ed[{i}] surfnr mismatch"
            assert ed2.name == nm, f"ed[{i}] name mismatch"
            assert ed2.singedge_left == sl
            assert ed2.singedge_right == sr
    finally:
        os.unlink(tmpfile)


def test_save_load_preserves_named_edges():
    """Named edges survive save/load via edgedescriptors section."""
    pytest.importorskip("netgen.occ")
    from netgen.occ import Box, Pnt, OCCGeometry
    from netgen.meshing import Mesh
    import tempfile, os

    box = Box(Pnt(0,0,0), Pnt(1,1,1))
    for i, e in enumerate(box.edges):
        e.name = f"wire_{i}"
    geo = OCCGeometry(box)
    mesh = geo.GenerateMesh(maxh=0.5)

    with tempfile.NamedTemporaryFile(suffix=".vol", delete=False) as f:
        tmpfile = f.name
    try:
        mesh.Save(tmpfile)
        mesh2 = Mesh()
        mesh2.Load(tmpfile)

        for i in range(mesh2.GetNED()):
            assert mesh2.EdgeDescriptor(i+1).name.startswith("wire_")
            assert mesh2.GetCD2Name(i).startswith("wire_")
    finally:
        os.unlink(tmpfile)


def test_all_segments_have_valid_edsi_save_load():
    """After save/load roundtrip, every segment has index >= 0."""
    pytest.importorskip("netgen.occ")
    from netgen.occ import Box, Pnt, OCCGeometry
    from netgen.meshing import Mesh
    import tempfile, os

    geo = OCCGeometry(Box(Pnt(0,0,0), Pnt(1,1,1)))
    mesh = geo.GenerateMesh(maxh=0.5)

    with tempfile.NamedTemporaryFile(suffix=".vol", delete=False) as f:
        tmpfile = f.name
    try:
        mesh.Save(tmpfile)
        mesh2 = Mesh()
        mesh2.Load(tmpfile)
        for seg in mesh2.Elements1D():
            assert seg.index >= 1, f"loaded segment has index={seg.index}"
    finally:
        os.unlink(tmpfile)


def test_save_load_boundary_names_survive():
    """OCC box with named faces: BC names and valid index survive save/load."""
    pytest.importorskip("netgen.occ")
    from netgen.occ import Box, Pnt, OCCGeometry
    from netgen.meshing import Mesh
    import tempfile, os

    box = Box(Pnt(0,0,0), Pnt(1,1,1))
    for i, f in enumerate(box.faces):
        f.name = f"face_{i}"
    geo = OCCGeometry(box)
    mesh = geo.GenerateMesh(maxh=0.5)

    nfd = mesh.GetNFaceDescriptors()
    orig_bcnames = [mesh.GetBCName(i) for i in range(nfd)]
    orig_cd2names = [mesh.GetCD2Name(i) for i in range(mesh.GetNED())]

    with tempfile.NamedTemporaryFile(suffix=".vol", delete=False) as f:
        tmpfile = f.name
    try:
        mesh.Save(tmpfile)
        mesh2 = Mesh()
        mesh2.Load(tmpfile)

        assert mesh2.GetNFaceDescriptors() == nfd
        for i in range(nfd):
            assert mesh2.GetBCName(i) == orig_bcnames[i], f"BCName mismatch at {i}"

        for i in range(mesh2.GetNED()):
            assert mesh2.GetCD2Name(i) == orig_cd2names[i], f"CD2Name mismatch at {i}"

        for seg in mesh2.Elements1D():
            assert seg.index >= 1
    finally:
        os.unlink(tmpfile)


def test_save_load_2d_spline_bc_names():
    """2D SplineGeometry with named boundaries: BC names survive save/load."""
    from netgen.geom2d import SplineGeometry
    from netgen.meshing import Mesh
    import tempfile, os

    geo = SplineGeometry()
    p = [geo.AppendPoint(x,y) for x,y in [(0,0),(1,0),(1,1),(0,1)]]
    names = ["bottom", "right", "top", "left"]
    geo.Append(["line", p[0], p[1]], leftdomain=1, rightdomain=0, bc=names[0])
    geo.Append(["line", p[1], p[2]], leftdomain=1, rightdomain=0, bc=names[1])
    geo.Append(["line", p[2], p[3]], leftdomain=1, rightdomain=0, bc=names[2])
    geo.Append(["line", p[3], p[0]], leftdomain=1, rightdomain=0, bc=names[3])
    mesh = geo.GenerateMesh(maxh=0.3)

    nseg_before = len(mesh.Elements1D())

    with tempfile.NamedTemporaryFile(suffix=".vol", delete=False) as f:
        tmpfile = f.name
    try:
        mesh.Save(tmpfile)
        mesh2 = Mesh()
        mesh2.Load(tmpfile)

        assert len(mesh2.Elements1D()) == nseg_before

        # Check that all BC names are present on edge descriptors after load
        loaded_ednames = set()
        for seg in mesh2.Elements1D():
            ed = mesh2.EdgeDescriptor(seg.index)
            loaded_ednames.add(ed.name)
        for nm in names:
            assert nm in loaded_ednames, f"BC name '{nm}' not found after load"
    finally:
        os.unlink(tmpfile)


def test_save_load_csg_segment_consistency():
    """CSG with two stacked boxes: segment count and edge descriptors survive save/load."""
    from netgen.csg import CSGeometry, OrthoBrick, Pnt
    from netgen.meshing import Mesh
    import tempfile, os

    geo = CSGeometry()
    geo.Add(OrthoBrick(Pnt(0,0,0), Pnt(1,1,0.5)).mat("lower"))
    geo.Add(OrthoBrick(Pnt(0,0,0.5), Pnt(1,1,1)).mat("upper"))
    mesh = geo.GenerateMesh(maxh=0.5)

    nseg_before = len(mesh.Elements1D())
    ned_before = mesh.GetNED()

    with tempfile.NamedTemporaryFile(suffix=".vol", delete=False) as f:
        tmpfile = f.name
    try:
        mesh.Save(tmpfile)
        mesh2 = Mesh()
        mesh2.Load(tmpfile)

        assert len(mesh2.Elements1D()) == nseg_before
        assert mesh2.GetNED() == ned_before
        for seg in mesh2.Elements1D():
            assert seg.index >= 1
    finally:
        os.unlink(tmpfile)


def test_medit_export_import_edsi():
    """Export to Medit format and reimport: segment count survives, index not preserved."""
    pytest.importorskip("netgen.occ")
    from netgen.occ import Box, Pnt, OCCGeometry
    from netgen import meshing
    import tempfile, os

    geo = OCCGeometry(Box(Pnt(0,0,0), Pnt(1,1,1)))
    mesh = geo.GenerateMesh(maxh=0.5)
    nseg_before = len(mesh.Elements1D())

    with tempfile.NamedTemporaryFile(suffix=".mesh", delete=False) as f:
        medit_file = f.name
    try:
        mesh.Export(medit_file, "Medit Format")
        mesh2 = meshing.ImportMesh(medit_file, "Medit Format")

        # Medit format preserves segments; EDs are now reconstructed on import
        assert len(mesh2.Elements1D()) == nseg_before
        assert mesh2.GetNED() > 0, "Medit import should reconstruct edge descriptors"
    finally:
        os.unlink(medit_file)


def test_backward_compat_load_without_edgedescriptors():
    """Loading a .vol file without edgedescriptors section triggers reconstruction."""
    pytest.importorskip("netgen.occ")
    from netgen.occ import Box, Pnt, OCCGeometry
    from netgen.meshing import Mesh
    import tempfile, os, re

    geo = OCCGeometry(Box(Pnt(0,0,0), Pnt(1,1,1)))
    mesh = geo.GenerateMesh(maxh=0.5)

    with tempfile.NamedTemporaryFile(suffix=".vol", delete=False) as f:
        tmpfile = f.name
    try:
        mesh.Save(tmpfile)

        # Strip the edgedescriptors section from the .vol file
        with open(tmpfile) as fh:
            content = fh.read()
        assert "edgedescriptors" in content

        # Remove the section: starts with edgedescriptors header, ends before next section
        content = re.sub(
            r"\n*edgedescriptors\n\d+\n(?:.*\n)*?(?=\n\w|\Z)",
            "\n",
            content,
        )
        assert "edgedescriptors" not in content

        with open(tmpfile, "w") as fh:
            fh.write(content)

        # Load the stripped file -- should trigger ReconstructEdgeDescriptors
        mesh2 = Mesh()
        mesh2.Load(tmpfile)

        assert mesh2.GetNED() > 0, "reconstruction should create edge descriptors"
        for seg in mesh2.Elements1D():
            assert seg.index >= 1, f"reconstructed segment has index={seg.index}"
    finally:
        os.unlink(tmpfile)


def test_setcd2name_syncs_edge_descriptor_after_load():
    """SetCD2Name updates EdgeDescriptor after load; persists through another save/load."""
    pytest.importorskip("netgen.occ")
    from netgen.occ import Box, Pnt, OCCGeometry
    from netgen.meshing import Mesh
    import tempfile, os

    geo = OCCGeometry(Box(Pnt(0,0,0), Pnt(1,1,1)))
    mesh = geo.GenerateMesh(maxh=0.5)

    with tempfile.NamedTemporaryFile(suffix=".vol", delete=False) as f:
        tmp1 = f.name
    with tempfile.NamedTemporaryFile(suffix=".vol", delete=False) as f:
        tmp2 = f.name
    try:
        mesh.Save(tmp1)
        mesh2 = Mesh()
        mesh2.Load(tmp1)

        # SetCD2Name is 1-based
        mesh2.SetCD2Name(1, "renamed_edge")
        assert mesh2.EdgeDescriptor(1).name == "renamed_edge"
        assert mesh2.GetCD2Name(0) == "renamed_edge"

        # Save again and reload -- name should persist
        mesh2.Save(tmp2)
        mesh3 = Mesh()
        mesh3.Load(tmp2)
        assert mesh3.EdgeDescriptor(1).name == "renamed_edge"
        assert mesh3.GetCD2Name(0) == "renamed_edge"
    finally:
        os.unlink(tmp1)
        os.unlink(tmp2)


def test_boundary_names_roundtrip_occ_3d():
    """OCC box with named edges: edge names available via GetCD2Name."""
    pytest.importorskip("netgen.occ")
    from netgen.occ import Box, Pnt, OCCGeometry
    box = Box(Pnt(0,0,0), Pnt(1,1,1))
    # OCC Box exposes 24 edge references (2 per geometric edge), name them all
    for i, e in enumerate(box.edges):
        e.name = f"ed_{i}"
    geo = OCCGeometry(box)
    mesh = geo.GenerateMesh(maxh=0.5)

    # After meshing there are 12 unique edges; each should have a name starting with "ed_"
    assert mesh.GetNED() == 12
    for i in range(mesh.GetNCD2Names()):
        name = mesh.GetCD2Name(i)
        assert name.startswith("ed_"), f"edge {i} has unexpected name: {name}"


def test_boundary_names_roundtrip_2d_spline():
    """2D SplineGeometry with named BCs: names accessible via GetBCName."""
    from netgen.geom2d import SplineGeometry
    geo = SplineGeometry()
    p = [geo.AppendPoint(x,y) for x,y in [(0,0),(1,0),(1,1),(0,1)]]
    geo.Append(["line", p[0], p[1]], leftdomain=1, rightdomain=0, bc="bottom")
    geo.Append(["line", p[1], p[2]], leftdomain=1, rightdomain=0, bc="right")
    geo.Append(["line", p[2], p[3]], leftdomain=1, rightdomain=0, bc="top")
    geo.Append(["line", p[3], p[0]], leftdomain=1, rightdomain=0, bc="left")
    mesh = geo.GenerateMesh(maxh=0.3)

    bc_names = set()
    for i in range(mesh.GetNED()):
        name = mesh.GetBCName(i)
        if name:
            bc_names.add(name)
    assert {"bottom", "right", "top", "left"}.issubset(bc_names)


def test_save_load_save_compare_occ():
    """OCC 3D: save -> load -> save produces identical .vol output."""
    pytest.importorskip("netgen.occ")
    from netgen.occ import Box, Pnt, OCCGeometry
    from netgen.meshing import Mesh
    geo = OCCGeometry(Box(Pnt(0,0,0), Pnt(1,1,1)))
    mesh = geo.GenerateMesh(maxh=0.5)

    with tempfile.TemporaryDirectory() as tmpdir:
        f1 = os.path.join(tmpdir, "mesh1.vol")
        f2 = os.path.join(tmpdir, "mesh2.vol")
        mesh.Save(f1)
        mesh2 = Mesh()
        mesh2.Load(f1)
        mesh2.Save(f2)

        with open(f1) as a, open(f2) as b:
            lines_a = a.readlines()
            lines_b = b.readlines()
        # exclude face_colours section
        for lines in (lines_a, lines_b):
            for i, line in enumerate(lines):
                if line.strip().startswith("face_colours"):
                    del lines[i:]
                    break
        assert lines_a == lines_b, "OCC save/load/save mismatch"


def test_save_load_save_compare_2d():
    """2D spline: save -> load -> save produces identical .vol output."""
    from netgen.geom2d import SplineGeometry
    from netgen.meshing import Mesh
    geo = SplineGeometry()
    p = [geo.AppendPoint(x, y) for x, y in [(0, 0), (1, 0), (1, 1), (0, 1)]]
    geo.Append(["line", p[0], p[1]], leftdomain=1, rightdomain=0, bc="bottom")
    geo.Append(["line", p[1], p[2]], leftdomain=1, rightdomain=0, bc="right")
    geo.Append(["line", p[2], p[3]], leftdomain=1, rightdomain=0, bc="top")
    geo.Append(["line", p[3], p[0]], leftdomain=1, rightdomain=0, bc="left")
    mesh = geo.GenerateMesh(maxh=0.3)

    with tempfile.TemporaryDirectory() as tmpdir:
        f1 = os.path.join(tmpdir, "mesh1.vol")
        f2 = os.path.join(tmpdir, "mesh2.vol")
        mesh.Save(f1)
        mesh2 = Mesh()
        mesh2.Load(f1)
        mesh2.Save(f2)

        with open(f1) as a, open(f2) as b:
            lines_a = a.readlines()
            lines_b = b.readlines()
        for lines in (lines_a, lines_b):
            for i, line in enumerate(lines):
                if line.strip().startswith("face_colours"):
                    del lines[i:]
                    break
        assert lines_a == lines_b, "2D save/load/save mismatch"


def test_save_load_roundtrip_ed_content_occ():
    """OCC: EdgeDescriptor content (edgenr, surfnr, name) survives save/load."""
    pytest.importorskip("netgen.occ")
    from netgen.occ import Box, Pnt, OCCGeometry
    from netgen.meshing import Mesh
    import tempfile, os

    geo = OCCGeometry(Box(Pnt(0,0,0), Pnt(1,1,1)))
    mesh = geo.GenerateMesh(maxh=0.5)

    orig_eds = []
    for i in range(mesh.GetNED()):
        ed = mesh.EdgeDescriptor(i+1)
        orig_eds.append((ed.edgenr, ed.surfnr, ed.name))

    with tempfile.NamedTemporaryFile(suffix=".vol", delete=False) as f:
        tmpfile = f.name
    try:
        mesh.Save(tmpfile)
        mesh2 = Mesh()
        mesh2.Load(tmpfile)

        assert mesh2.GetNED() == len(orig_eds), \
            f"NED mismatch: {mesh2.GetNED()} != {len(orig_eds)}"
        for i, (ednr, snr, nm) in enumerate(orig_eds):
            ed2 = mesh2.EdgeDescriptor(i+1)
            assert ed2.edgenr == ednr, f"ED[{i}] edgenr: {ed2.edgenr} != {ednr}"
            assert ed2.surfnr == snr, f"ED[{i}] surfnr: {ed2.surfnr} != {snr}"
            assert ed2.name == nm, f"ED[{i}] name: {ed2.name} != {nm}"
    finally:
        os.unlink(tmpfile)


def test_fdindex_occ_3d():
    """OCC 3D: every ED has fdindex that is a valid face descriptor index."""
    pytest.importorskip("netgen.occ")
    from netgen.occ import Box, Pnt, OCCGeometry
    geo = OCCGeometry(Box(Pnt(0,0,0), Pnt(1,1,1)))
    mesh = geo.GenerateMesh(maxh=0.5)
    nfd = mesh.GetNFaceDescriptors()
    for i in range(mesh.GetNED()):
        ed = mesh.EdgeDescriptor(i+1)
        assert ed.fdindex > 0, f"ED[{i}] fdindex={ed.fdindex} should be > 0"
        assert ed.fdindex <= nfd, f"ED[{i}] fdindex={ed.fdindex} exceeds nfd={nfd}"


def test_fdindex_csg_3d():
    """CSG 3D: every ED has fdindex > 0 after meshing (FindEdges sets it)."""
    from netgen.csg import CSGeometry, OrthoBrick, Pnt
    geo = CSGeometry()
    geo.Add(OrthoBrick(Pnt(0,0,0), Pnt(1,1,1)))
    mesh = geo.GenerateMesh(maxh=0.5)
    for i in range(mesh.GetNED()):
        ed = mesh.EdgeDescriptor(i+1)
        assert ed.fdindex > 0, f"ED[{i}] fdindex={ed.fdindex} should be > 0"


def test_fdindex_2d():
    """2D spline: every ED has fdindex > 0."""
    from netgen.geom2d import SplineGeometry
    geo = SplineGeometry()
    p1, p2, p3, p4 = [geo.AppendPoint(x, y) for x, y in [(0,0),(1,0),(1,1),(0,1)]]
    geo.Append(["line", p1, p2], leftdomain=1, rightdomain=0)
    geo.Append(["line", p2, p3], leftdomain=1, rightdomain=0)
    geo.Append(["line", p3, p4], leftdomain=1, rightdomain=0)
    geo.Append(["line", p4, p1], leftdomain=1, rightdomain=0)
    mesh = geo.GenerateMesh(maxh=0.3)
    nfd = mesh.GetNFaceDescriptors()
    for i in range(mesh.GetNED()):
        ed = mesh.EdgeDescriptor(i+1)
        assert ed.fdindex > 0, f"ED[{i}] fdindex={ed.fdindex} should be > 0"
        # In 2D, fdindex stores the BC number (spline.bc), not a FaceDescriptor index,
        # so it can exceed nfd. Only check upper bound for 3D.
    for seg in mesh.Elements1D():
        ed = mesh.EdgeDescriptor(seg.index)
        assert ed.fdindex > 0, \
            f"ED[{seg.index}] fdindex={ed.fdindex} should be > 0"


def test_save_load_second_order_mesh():
    """Save/load a second-order OCC box mesh preserves segments, indices, and ED edgenr values."""
    import tempfile, os
    pytest.importorskip("netgen.occ")
    from netgen.occ import Box, Pnt, OCCGeometry
    from netgen.meshing import Mesh

    geo = OCCGeometry(Box(Pnt(0,0,0), Pnt(1,1,1)))
    mesh = geo.GenerateMesh(maxh=0.5)
    mesh.SecondOrder()

    nseg_before = len(mesh.Elements1D())
    ned_before = mesh.GetNED()
    edgenrs_before = [mesh.EdgeDescriptor(i+1).edgenr for i in range(ned_before)]

    with tempfile.NamedTemporaryFile(suffix=".vol", delete=False) as f:
        tmpfile = f.name
    try:
        mesh.Save(tmpfile)
        mesh2 = Mesh()
        mesh2.Load(tmpfile)
        assert len(mesh2.Elements1D()) == nseg_before
        for seg in mesh2.Elements1D():
            assert seg.index >= 1, f"segment has index={seg.index}"
        assert mesh2.GetNED() == ned_before
        for i in range(ned_before):
            assert mesh2.EdgeDescriptor(i+1).edgenr == edgenrs_before[i], \
                f"ED[{i}] edgenr mismatch after second-order save/load"
    finally:
        os.unlink(tmpfile)


def test_pickle_roundtrip_full_ed_content():
    """Pickle roundtrip preserves all ED properties: edgenr, surfnr, name, singedge_left/right, domin, domout."""
    import pickle
    pytest.importorskip("netgen.occ")
    from netgen.occ import Box, Pnt, OCCGeometry

    geo = OCCGeometry(Box(Pnt(0,0,0), Pnt(1,1,1)))
    mesh = geo.GenerateMesh(maxh=0.5)

    orig_eds = []
    for i in range(mesh.GetNED()):
        ed = mesh.EdgeDescriptor(i+1)
        orig_eds.append({
            "edgenr": ed.edgenr,
            "surfnr": ed.surfnr,
            "name": ed.name,
            "singedge_left": ed.singedge_left,
            "singedge_right": ed.singedge_right,
            "domin": ed.domin,
            "domout": ed.domout,
        })

    data = pickle.dumps(mesh)
    mesh2 = pickle.loads(data)

    assert mesh2.GetNED() == len(orig_eds), \
        f"NED mismatch: {mesh2.GetNED()} != {len(orig_eds)}"
    for i, props in enumerate(orig_eds):
        ed2 = mesh2.EdgeDescriptor(i+1)
        assert ed2.edgenr == props["edgenr"], f"ED[{i}] edgenr mismatch"
        assert ed2.surfnr == props["surfnr"], f"ED[{i}] surfnr mismatch"
        assert ed2.name == props["name"], f"ED[{i}] name mismatch"
        assert ed2.singedge_left == pytest.approx(props["singedge_left"]), \
            f"ED[{i}] singedge_left mismatch"
        assert ed2.singedge_right == pytest.approx(props["singedge_right"]), \
            f"ED[{i}] singedge_right mismatch"
        assert ed2.domin == props["domin"], f"ED[{i}] domin mismatch"
        assert ed2.domout == props["domout"], f"ED[{i}] domout mismatch"


def test_save_load_boundarylayer_mesh():
    """Save/load a boundary-layer mesh preserves segment count, ED count, and ED content."""
    pytest.importorskip("netgen.occ")
    from netgen.occ import Box, Pnt, OCCGeometry
    from netgen.meshing import Mesh, BoundaryLayerParameters

    geo = OCCGeometry(Box(Pnt(0,0,0), Pnt(1,1,1)))
    blparams = [BoundaryLayerParameters(".*", [0.01, 0.01], "layer", outside=True)]
    mesh = geo.GenerateMesh(maxh=0.3, boundary_layers=blparams)

    nseg_before = len(mesh.Elements1D())
    ned_before = mesh.GetNED()
    edgenrs_before = [mesh.EdgeDescriptor(i+1).edgenr for i in range(ned_before)]

    with tempfile.NamedTemporaryFile(suffix=".vol", delete=False) as f:
        tmpfile = f.name
    try:
        mesh.Save(tmpfile)
        mesh2 = Mesh()
        mesh2.Load(tmpfile)
        assert len(mesh2.Elements1D()) == nseg_before
        assert mesh2.GetNED() == ned_before
        for i in range(mesh2.GetNED()):
            ed = mesh2.EdgeDescriptor(i+1)
            assert ed.edgenr == edgenrs_before[i], f"ED[{i}] edgenr mismatch after save/load"
        for seg in mesh2.Elements1D():
            assert seg.index >= 1, f"segment has index={seg.index} after save/load"
    finally:
        os.unlink(tmpfile)


def test_save_load_glue_compound():
    """Save/load a glued two-box compound preserves segment and ED counts."""
    pytest.importorskip("netgen.occ")
    from netgen.occ import Box, Pnt, OCCGeometry, Glue
    from netgen.meshing import Mesh

    box1 = Box(Pnt(0,0,0), Pnt(1,1,1))
    box2 = Box(Pnt(1,0,0), Pnt(2,1,1))
    geo = OCCGeometry(Glue([box1, box2]))
    mesh = geo.GenerateMesh(maxh=0.5)

    nseg_before = len(mesh.Elements1D())
    ned_before = mesh.GetNED()

    with tempfile.NamedTemporaryFile(suffix=".vol", delete=False) as f:
        tmpfile = f.name
    try:
        mesh.Save(tmpfile)
        mesh2 = Mesh()
        mesh2.Load(tmpfile)
        assert len(mesh2.Elements1D()) == nseg_before
        assert mesh2.GetNED() == ned_before
    finally:
        os.unlink(tmpfile)
