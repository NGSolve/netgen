import pyngcore
import netgen
import pytest
import tempfile

from meshes import unit_mesh_3d


def test_element_arrays(unit_mesh_3d):
    mesh = unit_mesh_3d
    el1 = mesh.Elements1D()
    el2 = mesh.Elements2D()
    el3 = mesh.Elements3D()
    p = mesh.Points()

    assert len(el1) > 0
    assert len(el2) > 0
    assert len(el3) > 0
    assert len(p) > 0

    for el in el2:
        assert len(el.vertices) == 3

    for el in el3:
        assert len(el.vertices) == 4


def test_copy_mesh():
    pytest.importorskip("netgen.occ")
    import netgen.occ as occ

    box1 = occ.Box((0, 0, 0), (1, 1, 1))
    box2 = occ.Box((1, 0, 0), (2, 1, 1))
    box1.faces.name = "bnd1"
    box1.name = "mat1"
    box2.faces.name = "bnd2"
    box2.name = "mat1"

    geo = occ.OCCGeometry(occ.Glue([box1, box2]))
    m3d = geo.GenerateMesh(maxh=99)

    plane1 = occ.WorkPlane(occ.Axes((0, 0, 0), occ.X, occ.Y)).Rectangle(1, 1).Face()
    plane1.name = "mat1"
    plane1.edges.name = "bnd1"

    plane2 = occ.WorkPlane(occ.Axes((0, 0, 0), occ.X, occ.Y)).Rectangle(2, 2).Face()
    plane2.name = "mat2"
    plane2.edges.name = "bnd2"

    geo2 = occ.OCCGeometry(occ.Glue([plane1, plane2]), dim=2)
    m2d = geo2.GenerateMesh(maxh=99)

    for mesh in [m2d, m3d]:
        copy = mesh.Copy()

        assert copy.dim == mesh.dim
        assert len(copy.Elements0D()) == len(mesh.Elements0D())
        assert len(copy.Elements1D()) == len(mesh.Elements1D())
        assert len(copy.Elements2D()) == len(mesh.Elements2D())
        assert len(copy.Elements3D()) == len(mesh.Elements3D())
        assert copy.GetNDomains() == mesh.GetNDomains()
        assert copy.GetNFaceDescriptors() == mesh.GetNFaceDescriptors()
        for dim in range(1, mesh.dim + 1):
            assert copy.GetRegionNames(dim) == mesh.GetRegionNames(dim)
        assert copy.GetIdentifications() == mesh.GetIdentifications()


def _brute_incidence(mesh, d1, d2):
    import itertools
    from collections import defaultdict
    els = {3: mesh.Elements3D, 2: mesh.Elements2D, 1: mesh.Elements1D, 0: mesh.Elements0D}
    hi, lo = max(d1, d2), min(d1, d2)
    # simplicial meshes: every vertex subset of an element is a sub-entity
    sub = defaultdict(set)
    for el in els[hi]():
        vs = [v.nr for v in el.vertices]
        for k in range(1, len(vs) + 1):
            for c in itertools.combinations(vs, k):
                sub[frozenset(c)].add(el.index - 1)
    inc = defaultdict(set)
    for el in els[lo]():
        for r in sub.get(frozenset(v.nr for v in el.vertices), ()):
            if d1 == hi: inc[r].add(el.index - 1)
            else: inc[el.index - 1].add(r)
    return inc

def _check_incidence(mesh):
    for d1 in range(mesh.dim + 1):
        for d2 in range(mesh.dim + 1):
            if d1 == d2: continue
            got = mesh.RegionIncidence(d1, d2)
            exp = _brute_incidence(mesh, d1, d2)
            assert all(r < len(got) for r in exp)
            for r, row in enumerate(got):
                assert set(row) == exp.get(r, set()), (d1, d2, r)

def test_region_incidence_3d():
    import netgen.occ as occ
    box = occ.Box((0,0,0), (1,1,1))
    box.edges.Min(occ.Z).name = "bottom"
    box.vertices.Min(occ.X+occ.Y+occ.Z).name = "corner"
    ball = occ.Sphere((0.5,0.5,0.5), 0.25)
    mesh = occ.OCCGeometry(occ.Glue([box - ball, ball])).GenerateMesh(maxh=0.3)
    _check_incidence(mesh)
    # two boxes touching along an edge only
    b2 = occ.Box((1,1,0), (2,2,1))
    _check_incidence(occ.OCCGeometry(occ.Glue([occ.Box((0,0,0),(1,1,1)), b2])).GenerateMesh(maxh=0.5))
    with pytest.raises(Exception):
        mesh.RegionIncidence(3, 3)

def test_region_incidence_2d():
    import netgen.occ as occ
    rect = occ.Rectangle(2, 1).Face()
    rect.vertices.Min(occ.X+occ.Y).name = "pin"
    circ = occ.Circle((1, 0.5), 0.2).Face()
    geo = occ.OCCGeometry(occ.Glue([rect - circ, circ]), dim=2)
    _check_incidence(geo.GenerateMesh(maxh=0.2))
    mesh = geo.GenerateMesh(maxh=0.2, quad_dominated=True)
    assert any(len(el.vertices) == 4 for el in mesh.Elements2D())
    _check_incidence(mesh)

def test_region_incidence_1d():
    import netgen.meshing as nm
    m = nm.Mesh(dim=1)
    pids = [m.Add(nm.MeshPoint(nm.Pnt(i/10, 0, 0))) for i in range(11)]
    for i in range(10):
        m.Add(nm.Element1D([pids[i], pids[i+1]], index=1 if i < 5 else 2))
    for p, idx in ((0, 1), (10, 2), (5, 3)):
        m.Add(nm.Element0D(pids[p], index=idx))
    assert m.RegionIncidence(0, 1) == [[0], [1], [0, 1]]
    _check_incidence(m)
    # no stale results after changing the mesh
    m.Add(nm.Element0D(pids[5], index=4))
    assert m.RegionIncidence(0, 1) == [[0], [1], [0, 1], [0, 1]]
    for el in m.Elements1D():
        el.index = 1
    assert m.RegionIncidence(0, 1) == [[0], [0], [0], [0]]
    _check_incidence(m)
