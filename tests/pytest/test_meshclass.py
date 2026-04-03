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
