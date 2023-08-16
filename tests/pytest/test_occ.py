import pytest
import math
from pytest import approx
from pytest_check import check

def check_volume(shape, volume, dim=3):
    from netgen.occ import OCCGeometry
    geo = OCCGeometry(shape, dim=dim)

    m = geo.GenerateMesh()
    assert len(m.Elements2D()) > 0
    assert len(m.Elements1D()) > 0
    if dim==3:
        assert len(m.Elements3D()) > 0

    ngs = pytest.importorskip("ngsolve")
    mesh = ngs.Mesh(m)
    mesh.Curve(5)
    with check: assert ngs.Integrate(1.0, mesh) == approx(volume)

def test_rect_with_two_holes():
    occ = pytest.importorskip("netgen.occ")

    face = occ.WorkPlane().Rectangle(7,4) \
    .Circle(2,2,1).Reverse() \
    .Circle(5,2,1).Reverse().Face()
    check_volume(face, 7*4-2*math.pi, 2)

def test_unit_square():
    occ = pytest.importorskip("netgen.occ")
    check_volume(occ.unit_square.shape, 1, dim=2)

def test_box_and_cyl():
    occ = pytest.importorskip("netgen.occ")
    box = occ.Box(occ.Pnt(0,0,0), occ.Pnt(1,1,1))
    check_volume(box, 1)
    r = 0.3
    h = 0.5
    vcyl = r*r*math.pi*h
    cyl = occ.Cylinder(occ.Pnt(1,0.5,0.5), occ.X, r=r, h=h)
    check_volume(cyl, vcyl)
    fused = box+cyl
    check_volume(fused, 1+vcyl)

def test_internal_face():
    occ = pytest.importorskip("netgen.occ")
    box = occ.Box((0,0,0), (3, 1, 10))

    face = occ.WorkPlane(occ.Axes((1.5,0,0), occ.X, occ.Y)).Rectangle(1, 6).Face()

    shape = occ.Glue([box, face])
    geo = occ.OCCGeometry(shape)
    mesh = geo.GenerateMesh(maxh=0.5)
    assert any(mesh.Elements2D().NumPy()['index'] == 8)

