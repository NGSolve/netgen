from netgen.geom2d import *
import pytest
import math
from pytest import approx
from pytest_check import check


def check_area(geo, area):
    if isinstance(geo, Solid2d):
        g = CSG2d()
        g.Add(geo)
        geo = g

    m = geo.GenerateMesh(maxh=0.2)
    ngs = pytest.importorskip("ngsolve")
    mesh = ngs.Mesh(m)
    mesh.Curve(5)
    with check: assert ngs.Integrate(1.0, mesh) == approx(area)

def test_two_circles():
    c1 = Circle(center=(0,0), radius=1)
    c2 = c1.Rotate(45)
    s = c1*c2
    geo = CSG2d()
    geo.Add(s)
    m = geo.GenerateMesh()
    assert len(m.Elements2D()) > 0

    ngs = pytest.importorskip("ngsolve")
    mesh = ngs.Mesh(m)
    mesh.Curve(5)
    assert ngs.Integrate(1.0, mesh) == approx(math.pi)
    ngs.Draw(mesh)

def test_two_edge():
    s = Solid2d( [(-1,0), cp(0,1), (1,0), cp(0,2)] )
    geo = CSG2d()
    geo.Add(s)
    m = geo.GenerateMesh()
    assert len(m.Elements2D()) > 0

    ngs = pytest.importorskip("ngsolve")
    g = geo.GenerateSplineGeometry()
    ngs.Draw(g)
    mesh = ngs.Mesh(m)
    mesh.Curve(5)
    ngs.Draw(mesh)

def test_trig_and_circle():
    g = CSG2d()

    trig = Solid2d( [(0,0), (1,1), (-1,1) ] ).BC("diamond")
    circle = Circle( center=(0,0.101), radius=0.1).BC("circle") # TODO: Failing with center=(0,0.1)

    d = trig-circle
    g.Add(d)
    g.Add(circle)

    m = g.GenerateMesh(maxh=0.1)
    assert len(m.Elements2D()) > 0

    ngs = pytest.importorskip("ngsolve")
    geo = g.GenerateSplineGeometry()
    ngs.Draw(geo)

    mesh = ngs.Mesh(m)
    mesh.Curve(3)
    ngs.Draw(mesh)


def test_circle_plus_rect():
    circle = Circle( center=(0,0), radius=1 )
    rect = Rectangle( pmin=(-0.5,0.0), pmax=(0.5,0.5) )

    geo = CSG2d()
    geo.Add(circle+rect)
    m = geo.GenerateMesh(maxh=0.2)


    ngs = pytest.importorskip("ngsolve")
    mesh = ngs.Mesh(m)
    mesh.Curve(5)
    assert ngs.Integrate(1.0, mesh) == approx(math.pi)

def test_circle_plus_rect1():
    circle = Circle( center=(0,0), radius=1 )
    rect = Rectangle( pmin=(-0.5,-0.5), pmax=(0.5,0.5) )

    geo = CSG2d()
    geo.Add(circle+rect)
    m = geo.GenerateMesh(maxh=0.2)


    ngs = pytest.importorskip("ngsolve")
    mesh = ngs.Mesh(m)
    mesh.Curve(5)
    assert ngs.Integrate(1.0, mesh) == approx(math.pi)

def test_circle_and_rect():
    c = Circle(center=(0,0),radius=1)
    r = Rectangle((0,0),(1,1))

    pi = math.pi
    check_area(c-r, 3/4*pi)
    check_area(c*r, 1/4*pi)
    check_area(c+r, 3/4*pi+1)
    check_area(r*c, 1/4*pi)
    check_area(r+c, 3/4*pi+1)
    check_area(r-c, 1-1/4*pi)


if __name__ == "__main__":
    test_two_circles()
    test_two_edge()
    test_trig_and_circle()
