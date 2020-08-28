from netgen.geom2d import *
import pytest
import math
from pytest import approx

def test_two_circles():
    c1 = Circle(center=(0,0), radius=1)
    c2 = c1.Rotate(deg=45)
    s = c1*c2
    geo = CSG2d()
    geo.Add(s)
    m = geo.GenerateMesh()
    assert m.ne > 0

    ngs = pytest.importorskip("ngsolve")
    mesh = ngs.Mesh(m)
    mesh.Curve(5)
    assert ngs.Integrate(1.0, mesh) == approx(math.pi)

def test_two_edge():
    s = Solid2d( [(-1,0), cp(0,1), (1,0), cp(0,2)] )
    geo = CSG2d()
    geo.Add(s)
    m = geo.GenerateMesh()
    assert m.ne > 0

    ngs = pytest.importorskip("ngsolve")
    g = geo.GenerateSplineGeometry()
    ngs.Draw(g)

if __name__ == "__main__":
    test_two_circles()
    test_two_edge()
