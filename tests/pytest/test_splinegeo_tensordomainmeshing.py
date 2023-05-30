from netgen.geom2d import *

def test_tensordomainmeshing():
    geo = SplineGeometry()
    w = 10
    h = 0.01

    p = [ (0, 0), (w, 0), (w, h), (0, h) ]
    p = [geo.AppendPoint(*px) for px in p]

    l0 = geo.Append ( ["line", p[0], p[1]], leftdomain=1, rightdomain=0 )
    l1 = geo.Append ( ["line", p[1], p[2]], leftdomain=1, rightdomain=0)
    geo.Append ( ["line", p[3], p[2]], leftdomain=0, rightdomain=1, copy=l0 )
    geo.Append ( ["line", p[0], p[3]], leftdomain=0, rightdomain=1, copy=l1 )

    geo._SetDomainTensorMeshing(1, True)

    mesh = geo.GenerateMesh(maxh=1)

    for el in mesh.Elements2D():
        print(el.vertices)
        assert len(el.vertices) == 4
