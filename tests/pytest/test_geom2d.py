from netgen.geom2d import SplineGeometry


def test_leftdom_equals_rightdom():
    geo = SplineGeometry()
    pnts = [(0,0), (1,0), (2,0), (2,1), (1,1), (0,1)]
    gp = [geo.AppendPoint(*p) for p in pnts]
    lines = [(0,1,0), (1,2,0), (2,3,0), (3,4,0), (4,5,0), (5,0,0), (1,4,1)]
    for p1, p2, rd in lines:
        geo.Append(["line", p1, p2], leftdomain=1, rightdomain=rd)

    mesh = geo.GenerateMesh()

