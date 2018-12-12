
import netgen.csg as csg
import pickle, numpy

def test_pickle_csg():
    geo = csg.CSGeometry()
    geo.Add(csg.Sphere(csg.Pnt(0,0,0), 2).bc("sphere"))
    brick = csg.OrthoBrick(csg.Pnt(-3,-3,-3), csg.Pnt(3,3,3))
    geo.Add(csg.Cylinder(csg.Pnt(0,0,0), csg.Pnt(1,0,0), 0.5) * brick)
    geo.Add(csg.Ellipsoid(csg.Pnt(0,0,0), csg.Vec(1,0,0), csg.Vec(0,1,0), csg.Vec(0,0,0.5)))
    geo.Add(csg.Cone(csg.Pnt(0,0,0), csg.Pnt(3,0,0), 1, 0.5) * brick)
    geo.Add(csg.EllipticCone(csg.Pnt(0,0,0), csg.Vec(2,0,0), csg.Vec(0,1,0), 3, 0.5) * brick)
    geo.Add(csg.Torus(csg.Pnt(0,0,0), csg.Vec(0,1,0), 0.3, 0.05))
    pts2d = [[1,1], [1,-1], [-1,-1], [-1,1]]
    segs = [[0,1], [1,2], [2,3], [3,0]]
    curve = csg.SplineCurve2d()
    pnrs = [curve.AddPoint(*p) for p in pts2d]
    for s in segs:
            curve.AddSegment(pnrs[s[0]], pnrs[s[1]])
    geo.Add(csg.Revolution(csg.Pnt(0,0,0), csg.Pnt(1,0,0), curve))
    path = csg.SplineCurve3d()
    pnts = [(0,0,0), (2,0,0), (2,2,0)]
    segs = [(0,1,2)]
    for pnt in pnts:
        path.AddPoint (*pnt)

    for seg in segs:
        path.AddSegment (*seg)
    geo.Add(csg.Extrusion(path, curve, csg.Vec(0,0,1)))

    geo_dump = pickle.dumps(geo)
    geo2 = pickle.loads(geo_dump)
    vd1 = geo._visualizationData()
    vd2 = geo2._visualizationData()
    for val1, val2 in zip(vd1.values(), vd2.values()):
        assert numpy.array_equal(val1, val2)

if __name__ == "__main__":
    test_pickle_csg()
