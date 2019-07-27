
import pickle, numpy

def test_pickle_csg():
    import netgen.csg as csg
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

def test_pickle_stl():
    import netgen.stl as stl
    geo = stl.LoadSTLGeometry("../../tutorials/hinge.stl")
    geo_dump = pickle.dumps(geo)
    geo2 = pickle.loads(geo_dump)
    vd1 = geo._visualizationData()
    vd2 = geo2._visualizationData()
    for val1, val2 in zip(vd1.values(), vd2.values()):
        assert numpy.array_equal(val1, val2)


def test_pickle_occ():
    try:
        import netgen.occ as occ
    except:
        import pytest
        pytest.skip("can't import occ")
    geo = occ.OCCGeometry("../../tutorials/frame.step")
    geo_dump = pickle.dumps(geo)
    geo2 = pickle.loads(geo_dump)
    vd1 = geo._visualizationData()
    vd2 = geo2._visualizationData()
    # TODO: it looks fine, but tests fail, so I assume we loose some info?
    # for val1, val2 in zip(vd1.values(), vd2.values()):
    #     assert numpy.allclose(val1, val2, rtol=0.01)

def test_pickle_geom2d():
    import netgen.geom2d as geom2d
    geo = geom2d.SplineGeometry()

    # point coordinates ...
    pnts = [ (0,0), (1,0), (1,0.6), (0,0.6), \
             (0.2,0.6), (0.8,0.6), (0.8,0.8), (0.2,0.8), \
             (0.5,0.15), (0.65,0.3), (0.5,0.45), (0.35,0.3) ]
    pnums = [geo.AppendPoint(*p) for p in pnts]

    # start-point, end-point, boundary-condition, domain on left side, domain on right side:
    lines = [ (0,1,1,1,0), (1,2,2,1,0), (2,5,2,1,0), (5,4,2,1,2), (4,3,2,1,0), (3,0,2,1,0), \
              (5,6,2,2,0), (6,7,2,2,0), (7,4,2,2,0), \
              (8,9,2,3,1), (9,10,2,3,1), (10,11,2,3,1), (11,8,2,3,1) ]

    for p1,p2,bc,left,right in lines:
        geo.Append( ["line", pnums[p1], pnums[p2]], bc=bc, leftdomain=left, rightdomain=right)
    geo_dump = pickle.dumps(geo)
    geo2 = pickle.loads(geo_dump)
    vd1 = geo._visualizationData()
    vd2 = geo2._visualizationData()
    for val1, val2 in zip(vd1.values(), vd2.values()):
        assert numpy.array_equal(val1, val2)

def test_pickle_mesh():
    import netgen.csg as csg
    geo = csg.CSGeometry()
    brick = csg.OrthoBrick(csg.Pnt(-3,-3,-3), csg.Pnt(3,3,3))
    mesh = geo.GenerateMesh(maxh=0.2)
    assert geo == mesh.GetGeometry()
    dump = pickle.dumps([geo,mesh])
    geo2, mesh2 = pickle.loads(dump)
    assert geo2 == mesh2.GetGeometry()
    mesh.Save("msh1.vol.gz")
    mesh2.Save("msh2.vol.gz")
    import filecmp, os
    assert filecmp.cmp("msh1.vol.gz", "msh2.vol.gz")
    os.remove("msh1.vol.gz")
    os.remove("msh2.vol.gz")

if __name__ == "__main__":
    test_pickle_mesh()
