from netgen.csg import *

def test_2_polyhedra():
    geo = CSGeometry()
    first = Polyhedron([(0,0,0), (0,1,0), (3,1,0), (3,0,0),
                        (0,1,1), (3,1,1), (3,0,1), (0,0,1)],
                        [(0,1,2,3), (1,4,5,2), (2,5,6,3), (3,6,7,0),
                         (0,7,4,1), (7,6,5,4)])
    # TODO: height = 0.1 not working!
    height = 0.3
    second = Polyhedron([(0,0,1), (0,1,1), (3,1,1), (3,0,1),
                         (0,1,1+height), (3,1,1+height),
                         (3,0,1+height), (0,0,1+height)],
                        [(0,1,2,3), (1,4,5,2), (2,5,6,3), (3,6,7,0),
                         (0,7,4,1), (7,6,5,4)])

    geo.Add(first)
    geo.Add(second)
    mesh = geo.GenerateMesh()
    return mesh


if __name__ == "__main__":
    from ngsolve import Mesh, Draw
    mesh = Mesh(test_2_polyhedra())
    Draw(mesh)
