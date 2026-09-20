import netgen
from pyngcore import *
from numpy import sort, array

def test_array_numpy():
    a = Array_I_S(5)
    a[:] = 0
    a[3:] = 2
    assert(sum(a) == 4)
    a[1] = 5
    b = sort(a)
    assert(all(b == array([0,0,2,2,5])))
    assert(all(a == array([0,5,0,2,2])))
    a.NumPy().sort()
    assert(all(a == array([0,0,2,2,5])))

def test_mesh_elements_numpy_array_access():
    from netgen.csg import unit_cube
    from netgen.meshing import PointId
    mesh = unit_cube.GenerateMesh()
    np_els = mesh.Elements3D().NumPy()
    vol_nodes = np_els["nodes"]
    indices = np_els["index"]
    types = np_els["type"]
    for nodes, el, index, typ in zip(vol_nodes, mesh.Elements3D(), indices, types):
        for n1, n2 in zip(nodes, el.vertices):
            assert n1 == n2
        for n in nodes[len(el.vertices):]:
            assert n == PointId.base - 1  # unused slot holds PointIndex::INVALID
        assert el.index == index
        assert typ == 20  # TET
        assert len(el.vertices) == 4

if __name__ == "__main__":
    test_array_numpy()
    test_mesh_elements_numpy_array_access()
