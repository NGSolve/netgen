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
    from netgen.meshing import PointId, ElementType
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
        assert typ == ElementType.TET.value
        assert el.type == ElementType.TET
        assert len(el.vertices) == 4

def test_element_type_tables():
    from netgen.csg import unit_cube
    from netgen.meshing import ElementNP, ElementNV, ElementNEdges, ElementNFaces, ElementType
    import numpy
    mesh = unit_cube.GenerateMesh()
    for els in (mesh.Elements3D(), mesh.Elements2D()):
        a = els.NumPy()
        nps = ElementNP[a["type"]]
        nvs = ElementNV[a["type"]]
        assert all(nps == numpy.array([len(el.points) for el in els]))
        assert all(nvs == numpy.array([len(el.vertices) for el in els]))
    types3d = mesh.Elements3D().NumPy()["type"]
    assert all(ElementNEdges[types3d] == 6) and all(ElementNFaces[types3d] == 4)
    assert not ElementNP.flags.writeable
    assert ElementNP[ElementType.TET.value] == 4 and ElementNP[ElementType.HEX20.value] == 20
    assert all(el.type == ElementType.SEGMENT for el in mesh.Elements1D())

if __name__ == "__main__":
    test_array_numpy()
    test_mesh_elements_numpy_array_access()
    test_element_type_tables()
