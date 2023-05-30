import pyngcore
import netgen

from meshes import unit_mesh_3d

def test_element_arrays(unit_mesh_3d):
    mesh = unit_mesh_3d
    el0 = mesh.Elements0D()
    el1 = mesh.Elements1D()
    el2 = mesh.Elements2D()
    el3 = mesh.Elements3D()
    p = mesh.Points()

    assert len(el2) > 0
    assert len(el3) > 0
    assert len(p) > 0

    for el in el2:
        assert len(el.vertices) == 3

    for el in el3:
        assert len(el.vertices) == 4
