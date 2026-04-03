import pytest

@pytest.fixture
def unit_mesh_2d():
    import netgen.geom2d as g2d
    return g2d.unit_square.GenerateMesh(maxh=0.2)

@pytest.fixture
def unit_mesh_3d():
    import netgen.csg as csg
    return csg.unit_cube.GenerateMesh(maxh=0.2)
