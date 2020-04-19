
import pytest
from netgen.csg import *

@pytest.mark.parametrize("outside", [True, False])
def test_boundarylayer(outside):
    mesh = unit_cube.GenerateMesh(maxh=0.3)
    ne_before = mesh.ne
    nse_in_layer = 0
    layer_surfacenames = ["right", "top"]
    for el in mesh.Elements2D():
        if mesh.GetBCName(el.index-1) in layer_surfacenames:
            nse_in_layer += 1
    mesh.BoundaryLayer("|".join(layer_surfacenames), [0.01, 0.02], "layer", outside=outside, grow_edges=True)
    assert mesh.ne == ne_before + 2 * nse_in_layer

    
