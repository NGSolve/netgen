
import pytest
from netgen.csg import *

def GetNSurfaceElements(mesh, boundary):
    nse_in_layer = 0
    for el in mesh.Elements2D():
        print(el.index)
        if mesh.GetBCName(el.index-1) == boundary:
            nse_in_layer += 1
    return nse_in_layer

@pytest.mark.parametrize("outside", [True, False])
def test_boundarylayer(outside, capfd):
    mesh = unit_cube.GenerateMesh(maxh=0.3)
    ne_before = mesh.ne
    layer_surfacenames = ["right", "top"]
    mesh.BoundaryLayer("|".join(layer_surfacenames), [0.01, 0.02], "layer", outside=outside, grow_edges=True)

    should_ne = ne_before + 2 * sum([GetNSurfaceElements(mesh, surf) for surf in layer_surfacenames])
    assert mesh.ne == should_ne
    capture = capfd.readouterr()
    assert not "elements are not matching" in capture.out

    for side in ["front"]:
        mesh.BoundaryLayer(side, [0.001], "layer", outside=outside, grow_edges=True)
        should_ne += GetNSurfaceElements(mesh, side)
        assert mesh.ne == should_ne
        capture = capfd.readouterr()
        assert not "elements are not matching" in capture.out
    
