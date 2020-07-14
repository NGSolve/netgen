
import pytest
from netgen.csg import *

def GetNSurfaceElements(mesh, boundaries):
    nse_in_layer = 0
    for el in mesh.Elements2D():
        if mesh.GetBCName(el.index-1) in boundaries:
            nse_in_layer += 1
    return nse_in_layer

@pytest.mark.parametrize("outside", [True, False])
def test_boundarylayer(outside, capfd):
    mesh = unit_cube.GenerateMesh(maxh=0.3)
    ne_before = mesh.ne
    layer_surfacenames = ["right", "top", "left", "back", "bottom"]
    mesh.BoundaryLayer("|".join(layer_surfacenames), [0.01, 0.02], "layer", outside=outside, grow_edges=True)

    should_ne = ne_before + 2 * GetNSurfaceElements(mesh, layer_surfacenames)
    assert mesh.ne == should_ne
    capture = capfd.readouterr()
    assert not "elements are not matching" in capture.out

    for side in ["front"]:
        mesh.BoundaryLayer(side, [0.001, 0.002], "layer", outside=outside, grow_edges=True)
        should_ne += 2 * GetNSurfaceElements(mesh, [side])
        assert mesh.ne == should_ne
        capture = capfd.readouterr()
        assert not "elements are not matching" in capture.out

@pytest.mark.parametrize("outside", [True, False])
@pytest.mark.parametrize("version", [1]) # version 2 not working yet
def test_boundarylayer2(outside, version, capfd):
    geo = CSGeometry()
    top = Plane(Pnt(0,0,0.5), Vec(0,0,1))
    bot = Plane(Pnt(0,0,0.1), Vec(0,0,-1))
    bigpart = OrthoBrick(Pnt(-5,-5,0), Pnt(1,1,1))
    part = bigpart* top * bot
    outer = ((OrthoBrick(Pnt(-1,-1,-1), Pnt(2,2,3)).bc("outer") * Plane(Pnt(2,2,2), Vec(0,0,1)).bc("outertop")))
    
    geo.Add((part * outer).mat("part"))
    if version == 1:
        geo.Add((outer-part).mat("rest"))
    else:
        geo.Add((outer-top).mat("rest"))
        geo.Add((outer-bot).mat("rest"))
        geo.Add((outer*top*bot-bigpart).mat("rest"))

    geo.CloseSurfaces(top, bot, [])
    mesh = geo.GenerateMesh()
    should_ne = mesh.ne + 2 * GetNSurfaceElements(mesh, ["default"])
    layersize = 0.05
    mesh.BoundaryLayer("default", [0.5 * layersize, layersize], "layer", domains="part", outside=outside, grow_edges=True)
    assert mesh.ne == should_ne
    assert not "elements are not matching" in capfd.readouterr().out
