
import pytest
from netgen.csg import *

def GetNSurfaceElements(mesh, boundaries, adjacent_domain=None):
    nse_in_layer = 0
    for el in mesh.Elements2D():
        if mesh.GetBCName(el.index-1) in boundaries:
            if adjacent_domain is None:
                nse_in_layer += 1
            else:
                if (mesh.FaceDescriptor(el.index).domin > 0 and mesh.GetMaterial(mesh.FaceDescriptor(el.index).domin) == adjacent_domain) or (mesh.FaceDescriptor(el.index).domout > 0 and mesh.GetMaterial(mesh.FaceDescriptor(el.index).domout) == adjacent_domain):
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
@pytest.mark.parametrize("version", [1, 2]) # version 2 not working yet
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
    should_ne = mesh.ne + 2 * GetNSurfaceElements(mesh, ["default"], "part")
    layersize = 0.05
    mesh.BoundaryLayer("default", [0.5 * layersize, layersize], "part", domains="part", outside=outside, grow_edges=True)
    assert mesh.ne == should_ne
    assert not "elements are not matching" in capfd.readouterr().out
    import netgen.gui
    ngs = pytest.importorskip("ngsolve")
    ngs.Draw(ngs.Mesh(mesh))
    mesh = ngs.Mesh(mesh)
    assert ngs.Integrate(1, mesh.Materials("part")) == pytest.approx(0.5*2.05*2.05 if outside else 0.4*2*2)
    assert ngs.Integrate(1, mesh) == pytest.approx(3**3)


@pytest.mark.parametrize("outside", [True, False])
def test_wrong_orientation(outside):
    geo = CSGeometry()
    brick = OrthoBrick((-1,0,0),(1,1,1)) - Plane((0,0,0), (1,0,0))
    geo.Add(brick.mat("air"))

    mesh = geo.GenerateMesh()

    mesh.BoundaryLayer(".*", 0.1, "air", domains="air", outside=outside,
                       grow_edges=True)
    ngs = pytest.importorskip("ngsolve")
    mesh = ngs.Mesh(mesh)
    assert ngs.Integrate(1, mesh) == pytest.approx(1.2**3 if outside else 1)

def test_splitted_surface():
    geo = CSGeometry()

    brick = OrthoBrick((0,0,0), (1,1,1))
    slots = OrthoBrick((0.2,0,-1), (0.4, 1, 2)) + OrthoBrick((0.6, 0,-1), (0.8, 1,2))
    geo.Add((brick-slots).mat("block"))
    geo.Add((brick*slots).mat("slot"))

    mesh = geo.GenerateMesh()
    mesh.BoundaryLayer(".*", [0.001, 0.002], "block", "block", outside=False,
                       grow_edges=True)
    ngs = pytest.importorskip("ngsolve")
    mesh = ngs.Mesh(mesh)
    assert ngs.Integrate(1, mesh) == pytest.approx(1)
    assert ngs.Integrate(1, mesh.Materials("slot")) == pytest.approx(0.4)
