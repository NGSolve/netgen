
import pytest
from netgen.csg import *

geometries=[unit_cube]

try:
    import netgen.occ as occ
    box = occ.Box( (0,0,0), (1,1,1) )
    box.faces.Min(occ.Y).name = "back"
    box.faces.Max(occ.Y).name = "front"
    box.faces.Min(occ.X).name = "left"
    box.faces.Max(occ.X).name = "right"
    box.faces.Min(occ.Z).name = "bottom"
    box.faces.Max(occ.Z).name = "top"
    geometries.append(occ.OCCGeometry(box))
except ImportError:
    pass

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
@pytest.mark.parametrize("geo", geometries)
def test_boundarylayer(outside, geo, capfd):
    mesh = geo.GenerateMesh(maxh=0.3)
    ne_before = mesh.ne
    layer_surfacenames = ["right", "top", "left", "back", "bottom"]
    mesh.BoundaryLayer("|".join(layer_surfacenames), [0.01, 0.01], "layer", outside=outside)

    should_ne = ne_before + 2 * GetNSurfaceElements(mesh, layer_surfacenames)
    assert mesh.ne == should_ne
    capture = capfd.readouterr()
    assert not "elements are not matching" in capture.out

    for side in ["front"]:
        mesh.BoundaryLayer(side, [0.001, 0.001], "layer", outside=outside)
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
    layersize = 0.025
    mesh.BoundaryLayer("default", [layersize, layersize], "part", domains="part", outside=outside)
    assert mesh.ne == should_ne
    assert not "elements are not matching" in capfd.readouterr().out
    # import netgen.gui
    ngs = pytest.importorskip("ngsolve")
    ngs.Draw(ngs.Mesh(mesh))
    mesh = ngs.Mesh(mesh)
    assert ngs.Integrate(1, mesh.Materials("part")) == pytest.approx(0.5*2.05*2.05 if outside else 0.4*2*2)
    assert ngs.Integrate(1, mesh) == pytest.approx(3**3)


@pytest.mark.parametrize("outside", [True, False])
def test_wrong_orientation(outside, capfd):
    geo = CSGeometry()
    brick = OrthoBrick((-1,0,0),(1,1,1)) - Plane((0,0,0), (1,0,0))
    geo.Add(brick.mat("air"))

    mesh = geo.GenerateMesh()

    mesh.BoundaryLayer(".*", 0.1, "air", domains="air", outside=outside)
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
    mesh.BoundaryLayer(".*", [0.001, 0.001], "block", "block", outside=False)
    ngs = pytest.importorskip("ngsolve")
    mesh = ngs.Mesh(mesh)
    assert ngs.Integrate(1, mesh) == pytest.approx(1)
    assert ngs.Integrate(1, mesh.Materials("slot")) == pytest.approx(0.4)

@pytest.mark.parametrize("outside", [True, False])
def test_pyramids(outside):
    geo = CSGeometry()
    box = OrthoBrick((0,0,0), (1,1,1))
    plate = OrthoBrick((0.3,0.3,0.4),(0.7,0.7,1)) * Plane((0,0,0.6), (0,0,1)).bc("top")
    geo.Add((box-plate).mat("air"))
    geo.Add(plate.mat("plate"))
    mesh = geo.GenerateMesh()
    mesh.BoundaryLayer("top", [0.01], "layer", "plate", outside=outside)
    ngs = pytest.importorskip("ngsolve")
    mesh = ngs.Mesh(mesh)
    assert ngs.Integrate(1, mesh.Materials("plate")) == pytest.approx(0.032 if outside else 0.0304)
    assert ngs.Integrate(1, mesh.Materials("layer")) == pytest.approx(0.0016)
    assert ngs.Integrate(1, mesh.Materials("air")) == pytest.approx(0.9664 if outside else 0.968)

# not working yet
@pytest.mark.parametrize("outside", [True, False])
def _test_with_inner_corner(outside, capfd):
    geo = CSGeometry()

    core_thickness = 0.1
    limb_distance = 0.5
    core_height = 0.5
    coil_r1 = 0.08
    coil_r2 = 0.16
    coil_h = 0.2
    domain_size = 1.2
    domain_size_y = 0.7

    def CreateCoil(x):
        outer = Cylinder((x, 0, -1), (x, 0, 1), coil_r2)
        inner = Cylinder((x, 0, -1), (x, 0, 1), coil_r1)
        top = Plane((0,0,coil_h/2), (0,0,1))
        bot = Plane((0,0,-coil_h/2), (0,0,-1))
        return ((outer - inner) * top * bot, (outer, inner, top, bot))

    core_front = Plane((0,-core_thickness/2, 0), (0,-1,0)).bc("core_front")
    core_back = Plane((0,core_thickness/2, 0), (0,1,0)).bc("core_front")
    core_limb1 = OrthoBrick((-limb_distance/2-core_thickness/2, -core_thickness, -core_height/2),(-limb_distance/2+core_thickness/2, core_thickness, core_height/2))
    core_limb2 = OrthoBrick((limb_distance/2-core_thickness/2, -core_thickness, -core_height/2),(limb_distance/2+core_thickness/2, core_thickness, core_height/2))
    core_top = OrthoBrick((-limb_distance/2-core_thickness/2, -core_thickness, core_height/2-core_thickness/2),(limb_distance/2+core_thickness/2, core_thickness, core_height/2+core_thickness/2))
    core_bot = OrthoBrick((-limb_distance/2-core_thickness/2, -core_thickness, -core_height/2-core_thickness/2),(limb_distance/2+core_thickness/2, core_thickness, -core_height/2+core_thickness/2))

    core = (core_limb1 + core_limb2 + core_top + core_bot).bc("core_rest")
    core = core * core_front * core_back
    core.maxh(core_thickness * 0.4)

    coil1, (outer1, inner1, top1, bot1) = CreateCoil(-limb_distance/2)
    coil1.mat("coil_1")
    coil2, (outer2, inner2, top2, bot2) = CreateCoil(limb_distance/2)
    coil2.mat("coil_2")

    oil = OrthoBrick((-domain_size/2, -domain_size_y/2, -domain_size/2), (domain_size/2, domain_size_y/2, domain_size/2)).bc("tank") - core # - coil1 - coil2

    geo.Add(core.mat("core"), col=(0.4,0.4,0.4))
    geo.Add(coil1, col=(0.72, 0.45, 0.2))
    geo.Add(coil2, col=(0.72, 0.45, 0.2))
    geo.Add(oil.mat("oil"), transparent=True)
    mesh = geo.GenerateMesh()
    mesh.BoundaryLayer("core_front", [0.001, 0.002], "core", "core", outside=outside)
    ngs = pytest.importorskip("ngsolve")
    mesh = ngs.Mesh(mesh)
    capture = capfd.readouterr()
    assert not "elements are not matching" in capture.out
    assert ngs.Integrate(1, mesh.Materials("core")) == pytest.approx(0.0212 if outside else 0.02)
    assert ngs.Integrate(1, mesh.Materials("oil")) == pytest.approx(0.9868 if outside else 0.988)
