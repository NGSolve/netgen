import pytest

from netgen.meshing import IdentificationType
idtype = IdentificationType.CLOSESURFACES

def test_two_boxes():
    occ = pytest.importorskip("netgen.occ")
    inner = occ.Box((0,0,0), (1,1,1))
    trafo = occ.gp_Trsf().Scale(inner.center, 1.1)
    outer = trafo(inner)

    inner.Identify(outer, "", idtype, trafo)
    shape = occ.Glue([outer-inner, inner])

    geo = occ.OCCGeometry(shape)
    mesh = geo.GenerateMesh(maxh=0.3)
    have_prisms = False

    for el in mesh.Elements3D():
        if len(el.vertices)==6:
            have_prisms = True
            break

    assert have_prisms

def test_two_circles():
    occ = pytest.importorskip("netgen.occ")
    circ1 = occ.WorkPlane().Circle(1).Face()
    trafo = occ.gp_Trsf().Scale(circ1.center, 1.1)

    circ2 = trafo(circ1)
    circ1.edges[0].Identify(circ2.edges[0], "", idtype, trafo)
    circ2 -= circ1
    shape = occ.Glue([circ1, circ2])

    geo = occ.OCCGeometry(shape, 2)
    mesh = geo.GenerateMesh(maxh=0.2)
    have_quads = False

    for el in mesh.Elements2D():
        if len(el.vertices)==4:
            have_quads = True
            break

    assert have_quads

def test_cut_identified_face():
    occ = pytest.importorskip("netgen.occ")
    from netgen.occ import Z, Box, Cylinder, Glue, OCCGeometry
    box = Box((-1,-1,0), (1,1,1))
    cyl = Cylinder( (0,0,0), Z, 0.5, 1 )

    box.faces.Min(Z).Identify(box.faces.Max(Z), "", idtype)
    shape = Glue([cyl, box])
    geo = OCCGeometry(shape)
    mesh = geo.GenerateMesh(maxh=0.5)

    for el in mesh.Elements3D():
        assert len(el.vertices)==6

def test_identify_multiple_faces():
    occ = pytest.importorskip("netgen.occ")
    from netgen.occ import Z, Box, Cylinder, Glue, OCCGeometry, gp_Trsf
    box = Box((-1,-1,0), (1,1,1))
    cyl = Cylinder( (0,0,0), Z, 0.5, 1 )

    shape = Glue([box, cyl])
    bot_faces = shape.faces[Z < 0.1]
    top_faces = shape.faces[Z > 0.1]
    bot_faces.Identify(top_faces, "", idtype, gp_Trsf.Translation((0,0,1)))

    geo = OCCGeometry(shape)
    mesh = geo.GenerateMesh(maxh=0.3)

    for el in mesh.Elements3D():
        assert len(el.vertices)==6
