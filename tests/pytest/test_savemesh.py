
from netgen.csg import *
from netgen import meshing
import filecmp
import difflib
from math import sqrt, cos, sin

def CreateQuad():
    base = Plane(Pnt(0,0,0),Vec(0,0,1))
    surface = SplineSurface(base)
    pts = [(-0.2,-0.2,0),(-0.2,0.2,0),(0.2,0.2,0),(0.2,-0.2,0)]
    geopts = [surface.AddPoint(*p) for p in pts]
    for p1,p2,bc in [(0,1,"wire"), (1, 2,"contact"),(2,3,"wire"),(3,0,"wire")]:
        surface.AddSegment(geopts[p1],geopts[p2],bc)
    return surface

Cross = lambda a,b: [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-b[0]*a[1]]


def CreateGeo():
    geo = CSGeometry()
    air = OrthoBrick(Pnt(-1,-1,-1),Pnt(1,1,1))
    geo.Add(air.mat("air"))
    surface = CreateQuad()
    geo.AddSplineSurface(surface)
    return geo

def test_BBNDsave():
    mesh = CreateGeo().GenerateMesh(maxh=0.4,perfstepsend = meshing.MeshingStep.MESHSURFACE)
    for i in range(2):
        mesh.GenerateVolumeMesh(mp = MeshingParameters(only3D_domain=i+1,maxh=0.4))
    mesh.SetGeometry(None)
    mesh.Save("test.vol")
    mesh2 = meshing.Mesh()
    mesh2.Load("test.vol")
    mesh2.Save("test2.vol")
    with open("test.vol","r") as f:
        first = f.readlines()
    with open("test2.vol","r") as f:
        second = f.readlines()
    # exclude the face colours section (because they aren't in the same order)
    for i,line in enumerate(first):
        if line[0:12] == "face_colours":
            first = first[0:i]
            second = second[0:i]
            break
    diff = difflib.context_diff(first,second)
    print("File diff:")
    l = list(diff)
    print(*l)
    assert len(l)==0

