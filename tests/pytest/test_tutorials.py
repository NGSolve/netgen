
import os, pytest
from netgen.meshing import meshsize, MeshingParameters, SetMessageImportance
import netgen.csg as csg
import netgen.stl as stl
from pyngcore import TaskManager
import json
try:
    import netgen.occ as occ
    has_occ = True
except ImportError:
    has_occ = False

SetMessageImportance(0)

def getData(mesh, mp):
    out = {}
    out['ne1d'] = len(mesh.Elements1D())
    out['ne2d'] = len(mesh.Elements2D())
    out['ne3d'] = len(mesh.Elements3D())
    # round badness to avoid fluctuations in last digits
    out["total_badness"] = float("{:.11g}".format(mesh.CalcTotalBadness(mp)))
    out["quality_histogram"] = str(list(mesh.GetQualityHistogram()))
    return out

def checkData(mesh, mp, ref):
    data = getData(mesh, mp)
    assert ref['ne1d'] == data['ne1d']
    assert ref['ne2d'] == data['ne2d']
    assert ref['ne3d'] == data['ne3d']
    assert ref['quality_histogram'] == data['quality_histogram']
    assert ref['total_badness'] == pytest.approx(data['total_badness'], rel=1e-5)

# get tutorials
def getFiles(fileEnding):
    r, d, files = next(os.walk(os.path.join("..","..","tutorials")))
    return (f for f in files if f.endswith(fileEnding))

# get additional tests
def getAdditionalFiles(fileEnding):
    r, d, files = next(os.walk("geofiles"))
    return (f for f in files if f.endswith(fileEnding))

@pytest.fixture
def refdata():
    return json.load(open('results.json','r'))


def getMeshingparameters(filename):
    standard = [MeshingParameters()] + [MeshingParameters(ms)  for ms in (meshsize.very_coarse, meshsize.coarse, meshsize.moderate, meshsize.fine, meshsize.very_fine)]
    if filename == "shell.geo":
        return [] # do not test this example cause it needs so long...
    if filename == "extrusion.geo":
        return [] # this segfaults right now
    if filename == "manyholes2.geo":
        return [standard[1]] # this gets too big for finer meshsizes
    if filename in ("manyholes.geo", "frame.step"):
        return standard[:3] # this gets too big for finer meshsizes
    if filename == "screw.step":
        return standard[3:] # coarser meshes don't work here
    return standard

_geofiles = [f for f in getFiles(".geo")] + [f for f in getFiles(".stl")]
if has_occ:
    _geofiles += [f for f in getFiles(".step")]
_geofiles.sort()
_additional_testfiles = [f for f in getAdditionalFiles(".stl")]
if has_occ:
    _additional_testfiles += [f for f in getAdditionalFiles(".step")]
_additional_testfiles.sort()

def generateMesh(filename, mp):
    folder = os.path.join("..","..","tutorials") if filename in _geofiles else "geofiles"
    if filename.endswith(".geo"):
        geo = csg.CSGeometry(os.path.join(folder, filename))
    elif filename.endswith(".stl"):
        geo = stl.STLGeometry(os.path.join(folder, filename))
    elif filename.endswith(".step"):
        geo = occ.OCCGeometry(os.path.join(folder, filename))
    return geo.GenerateMesh(mp)

def isSlowTest(filename):
    return filename in ["cubemcyl.geo", "frame.step", "revolution.geo", "manyholes.geo", "torus.geo",
                        "cubemsphere.geo", "manyholes2.geo", "matrix.geo", "trafo.geo", "ellipticcone.geo",
                        "period.geo", "shaft.geo", "cubeandring.geo", "ellipticcyl.geo",
                        "ellipsoid.geo", "cone.geo", "plane.stl"]

def getParameters():
    res = []
    for f in _geofiles + _additional_testfiles:
        for i,mp in enumerate(getMeshingparameters(f)):
            if isSlowTest(f):
                res.append( pytest.param(f, mp, i, marks=pytest.mark.slow ) )
            else:
                res.append( (f, mp, i) )
    return res

@pytest.mark.parametrize(("filename", "mp", "i"), getParameters())
def test_geoFiles(filename, mp, i, refdata):
    ref = refdata[filename]
    import filecmp
    print("load geo", filename)
    mp = MeshingParameters(mp, parallel_meshing=False)
    mesh = generateMesh(filename, mp)
    mesh.Save(filename+'_seq.vol.gz')
    with TaskManager():
        mesh_par = generateMesh(filename, mp)
        mesh_par.Save(filename+'_par.vol.gz')

    assert filecmp.cmp(filename+'_seq.vol.gz', filename+'_par.vol.gz')
    checkData(mesh, mp, ref[i])


def generateResultFile():
    import re, time
    data = {}
    with TaskManager():
        for _file in _geofiles + _additional_testfiles:
            print("generate "+_file)
            start = time.time()
            mps = getMeshingparameters(_file)
            if not mps:
                continue
            meshdata = []
            for mp in mps:
                mesh = generateMesh(_file, mp)
                meshdata.append( getData(mesh, mp) )
            data[_file] = meshdata
            print("needed", time.time() - start, "seconds")
        
    s = json.dumps(data, sort_keys=True, indent=4)
    open("results.json", "w").write(s)
    print("done")

if __name__ == "__main__":
    generateResultFile()
