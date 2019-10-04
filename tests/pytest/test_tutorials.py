
import os, pytest
from netgen.meshing import meshsize, MeshingParameters, SetMessageImportance
import netgen.csg as csg
import netgen.stl as stl
from pyngcore import TaskManager
try:
    import netgen.occ as occ
    has_occ = True
except ImportError:
    has_occ = False
from results import *

SetMessageImportance(0)

def getFiles(fileEnding):
    r, d, files = next(os.walk(os.path.join("..","..","tutorials")))
    return (f for f in files if f.endswith(fileEnding))


def getCheckFunc(filename):
    def func(mesh,mp,i):
        if filename in number_elements:
            # number of elements should be in 2% range of expected value
            assert mesh.ne == number_elements[filename][i]
            badness = mesh.CalcTotalBadness(mp)
            qual_classes = list(mesh.GetQualityHistogram())
            assert badness == pytest.approx(total_badness[filename][i], rel=1e-6)
            assert qual_classes == quality_histogram[filename][i]
    return func
            
def getResultFunc(filename):
    def resultFunc(mesh, mp):
        results = {}
        results["number_elements"] = mesh.ne
        results["total_badness"] = mesh.CalcTotalBadness(mp)
        results["quality_histogram"] = list(mesh.GetQualityHistogram())
        return results
    return resultFunc

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

def generateMesh(filename, mp):
    if filename.endswith(".geo"):
        geo = csg.CSGeometry(os.path.join("..","..","tutorials", filename))
    elif filename.endswith(".stl"):
        geo = stl.STLGeometry(os.path.join("..","..","tutorials", filename))
    elif filename.endswith(".step"):
        geo = occ.OCCGeometry(os.path.join("..","..","tutorials", filename))
    return geo.GenerateMesh(mp)

def isSlowTest(filename):
    return filename in ["cubemcyl.geo", "frame.step", "revolution.geo", "manyholes.geo", "torus.geo",
                        "cubemsphere.geo", "manyholes2.geo", "matrix.geo", "trafo.geo", "ellipticcone.geo",
                        "period.geo", "shaft.geo", "cubeandring.geo", "ellipticcyl.geo",
                        "ellipsoid.geo", "cone.geo"]

def getParamForTest(filename):
    return pytest.param(filename, getCheckFunc(filename), marks=pytest.mark.slow) if isSlowTest(filename) \
        else (filename, getCheckFunc(filename))

@pytest.mark.parametrize(("filename, checkFunc"), [getParamForTest(f) for f in _geofiles])
def test_geoFiles(filename, checkFunc):
    import filecmp
    for i, mp_ in enumerate(getMeshingparameters(filename)):
        print("load geo", filename)
        mp = MeshingParameters(mp_, parallel_meshing=False)
        mesh = generateMesh(filename, mp)
        if checkFunc is not None:
            checkFunc(mesh,mp,i)
        mesh.Save(filename+'_seq.vol.gz')

        with TaskManager():
            mesh_par = generateMesh(filename, mp)
            mesh_par.Save(filename+'_par.vol.gz')

        assert filecmp.cmp(filename+'_seq.vol.gz', filename+'_par.vol.gz')

import time
def generateResultFile():
  with TaskManager():
    with open("results.py", "w") as f:
        print("number_elements = {}", file=f)
        print("total_badness = {}", file=f)
        print("quality_histogram = {}", file=f)
        for _file, _func in ((gf, getResultFunc(gf)) for gf in _geofiles):
            start = time.time()
            print("write", _file)
            mps = getMeshingparameters(_file)
            if not mps:
                continue
            results = [_func(generateMesh(_file, mp), mp) for mp in mps]
            print("number_elements['{}'] = {}".format(_file, "(" + ",".join((str(r["number_elements"]) for r in results)) + ")"), file=f)
            print("total_badness['{}'] = {}".format(_file, "(" + ",".join((str(r["total_badness"]) for r in results)) + ")"), file=f)
            print("quality_histogram['{}'] = {}".format(_file, "(" + ",".join((str(r["quality_histogram"]) for r in results)) + ")"), file=f)
            print("needed", time.time() - start, "seconds")
        

if __name__ == "__main__":
    generateResultFile()
