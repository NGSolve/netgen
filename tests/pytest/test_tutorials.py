
import os, pytest
from netgen.meshing import meshsize, MeshingParameters, SetMessageImportance
import netgen.gui
import netgen.csg as csg
import netgen.stl as stl
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
    def func(mesh,i):
        if filename in number_elements:
            # number of elements should be in 2% range of expected value
            assert mesh.ne == pytest.approx(number_elements[filename][i], rel=0.02)
        return func
            
def getResultFunc(filename):
    def resultFunc(mesh):
        results = {}
        results["number_elements"] = mesh.ne
        return results
    return resultFunc

def getMeshingparameters(filename):
    standard = [{}] + [{ "mp" : ms } for ms in (meshsize.very_coarse, meshsize.coarse, meshsize.moderate, meshsize.fine, meshsize.very_fine)]
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

# TODO: step files do not respect gui meshsizes yet.
_geofiles = [f for f in getFiles(".geo")] + [f for f in getFiles(".stl")]
if has_occ:
    _geofiles += [f for f in getFiles(".step")]


def generateMesh(filename, mp):
    if filename.endswith(".geo"):
        geo = csg.CSGeometry(os.path.join("..","..","tutorials", filename))
    elif filename.endswith(".stl"):
        geo = stl.STLGeometry(os.path.join("..","..","tutorials", filename))
    elif filename.endswith(".step"):
        geo = occ.OCCGeometry(os.path.join("..","..","tutorials", filename))
    return geo.GenerateMesh(**mp)

@pytest.mark.parametrize("filename, checkFunc", [(f, getCheckFunc(f)) for f in _geofiles])
def test_geoFiles(filename, checkFunc):
    for i, mp in enumerate(getMeshingparameters(filename)):
        print("load geo", filename)
        mesh = generateMesh(filename, mp)
        if checkFunc is not None:
            checkFunc(mesh,i)

import time
def generateResultFile():
    with open("results.py", "w") as f:
        print("number_elements = {}", file=f)
        for _file, _func in ((gf, getResultFunc(gf)) for gf in _geofiles):
            start = time.time()
            print("write", _file)
            mps = getMeshingparameters(_file)
            if not mps:
                continue
            results = [_func(generateMesh(_file, mp)) for mp in mps]
            print("number_elements['{}'] = {}".format(_file, "(" + ",".join((str(r["number_elements"]) for r in results)) + ")"), file=f)
            print("needed", time.time() - start, "seconds")
        

if __name__ == "__main__":
    generateResultFile()
