
import os, pytest
from netgen.meshing import meshsize, MeshingParameters, SetMessageImportance
import netgen.csg as csg
import netgen.stl as stl
import netgen.geom2d as geom2d
from pyngcore import TaskManager
import json
import hashlib
import numpy as np
try:
    import netgen.occ as occ
    has_occ = occ.occ_version >= "7.4.0"
    has_occ = False # skip occ tests (current development on interface)
except ImportError:
    has_occ = False

SetMessageImportance(0)

def round(x, digits=11):
    try:
        return float(("{:."+str(digits)+"g}").format(x))
    except: #list
        return [float(("{:."+str(digits)+"g}").format(y)) for y in x]


def saveNormalizedMesh(mesh, filename):
    mesh.Save(filename)
    with open(filename, 'rb') as f:
        lines = f.read().replace(b'\r\n', b'\n').split(b'\n')
    while lines and (lines[0].startswith(b'#') or lines[0].strip() == b''):
        lines.pop(0)
    data = b'\n'.join(lines)
    with open(filename, 'wb') as f:
        f.write(data)
    return hashlib.sha256(data).hexdigest()

def hashArray(a):
    return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()

def elementArray(elements, ncols=None):
    # copy fields, the structured array contains uninitialized padding bytes
    a = elements.NumPy()
    nodes = a['nodes']
    if ncols is not None and nodes.shape[1] < ncols:
        # elements store only as many node slots as needed, pad to the reference width
        from netgen.meshing import PointId
        nodes = np.pad(nodes, ((0, 0), (0, ncols - nodes.shape[1])), constant_values=PointId.base - 1)
    return np.column_stack([nodes, a['index']]).astype(np.int32)

def getData(mesh, mp, vol_filename):
    out = {}
    out['hash'] = saveNormalizedMesh(mesh, vol_filename)
    out['hash_el1d'] = hashArray(elementArray(mesh.Elements1D()))
    out['hash_el2d'] = hashArray(elementArray(mesh.Elements2D(), 8))
    out['hash_el3d'] = hashArray(elementArray(mesh.Elements3D(), 20))
    out['hash_points'] = hashArray(np.array(mesh.Coordinates(), dtype=np.float64))
    out['ne1d'] = len(mesh.Elements1D())
    out['ne2d'] = len(mesh.Elements2D())
    out['ne3d'] = len(mesh.Elements3D())
    # round badness to avoid fluctuations in last digits
    out["total_badness"] = round(mesh.CalcTotalBadness(mp))
    angles = mesh.CalcMinMaxAngle()
    out["angles_trig"] = round(angles["trig"], 5)
    out["angles_tet"] = round(angles["tet"], 5)
    out["quality_histogram"] = str(list(mesh.GetQualityHistogram()))
    return out

def checkData(mesh, mp, ref, vol_filename):
    data = getData(mesh, mp, vol_filename)
    assert ref['ne1d'] == data['ne1d']
    assert ref['ne2d'] == data['ne2d']
    assert ref['ne3d'] == data['ne3d']
    assert json.loads(ref['quality_histogram']) == json.loads(data['quality_histogram'])
    assert ref['total_badness'] == pytest.approx(data['total_badness'], rel=1e-10)
    assert ref['angles_trig'] == pytest.approx(data['angles_trig'], rel=1e-10)
    assert ref['angles_tet'] == pytest.approx(data['angles_tet'], rel=1e-10)
    hash_keys = ['hash_el1d', 'hash_el2d', 'hash_el3d', 'hash_points', 'hash']
    assert {k: ref[k] for k in hash_keys} == {k: data[k] for k in hash_keys}

# get tutorials
def getFiles(fileEnding):
    r, d, files = next(os.walk(os.path.join("..","..","tutorials")))
    return [f for f in files if f.endswith(fileEnding)]

# get additional tests
def getAdditionalFiles(fileEnding):
    r, d, files = next(os.walk("geofiles"))
    return [f for f in files if f.endswith(fileEnding)]

@pytest.fixture
def refdata():
    return json.load(open('results.json','r'))


def getMeshingparameters(filename):
    standard = [MeshingParameters()] + [MeshingParameters(ms)  for ms in (meshsize.very_coarse, meshsize.coarse, meshsize.moderate, meshsize.fine, meshsize.very_fine)]
    if filename == "plane.stl":
        return [] # do not test this example cause it needs so long...
    if filename == "shell.geo":
        return [] # do not test this example cause it needs so long...
    if filename == "manyholes2.geo":
        return [standard[1]] # this gets too big for finer meshsizes
    if filename in ("manyholes.geo", "frame.step"):
        return standard[:3] # this gets too big for finer meshsizes
    if filename == "extrusion.geo":
        return standard[:-1]
    if filename == "screw.step":
        return standard[3:] # coarser meshes don't work here
    if filename == "cylinder.geo":
        return [] # gives inconsistent reults
    if filename == "cylsphere.geo":
        return standard[0:2] + standard[3:] # coarse gives inconsistent reults (other mesh on MacOS)
    if filename == "part1.stl":
        return standard[0:1] + standard[2:] # very coarse does not work
    return standard

_geofiles =  getFiles(".in2d") + getFiles(".geo") + getFiles(".stl")
if has_occ:
    _geofiles += getFiles(".step")
_geofiles.sort()
_additional_testfiles = getAdditionalFiles(".stl")
if has_occ:
    _additional_testfiles += getAdditionalFiles(".step")
_additional_testfiles.sort()

def generateMesh(filename, mp):
    folder = os.path.join("..","..","tutorials") if filename in _geofiles else "geofiles"
    if filename.endswith(".geo"):
        geo = csg.CSGeometry(os.path.join(folder, filename))
    elif filename.endswith(".stl"):
        geo = stl.STLGeometry(os.path.join(folder, filename))
    elif filename.endswith(".step"):
        geo = occ.OCCGeometry(os.path.join(folder, filename))
    elif filename.endswith(".in2d"):
        geo = geom2d.SplineGeometry(os.path.join(folder, filename))
    return geo.GenerateMesh(mp)

def isSlowTest(filename):
    return filename in ["cubemcyl.geo", "frame.step", "revolution.geo", "manyholes.geo", "torus.geo",
                        "cubemsphere.geo", "manyholes2.geo", "matrix.geo", "trafo.geo", "ellipticcone.geo",
                        "period.geo", "shaft.geo", "cubeandring.geo", "ellipticcyl.geo",
                        "ellipsoid.geo", "cone.geo", "plane.stl"]

def getParameters():
    res = []
    for f in _geofiles + _additional_testfiles:
        parameters = getMeshingparameters(f)
        if os.environ.get("BUILD_TYPE") == "Debug":
            parameters = parameters[:1]
        for i,mp in enumerate(parameters):
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
    checkData(mesh, mp, ref[i], f'{filename}_{i}.vol')


def generateResultFile(output_file='results.json'):
    import time
    data = {}
    with TaskManager():
        for _file in _geofiles + _additional_testfiles:
            print("generate "+_file)
            start = time.time()
            mps = getMeshingparameters(_file)
            if not mps:
                continue
            meshdata = []
            for i, mp in enumerate(mps):
                try:
                    mesh = generateMesh(_file, mp)
                except Exception as e:
                    print("Meshingparameters: ", mp)
                    raise e
                meshdata.append( getData(mesh, mp, f'{_file}_{i}.vol') )
            data[_file] = meshdata
            print("needed", time.time() - start, "seconds")
        
    s = json.dumps(data, sort_keys=True, indent=4)
    open(output_file, "w").write(s)
    print("done")

if __name__ == "__main__":
    import sys
    generateResultFile(*sys.argv[1:])
