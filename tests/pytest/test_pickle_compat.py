"""
Mesh pickles written by older netgen releases must stay readable.

pickle_compat/<release>/<case>.pkl.gz are netgen meshes pickled by the PyPI
wheel of that release, pickle_compat/expected.json holds what the mesh
contains. New releases are added with pickle_compat/add_release.py.
"""

import gzip, hashlib, json, pickle
from pathlib import Path

import pytest

datadir = Path(__file__).parent / "pickle_compat"
expected = json.loads((datadir / "expected.json").read_text())
cases = [(release, case) for release in sorted(expected) for case in sorted(expected[release])]


def fingerprint(mesh):
    """content of a netgen mesh, independent of the numbering base of points"""
    dim = mesh.dim
    points = list(mesh.Points())
    elements = {3: mesh.Elements3D(), 2: mesh.Elements2D(), 1: mesh.Elements1D(), 0: mesh.Elements0D()}

    def elhash(els):
        h = hashlib.md5()
        for el in els:
            h.update(("%d:" % el.index).encode())
            for v in el.vertices:
                p = mesh[v]
                h.update(("%.9g,%.9g,%.9g;" % (p[0], p[1], p[2])).encode())
        return h.hexdigest()

    fp = {"dim": dim, "np": len(points),
          "coordsum": [round(sum(p[i] for p in points), 9) for i in range(3)]}
    for codim in range(dim + 1):
        els = elements[dim - codim]
        fp["nel_codim%d" % codim] = len(els)
        fp["elhash_codim%d" % codim] = elhash(els)
        fp["names_codim%d" % codim] = list(mesh.GetRegionNames(codim=codim))
    fp["identifications"] = len(mesh.GetIdentifications())
    fp["curveorder"] = mesh.GetCurveOrder()
    fp["has_geometry"] = mesh.GetGeometry() is not None
    return fp


def names(lst):
    """region names: unnamed is 'default', and trailing unnamed regions do not count"""
    lst = [name or "default" for name in lst]
    while lst and lst[-1] == "default":
        lst.pop()
    return lst


def same(a, b):
    if isinstance(a, float) or isinstance(b, float):
        return a == pytest.approx(b, abs=1e-8)
    if isinstance(a, list) and isinstance(b, list):
        if all(isinstance(x, str) for x in a + b):
            return names(a) == names(b)
        return len(a) == len(b) and all(same(x, y) for x, y in zip(a, b))
    return a == b


def test_golden_files_present():
    assert len(cases) >= 84
    for release, case in cases:
        assert (datadir / release / (case + ".pkl.gz")).exists(), (release, case)


@pytest.mark.parametrize("release,case", cases)
def test_read_old_pickle(release, case):
    exp = expected[release][case]
    if exp.get("needs_occ"):
        pytest.importorskip("netgen.occ")
    with gzip.open(datadir / release / (case + ".pkl.gz"), "rb") as f:
        mesh = pickle.load(f)
    fp = fingerprint(mesh)
    known_bugs = exp.get("known_bugs", {})
    wrong = [key for key, val in exp["fingerprint"].items() if key not in known_bugs and not same(fp[key], val)]
    assert not wrong, {key: (fp[key], exp["fingerprint"][key]) for key in wrong}
    still_wrong = [key for key in known_bugs if not same(fp[key], exp["fingerprint"][key])]
    if still_wrong:
        pytest.xfail("; ".join(key + ": " + known_bugs[key] for key in still_wrong))


def test_truncated_archive_raises():
    """an unreadable archive must end in an exception, not in allocation without bound"""
    from netgen.csg import unit_cube
    from netgen.meshing import Mesh
    state = unit_cube.GenerateMesh(maxh=0.5).__getstate__()[0]
    big = max((x for x in state if isinstance(x, bytes)), key=len)
    for cut in (len(big) // 2, len(big) // 3, len(big) - 100):
        broken = [x[:cut] if x is big else x for x in state]
        mesh = Mesh.__new__(Mesh)
        with pytest.raises(Exception):
            mesh.__setstate__((broken,))
