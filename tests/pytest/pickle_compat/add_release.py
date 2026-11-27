"""
Add the golden pickles of one netgen release.

  1. write the pickles with the release itself (any python that has the wheel
     and ngsolve of the same version):
         python add_release.py write <outdir>
  2. take them over with the current build, which checks its reading against
     what the writer saw:
         python add_release.py add <release> <outdir>

A field that the writer already lost in its own pickle round trip is stored
as the file has it. Any other difference is listed under known_bugs with the
writer's value as the expected one, and wants a look.
"""

import gzip, json, pickle, sys
from pathlib import Path

here = Path(__file__).parent


def meshes():
    from netgen.geom2d import unit_square, SplineGeometry
    from netgen.csg import CSGeometry, OrthoBrick, Sphere, Plane, Pnt, Vec, unit_cube
    from ngsolve import Mesh
    from ngsolve.meshes import MakeStructured3DMesh, MakeStructured2DMesh, Make1DMesh

    def circgeo():
        geo = SplineGeometry()
        geo.AddRectangle((-1, -1), (1, 1), bcs=["b", "r", "t", "l"], leftdomain=1, rightdomain=0)
        geo.AddCircle((0, 0), 0.4, bc="circ", leftdomain=2, rightdomain=1)
        geo.SetMaterial(1, "outer")
        geo.SetMaterial(2, "inner")
        return geo

    def curved(ngmesh, order=3):
        m = Mesh(ngmesh); m.Curve(order); return m

    def csg3d():
        cube = OrthoBrick(Pnt(-1, -1, -1), Pnt(1, 1, 1)).bc("outer")
        sph = Sphere(Pnt(0, 0, 0), 0.5).bc("sph")
        geo = CSGeometry()
        geo.Add((cube - sph).mat("shell"))
        geo.Add(sph.mat("ball"))
        return curved(geo.GenerateMesh(maxh=0.6))

    def csgper():
        left = Plane(Pnt(0, 0, 0), Vec(-1, 0, 0)).bc("left")
        right = Plane(Pnt(1, 0, 0), Vec(1, 0, 0)).bc("right")
        front = Plane(Pnt(0, 0, 0), Vec(0, -1, 0)).bc("front")
        back = Plane(Pnt(0, 1, 0), Vec(0, 1, 0)).bc("back")
        bot = Plane(Pnt(0, 0, 0), Vec(0, 0, -1)).bc("bot")
        top = Plane(Pnt(0, 0, 1), Vec(0, 0, 1)).bc("top")
        geo = CSGeometry()
        geo.Add((left * right * front * back * bot * top).mat("cube"))
        geo.PeriodicSurfaces(left, right)
        return Mesh(geo.GenerateMesh(maxh=0.4))

    def occ3d():
        from netgen.occ import Box, Cylinder, Pnt, Z, OCCGeometry, Glue
        box = Box(Pnt(0, 0, 0), Pnt(2, 1, 1))
        box.faces.name = "wall"
        box.faces.Min(Z).name = "floor"
        box.edges.name = "boxedge"
        box.vertices.name = "corner"
        cyl = Cylinder(Pnt(1, 0.5, 0), Z, r=0.3, h=1)
        cyl.faces.name = "cyl"
        cyl.mat("rod")
        rest = box - cyl
        rest.mat("block")
        return curved(OCCGeometry(Glue([rest, cyl])).GenerateMesh(maxh=0.5))

    def occ2d():
        from netgen.occ import WorkPlane, OCCGeometry, Glue
        rect = WorkPlane().Rectangle(2, 1).Face()
        rect.edges.name = "rim"
        circ = WorkPlane().Circle(1, 0.5, 0.25).Face()
        circ.edges.name = "circ"
        circ.faces.name = "disk"
        rest = rect - circ
        rest.faces.name = "plate"
        return curved(OCCGeometry(Glue([rest, circ]), dim=2).GenerateMesh(maxh=0.3))

    def occsurf():
        from netgen.occ import Sphere, Pnt, OCCGeometry
        sp = Sphere(Pnt(0, 0, 0), 1)
        sp.faces.name = "sphere"
        return curved(OCCGeometry(sp.faces[0]).GenerateMesh(maxh=0.5), 2)

    def refined():
        m = Mesh(unit_cube.GenerateMesh(maxh=0.5)); m.Refine(); return m

    def refined2d():
        m = curved(circgeo().GenerateMesh(maxh=0.4)); m.Refine(); return m

    return {"sq2d": lambda: Mesh(unit_square.GenerateMesh(maxh=0.3)),
            "circ2d": lambda: curved(circgeo().GenerateMesh(maxh=0.4)),
            "quad2d": lambda: curved(circgeo().GenerateMesh(maxh=0.4, quad_dominated=True)),
            "csg3d": csg3d, "csgper": csgper, "occ3d": occ3d, "occ2d": occ2d, "occsurf": occsurf,
            "hex3d": lambda: MakeStructured3DMesh(hexes=True, nx=2, ny=2, nz=2),
            "prism3d": lambda: MakeStructured3DMesh(hexes=False, prism=True, nx=2, ny=2, nz=2),
            "quadstruct": lambda: MakeStructured2DMesh(quads=True, nx=3, ny=3),
            "seg1d": lambda: Make1DMesh(n=5),
            "refined": refined, "refined2d": refined2d}


def writer_view(mesh):
    """what the writing release sees, through ngsolve (stable across the releases)"""
    from ngsolve import VOL, BND, BBND
    d = {"dim": mesh.dim, "np": mesh.nv,
         "nel_codim0": mesh.ne, "nel_codim1": mesh.GetNE(BND), "nel_codim2": mesh.GetNE(BBND),
         "names_codim0": list(mesh.GetMaterials()), "names_codim1": list(mesh.GetBoundaries()),
         "names_codim2": list(mesh.GetBBoundaries()),
         "coordsum": [round(sum(v.point[i] for v in mesh.vertices), 9) for i in range(mesh.dim)]}
    return d


def write(outdir):
    from ngsolve import SetNumThreads
    SetNumThreads(1)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    info = {}
    for name, make in meshes().items():
        mesh = make()
        view = writer_view(mesh)
        raw = pickle.dumps(mesh.ngmesh)
        (outdir / (name + ".pkl")).write_bytes(raw)
        from ngsolve import Mesh
        view2 = writer_view(Mesh(pickle.loads(raw)))
        info[name] = {"view": view, "lost_in_roundtrip": sorted(k for k in view if view[k] != view2[k])}
    (outdir / "writer.json").write_text(json.dumps(info, indent=1))


def add(release, outdir):
    sys.path.insert(0, str(here.parent))
    from test_pickle_compat import fingerprint, same
    outdir = Path(outdir)
    info = json.loads((outdir / "writer.json").read_text())
    expfile = here / "expected.json"
    expected = json.loads(expfile.read_text()) if expfile.exists() else {}
    expected[release] = {}
    (here / release).mkdir(exist_ok=True)
    for name, rec in info.items():
        raw = (outdir / (name + ".pkl")).read_bytes()
        fp = fingerprint(pickle.loads(raw))
        entry = {"fingerprint": fp}
        if name.startswith("occ"):
            entry["needs_occ"] = True
        bugs = {}
        for key, val in rec["view"].items():
            if key not in fp:
                continue
            if key == "coordsum":
                val = val + [0.0] * (3 - len(val))
            if same(fp[key], val) or key in rec["lost_in_roundtrip"]:
                continue
            fp[key] = val
            bugs[key] = "reader differs from the writing release"
        if bugs:
            entry["known_bugs"] = bugs
            print(release, name, "known_bugs:", sorted(bugs))
        expected[release][name] = entry
        with gzip.GzipFile(here / release / (name + ".pkl.gz"), "wb", mtime=0) as f:
            f.write(raw)
    expfile.write_text(json.dumps(expected, indent=1, sort_keys=True) + "\n")


if __name__ == "__main__":
    if sys.argv[1] == "write":
        write(sys.argv[2])
    else:
        add(sys.argv[2], sys.argv[3])
