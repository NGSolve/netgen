"""
Mesh distribution over an MPI communicator: run with
    mpirun -np 3 python3 -m pytest test_mpi_distribute.py
Needs mpi4py only; with a single rank all tests are skipped.
"""
import pickle, warnings
import pytest

mpi4py = pytest.importorskip("mpi4py")
from mpi4py import MPI
import netgen.meshing as nm

comm = MPI.COMM_WORLD
pytestmark = pytest.mark.skipif(comm.size == 1, reason="needs more than one MPI rank")

try:
    from netgen.occ import Box, Pnt, OCCGeometry, IdentificationType, X, Y
    has_occ = True
except ImportError:
    has_occ = False


def unit_cube_mesh(maxh=0.3):
    from netgen.csg import unit_cube
    return unit_cube.GenerateMesh(maxh=maxh)


def distribute(generate, root_participates):
    """rank 0 generates, everybody gets its part; returns (mesh, serial copy or None)"""
    if comm.rank == 0:
        m = generate()
        ref = pickle.loads(pickle.dumps(m))
        return m.Distribute(comm, root_participates=root_participates), ref
    return nm.Mesh.Receive(comm), None


def signature(m):
    """point set and element vertex coordinates, independent of numbering"""
    pts = {i + 1: tuple(round(c, 9) for c in p.p) for i, p in enumerate(m.Points())}
    def els(elements):
        return sorted(tuple(sorted(pts[v.nr] for v in el.vertices)) for el in elements)
    return sorted(pts.values()), els(m.Elements3D()), els(m.Elements2D()), els(m.Elements1D())


@pytest.mark.parametrize("root_participates", [True, False])
def test_distribute_refine_curve(root_participates):
    mesh, ref = distribute(unit_cube_mesh, root_participates)
    ne_ref = comm.bcast(ref.ne if ref else None, root=0)

    assert mesh.comm.rank == comm.rank and mesh.comm.size == comm.size
    nes = comm.allgather(mesh.ne)
    assert sum(nes) == ne_ref
    assert (nes[0] > 0) == root_participates
    assert all(n > 0 for n in nes[1:])

    mesh.Refine()
    assert comm.allreduce(mesh.ne, op=MPI.SUM) == 8 * ne_ref
    mesh.Curve(2)


@pytest.mark.parametrize("root_participates", [True, False])
def test_gather_matches_serial_refinement(root_participates):
    mesh, ref = distribute(unit_cube_mesh, root_participates)
    mesh.Refine()
    gathered = mesh.Gather()
    if comm.rank == 0:
        ref.Refine()
        assert gathered.comm.size == 1
        assert signature(gathered) == signature(ref)
        assert gathered.GetRegionNames(dim=2) == ref.GetRegionNames(dim=2)
        # the gathered mesh is a plain serial mesh
        again = pickle.loads(pickle.dumps(gathered))
        assert again.ne == gathered.ne
    else:
        assert gathered is None


@pytest.mark.skipif(not has_occ, reason="OCC not available")
def test_periodic_occ_descriptors_names_identifications():
    def generate():
        box = Box(Pnt(0, 0, 0), Pnt(1, 1, 1))
        box.faces.Min(X).name = "left"
        box.faces.Max(X).name = "right"
        box.faces.Max(X).Identify(box.faces.Min(X), "periodic", IdentificationType.PERIODIC)
        box.edges.Min(X + Y).name = "myedge"
        box.solids.name = "mat"
        return OCCGeometry(box).GenerateMesh(maxh=0.3)

    mesh, ref = distribute(generate, True)
    info = comm.bcast((ref.ne, ref.GetNFaceDescriptors(), ref.GetNED(), ref.GetNrIdentifications(),
                       ref.GetRegionNames(dim=3), ref.GetRegionNames(dim=2), ref.GetRegionNames(dim=1))
                      if ref else None, root=0)
    assert comm.allreduce(mesh.ne, op=MPI.SUM) == info[0]
    assert mesh.GetNFaceDescriptors() == info[1]
    assert mesh.GetNED() == info[2]
    assert mesh.GetNrIdentifications() == info[3]
    assert mesh.GetRegionNames(dim=3) == info[4]
    assert mesh.GetRegionNames(dim=2) == info[5]
    assert mesh.GetRegionNames(dim=1) == info[6]
    for seg in mesh.Elements1D():
        assert seg.index >= 1
    mesh.Refine()
    mesh.Curve(3)
    assert comm.allreduce(mesh.ne, op=MPI.SUM) == 8 * info[0]


def test_pickling_is_local_and_collective_is_deprecated():
    mesh, ref = distribute(unit_cube_mesh, True)
    ne_ref = comm.bcast(ref.ne if ref else None, root=0)

    with warnings.catch_warnings():
        warnings.simplefilter("error")       # default: local, no warning
        local = pickle.loads(pickle.dumps(mesh))
    assert local.ne == mesh.ne

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        nm.SetParallelPickling(True)
        assert any(issubclass(x.category, DeprecationWarning) for x in w)
    try:
        collective = pickle.loads(pickle.dumps(mesh))   # collective: whole mesh on rank 0
        assert collective.ne == (ne_ref if comm.rank == 0 else mesh.ne)
    finally:
        nm.SetParallelPickling(False)
    comm.Barrier()
