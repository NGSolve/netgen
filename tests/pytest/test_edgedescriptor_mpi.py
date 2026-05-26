import pytest
import netgen.meshing

mpi4py = pytest.importorskip("mpi4py")
_ = pytest.importorskip("pytest_mpi")

try:
    from netgen.occ import Box, Pnt, OCCGeometry
    has_occ = True
except ImportError:
    has_occ = False


@pytest.mark.mpi
@pytest.mark.skipif(not has_occ, reason="OCC not available")
def test_mpi_segments_have_valid_ed_after_load():
    """Check that every local segment has a valid edge descriptor index (>= 0) after MPI distribution."""
    comm = mpi4py.MPI.COMM_WORLD

    if comm.rank == 0:
        geo = OCCGeometry(Box(Pnt(0, 0, 0), Pnt(1, 1, 1)))
        m = geo.GenerateMesh(maxh=0.3)
        m.Save("test_ed_segments")

    comm.Barrier()

    mesh = netgen.meshing.Mesh(3, comm)
    mesh.Load("test_ed_segments.vol.gz")

    if comm.rank == 0:
        assert mesh.ne == 0
    else:
        for seg in mesh.Elements1D():
            assert seg.index >= 0, f"segment has index={seg.index}"


@pytest.mark.mpi
@pytest.mark.skipif(not has_occ, reason="OCC not available")
def test_mpi_edge_descriptors_distributed():
    """Check that non-root ranks have edge descriptors with valid edgenr after MPI distribution."""
    comm = mpi4py.MPI.COMM_WORLD

    if comm.rank == 0:
        geo = OCCGeometry(Box(Pnt(0, 0, 0), Pnt(1, 1, 1)))
        m = geo.GenerateMesh(maxh=0.3)
        m.Save("test_ed_distributed")

    comm.Barrier()

    mesh = netgen.meshing.Mesh(3, comm)
    mesh.Load("test_ed_distributed.vol.gz")

    if comm.rank > 0:
        ned = mesh.GetNED()
        assert ned > 0, "non-root rank should have edge descriptors"
        for i in range(ned):
            ed = mesh.EdgeDescriptor(i)
            assert ed.edgenr > 0, f"ED[{i}] edgenr={ed.edgenr} should be > 0"


@pytest.mark.mpi
@pytest.mark.skipif(not has_occ, reason="OCC not available")
def test_mpi_cd2names_survive_distribution():
    """Check that named edges (CD2 names) survive MPI distribution."""
    comm = mpi4py.MPI.COMM_WORLD

    if comm.rank == 0:
        box = Box(Pnt(0, 0, 0), Pnt(1, 1, 1))
        box.edges[0].name = "test_edge"
        geo = OCCGeometry(box)
        m = geo.GenerateMesh(maxh=0.3)
        m.Save("test_ed_cd2names")

    comm.Barrier()

    mesh = netgen.meshing.Mesh(3, comm)
    mesh.Load("test_ed_cd2names.vol.gz")

    if comm.rank > 0:
        found_named = False
        for i in range(mesh.GetNCD2Names()):
            name = mesh.GetCD2Name(i)
            if name and name != "default":
                found_named = True
                break
        assert found_named, "No non-default CD2 names found after distribution"