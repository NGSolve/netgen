import pytest
import netgen.meshing

MPI = pytest.importorskip("mpi4py.MPI")
pytestmark = pytest.mark.skipif(MPI.COMM_WORLD.size == 1, reason="needs more than one MPI rank")

def test_mpi4py():
    comm = MPI.COMM_WORLD

    if comm.rank==0:
        from netgen.csg import unit_cube
        m = unit_cube.GenerateMesh(maxh=0.1)
        m.Save("mpimesh")

    comm.Barrier()

    mesh = netgen.meshing.Mesh(3, comm)
    mesh.Load("mpimesh.vol.gz")

    assert mesh.ne > 0   # every rank, including 0, holds a part
