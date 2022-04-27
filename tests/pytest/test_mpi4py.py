import pytest
import netgen.meshing

mpi4py = pytest.importorskip("mpi4py")
_ = pytest.importorskip("pytest_mpi")

@pytest.mark.mpi
def test_mpi4py():
    comm = mpi4py.MPI.COMM_WORLD

    if comm.rank==0:
        from netgen.csg import unit_cube
        m = unit_cube.GenerateMesh(maxh=0.1)
        m.Save("mpimesh")

    comm.Barrier()

    mesh = netgen.meshing.Mesh(3, comm)
    mesh.Load("mpimesh.vol.gz")

    if comm.rank==0:
        assert mesh.ne==0
