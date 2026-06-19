set -e
ulimit -n 1024000 # lower open file limit, seems to affect performance

# This script runs *inside* the manylinux build image (see pip_linux in
# .gitlab-ci.yml) and builds one wheel per supported CPython using
# scikit-build-core (PEP 517) + auditwheel. The MPI headers are staged where the
# CMake build expects them (Netgen is built with the run-time MPI wrapper).
mkdir -p /opt/openmpi/include /opt/mpich/include
cp -a /usr/include/openmpi-x86_64/* /opt/openmpi/include/
cp -a /usr/include/mpich-x86_64/* /opt/mpich/include/

rm -rf wheelhouse dist
export NETGEN_CCACHE=1
export NETGEN_ARCH=avx2
export NETGEN_MPI=ON

for pyversion in 314 313 312 311 310
do
    export PYDIR="/opt/python/cp${pyversion}-cp${pyversion}/bin"
    echo $PYDIR
    $PYDIR/pip install -U pip requests packaging
    # skip this version if it is already published on PyPI
    $PYDIR/python3 ./tests/utils.py --check-pip || continue
    # Build deps are installed directly into the (header-equipped) manylinux
    # interpreter; we build with --no-isolation so CMake's FindPython sees the
    # real Python headers (an isolated build venv has none).
    $PYDIR/pip install -U build auditwheel twine ninja cmake
    $PYDIR/pip install -U "scikit-build-core>=0.10" pybind11-stubgen netgen-occt==7.8.1 netgen-occt-devel==7.8.1

    # Point CMake at the Python prefix for Tcl/Tk discovery (as the old setup.py did)
    PYPREFIX=$($PYDIR/python3 -c 'import sys; print(sys.prefix)')
    export CMAKE_ARGS="-DCMAKE_PREFIX_PATH=${PYPREFIX} -DPython3_ROOT_DIR=${PYPREFIX}"

    rm -rf dist
    $PYDIR/python3 -m build --wheel --no-isolation --outdir dist .

    # Repair to a manylinux wheel. OCC (libTK*) is shipped by the separate
    # netgen-occt wheel and pre-loaded at run time, so exclude it from vendoring.
    mkdir -p wheelhouse
    $PYDIR/python3 -m auditwheel repair --exclude 'libTK*' -w wheelhouse dist/*.whl

    # ---- debug: show what ended up in the wheel ----
    echo "==================== wheel contents ===================="
    $PYDIR/python3 -m zipfile -l wheelhouse/netgen_mesher*-cp${pyversion}-*manylinux*.whl || true
    echo "========================================================"

    $PYDIR/pip install wheelhouse/netgen_mesher*-cp${pyversion}-*.whl

    # ---- debug: show the installed package tree and version info ----
    echo "================= installed package tree ==============="
    $PYDIR/python3 - <<'PYEOF'
import os, netgen
pkg = os.path.dirname(netgen.__file__)
print("netgen package dir:", pkg)
for root, _, files in os.walk(pkg):
    for f in sorted(files):
        print("  ", os.path.relpath(os.path.join(root, f), pkg))
print("netgen.__version__         =", getattr(netgen, "__version__", "<MISSING>"))
print("netgen.config.NETGEN_VERSION =", netgen.config.NETGEN_VERSION)
PYEOF
    echo "========================================================"

    $PYDIR/python3 -c 'import netgen; print(netgen.__version__)'
    $PYDIR/twine upload --skip-existing wheelhouse/netgen_mesher*-cp${pyversion}*manylinux*.whl
done
