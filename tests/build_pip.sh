set -e
ulimit -n 1024000 # lower open file limit, seems to affect performance
yum -y update
yum -y install ninja-build fontconfig-devel tk-devel tcl-devel libXmu-devel mesa-libGLU-devel ccache mpich-3.2-devel openmpi-devel

rm -rf wheelhouse
export NETGEN_CCACHE=1

/opt/python/cp39-cp39/bin/python tests/fix_auditwheel_policy.py

for pyversion in 312
do
    export PYDIR="/opt/python/cp${pyversion}-cp${pyversion}/bin"
    echo $PYDIR
    $PYDIR/pip install -U pytest-check numpy wheel scikit-build pybind11-stubgen

    rm -rf _skbuild
    NETGEN_ARCH=avx2 $PYDIR/pip wheel .
    auditwheel repair netgen_mesher*-cp${pyversion}-*.whl
    rm netgen_mesher-*.whl

    $PYDIR/pip install wheelhouse/netgen_mesher*-cp${pyversion}-*.whl
    $PYDIR/python3 -c 'import netgen'
    #cd ../tests/pytest
    #$PYDIR/python3 -m pytest
done

$PYDIR/pip install -U twine
$PYDIR/twine upload wheelhouse/*manylinux*.whl
