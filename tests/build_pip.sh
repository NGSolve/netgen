set -e
yum -y update
yum -y install ninja-build fontconfig-devel tk-devel tcl-devel libXmu-devel mesa-libGLU-devel ccache

rm -rf wheelhouse
export NETGEN_CCACHE=1

/opt/python/cp39-cp39/bin/python tests/fix_auditwheel_policy.py

for pyversion in 38 39 310
do
    export PYDIR="/opt/python/cp${pyversion}-cp${pyversion}/bin"
    echo $PYDIR
    $PYDIR/pip install -U pytest-check numpy wheel scikit-build pybind11-stubgen

    rm -rf _skbuild
    $PYDIR/pip wheel .
    auditwheel repair netgen_mesher*-cp${pyversion}-*.whl
    rm netgen_mesher-*.whl

    rm -rf _skbuild
    NETGEN_ARCH=avx2 $PYDIR/pip wheel .
    auditwheel repair netgen_mesher_avx2*-cp${pyversion}-*.whl
    rm netgen_mesher_avx2-*.whl

    $PYDIR/pip install wheelhouse/netgen_mesher*-cp${pyversion}-*.whl
    $PYDIR/python3 -c 'import netgen'
    #cd ../tests/pytest
    #$PYDIR/python3 -m pytest
done

$PYDIR/pip install -U twine
$PYDIR/twine upload wheelhouse/*manylinux*.whl
