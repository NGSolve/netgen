set -e
ulimit -n 1024000 # lower open file limit, seems to affect performance
yum -y update
yum -y install ninja-build fontconfig-devel tk-devel tcl-devel libXmu-devel mesa-libGLU-devel openmpi-devel mpich-devel

mkdir -p /opt/openmpi/include /opt/mpich/include
cp -a /usr/include/openmpi-x86_64/* /opt/openmpi/include/
cp -a /usr/include/mpich-x86_64/* /opt/mpich/include/

curl https://dl.fedoraproject.org/pub/epel/8/Everything/x86_64/Packages/c/ccache-3.7.7-1.el8.x86_64.rpm -o ccache.rpm
dnf -y install ccache.rpm

rm -rf wheelhouse
export NETGEN_CCACHE=1

for pyversion in 313 312 311 310 39 38
do
    export PYDIR="/opt/python/cp${pyversion}-cp${pyversion}/bin"
    echo $PYDIR
    $PYDIR/pip install requests packaging
    $PYDIR/python3 ./tests/utils.py --check-pip || continue
    $PYDIR/pip install -U pytest-check numpy wheel scikit-build pybind11-stubgen netgen-occt==7.8.1 netgen-occt-devel==7.8.1
    $PYDIR/pip install -i https://pypi.anaconda.org/mpi4py/simple/ --pre mpi4py

    rm -rf _skbuild
    NETGEN_ARCH=avx2 $PYDIR/pip wheel .
    mkdir -p wheelhouse
    rename linux_x86_64 manylinux_2_17_x86_64.manylinux2014_x86_64 netgen_mesher*-cp${pyversion}-*.whl
    mv netgen_mesher*-cp${pyversion}-*.whl wheelhouse/

    $PYDIR/pip install wheelhouse/netgen_mesher*-cp${pyversion}-*.whl
    $PYDIR/python3 -c 'import netgen'
    $PYDIR/pip install -U twine
    $PYDIR/twine upload --skip-existing wheelhouse/netgen_mesher*-cp${pyversion}*manylinux*.whl
    #cd ../tests/pytest
    #$PYDIR/python3 -m pytest
done

