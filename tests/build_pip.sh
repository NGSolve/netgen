set -e
ulimit -n 1024000 # lower open file limit, seems to affect performance
yum -y update
yum -y install ninja-build fontconfig-devel tk-devel tcl-devel libXmu-devel mesa-libGLU-devel ccache dpkg


curl http://ftp.de.debian.org/debian/pool/main/o/openmpi/libopenmpi-dev_4.1.6-13.3_amd64.deb -o openmpi-dev.deb
dpkg-deb -R openmpi-dev.deb /opt/openmpi
mv /opt/openmpi/usr/lib/x86_64-linux-gnu/openmpi/include /opt/openmpi/include


curl http://ftp.de.debian.org/debian/pool/main/m/mpich/libmpich-dev_4.2.0-5.1_amd64.deb -o mpich.deb
dpkg-deb -R mpich.deb /opt/mpich
mv /opt/mpich/usr/lib/x86_64-linux-gnu/mpich/include /opt/mpich/include


rm -rf wheelhouse
export NETGEN_CCACHE=1

/opt/python/cp39-cp39/bin/python tests/fix_auditwheel_policy.py

for pyversion in 38 39 310 311 312
do
    export PYDIR="/opt/python/cp${pyversion}-cp${pyversion}/bin"
    echo $PYDIR
    $PYDIR/pip install -U pytest-check numpy wheel scikit-build pybind11-stubgen
    $PYDIR/pip install -i https://pypi.anaconda.org/mpi4py/simple/ --pre mpi4py

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
