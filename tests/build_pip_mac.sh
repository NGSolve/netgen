set -e

rm -rf _skbuild dist

export PYDIR=/Library/Frameworks/Python.framework/Versions/$1/bin
export PATH=$PYDIR:/Applications/CMake.app/Contents/bin:$PATH
export NETGEN_CCACHE=1

$PYDIR/python3 --version
$PYDIR/python3 -m pip install packaging requests
$PYDIR/python3 tests/utils.py --check-pip || exit 0
$PYDIR/python3 -m pip install --user numpy twine scikit-build wheel pybind11-stubgen
$PYDIR/python3 -m pip install --user -U netgen-occt==7.8.1 netgen-occt-devel==7.8.1

export CMAKE_OSX_ARCHITECTURES='arm64;x86_64'
export NETGEN_ARCH='avx2'
$PYDIR/python3 setup.py bdist_wheel --plat-name macosx-10.15-universal2 -j10
$PYDIR/python3 -m twine upload dist/*.whl
