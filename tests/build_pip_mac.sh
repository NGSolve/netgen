set -e
rm -rf _skbuild dist

export PATH=/Applications/CMake.app/Contents/bin:$PATH
export NETGEN_CCACHE=1

export PYDIR=/Library/Frameworks/Python.framework/Versions/$1/bin
$PYDIR/python3 --version
$PYDIR/pip3 install --user numpy twine scikit-build wheel pybind11-stubgen

export CMAKE_OSX_ARCHITECTURES='arm64;x86_64'
$PYDIR/python3 setup.py bdist_wheel --plat-name macosx-10.15-universal2 -j10
$PYDIR/python3 -m twine upload dist/*.whl
