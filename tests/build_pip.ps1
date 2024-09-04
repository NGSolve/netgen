if (test-path _skbuild) {
    cmd.exe /c rd /s /q _skbuild
}
if (test-path dist) {
    cmd.exe /c rd /s /q dist
}

$env:NETGEN_CCACHE = 1
$env:NETGEN_ARCH = 'avx2'

$pydir=$args[0]
& $pydir\python.exe --version
& $pydir\python.exe -m pip install packaging
& $pydir\python.exe tests\utils.py --check-pip
if ($LASTEXITCODE -ne 0) {
    exit 0
}
& $pydir\python.exe -m pip install scikit-build wheel numpy twine pybind11-stubgen
& $pydir\python.exe -m pip install --upgrade netgen-occt==7.8.1 netgen-occt-devel==7.8.1
& $pydir\python setup.py bdist_wheel -G"Visual Studio 16 2019"
& $pydir\python -m twine upload dist\*.whl
