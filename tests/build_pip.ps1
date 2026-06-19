$ErrorActionPreference = "Stop"

# Build and upload the Windows wheels for all supported CPython versions.
#
# We build directly against the python.org installations (C:\Python3XX) rather
# than via cibuildwheel: Netgen's GUI needs Tcl/Tk + the tcl/tk stub libraries
# and tkinter, none of which ship in the nuget CPython packages cibuildwheel
# uses on Windows. The python.org installers include all of them. We still use
# scikit-build-core (PEP 517) for the build and delvewheel for the repair step.

# Run a native command and fail the whole script (and the CI job) if it returns
# a non-zero exit code. PowerShell's $ErrorActionPreference does NOT do this for
# external programs, so it has to be checked explicitly.
function Invoke-Step([scriptblock]$Block) {
    & $Block
    if ($LASTEXITCODE -ne 0) {
        throw "Command failed with exit code ${LASTEXITCODE}: $Block"
    }
}

if (test-path wheelhouse) {
    cmd.exe /c rd /s /q wheelhouse
}

$env:NETGEN_CCACHE = "1"
$env:NETGEN_ARCH = "avx2"

$pythons = @(
    "C:\Python314",
    "C:\Python313",
    "C:\Python312",
    "C:\Python311",
    "C:\Python310"
)

foreach ($pydir in $pythons) {
    $python = Join-Path $pydir "python.exe"
    Invoke-Step { & $python --version }
    Invoke-Step { & $python -m pip install -U pip requests packaging }
    # skip this version if it is already published on PyPI (non-zero = skip)
    & $python tests\utils.py --check-pip
    if ($LASTEXITCODE -ne 0) { continue }

    # The runner's system site-packages is shared and occasionally contains
    # packages with broken/old metadata (e.g. 'bleach', 'requests-toolbelt')
    # that newer pip/urllib3 refuse to process. Remove the known offenders so
    # fresh, valid versions get pulled in below.
    & $python -m pip uninstall -y bleach requests-toolbelt

    # Build-critical deps only (twine is installed separately before upload, so a
    # twine dependency hiccup can never block the build itself).
    Invoke-Step { & $python -m pip install -U build delvewheel "scikit-build-core>=0.10" pybind11-stubgen netgen-occt==7.8.1 netgen-occt-devel==7.8.1 ninja cmake }

    if (test-path dist) { cmd.exe /c rd /s /q dist }

    # Point CMake at the Python prefix for Tcl/Tk discovery (as the old setup.py did)
    $prefix = $pydir -replace '\\', '/'
    $env:CMAKE_ARGS = "-DCMAKE_PREFIX_PATH=$prefix -DPython3_ROOT_DIR=$prefix"

    Invoke-Step { & $python -m build --wheel --no-isolation --outdir dist . }

    # Repair. OCC DLLs (TK*.dll) are shipped by the separate netgen-occt wheel and
    # pre-loaded at run time, so exclude them from vendoring.
    New-Item -ItemType Directory -Force -Path wheelhouse | Out-Null
    Get-ChildItem dist\*.whl | ForEach-Object {
        Invoke-Step { & $python -m delvewheel repair --ignore-existing --exclude "TK*.dll" -w wheelhouse $_.FullName }
    }

    # Upgrade twine together with the transitive deps that tend to be stale on the
    # runner (requests-toolbelt importing a urllib3 module removed in urllib3 2.x).
    Invoke-Step { & $python -m pip install -U twine requests-toolbelt urllib3 }
    Invoke-Step { & $python -m twine upload --skip-existing wheelhouse\*.whl }
}
