set -e

# Build and upload the macOS wheels for all supported CPython versions using
# cibuildwheel (it downloads the CPython runtimes and runs delocate itself).
# The build configuration lives in pyproject.toml ([tool.cibuildwheel.*]).
rm -rf wheelhouse dist

export PATH=/Applications/CMake.app/Contents/bin:$PATH

python3 -m pip install -U pip cibuildwheel twine
python3 -m cibuildwheel --platform macos --output-dir wheelhouse
python3 -m twine upload --skip-existing wheelhouse/*.whl
