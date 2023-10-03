import glob
import os.path
import sys
import pathlib
import sysconfig

from skbuild import setup
import skbuild.cmaker
from subprocess import check_output

setup_requires = []

def install_filter(cmake_manifest):
    print(cmake_manifest)
    return cmake_manifest

def _patched_parse_manifests(self):
    paths = \
        glob.glob(os.path.join(skbuild.cmaker.CMAKE_BUILD_DIR(), "netgen", "install_manifest*.txt"))
    try:
        return [self._parse_manifest(path) for path in paths][0]
    except IndexError:
        return []
   
# we are using the netgen superbuild (to download and build some dependencies)
# patch the parse_manifests function to point to the actual netgen cmake project within the superbuild
skbuild.cmaker.CMaker._parse_manifests = _patched_parse_manifests

git_version = check_output(['git', 'describe', '--tags']).decode('utf-8').strip()
version = git_version[1:].split('-')
if len(version)>2:
    version = version[:2]
if len(version)>1:
    version = '.post'.join(version) + '.dev'
else:
    version = version[0]

py_install_dir = os.path.relpath(sysconfig.get_path('platlib'), sysconfig.get_path('data')).replace('\\','/')

name = "netgen-mesher"
arch = None
cmake_args = [
        f'-DNETGEN_VERSION_GIT={git_version}',
        f'-DNETGEN_VERSION_PYTHON={version}',
    ]

if 'NETGEN_ARCH' in os.environ and os.environ['NETGEN_ARCH'] == 'avx2':
    # build for avx2 architecture
    if 'darwin' in sys.platform:
        flag = "'-Xarch_x86_64;-march=core-avx2'"
    elif 'win' in sys.platform:
        flag = '/AVX2'
    else:
        flag = '-march=core-avx2'
    cmake_args += [f'-DNG_COMPILE_FLAGS={flag}']

if 'NETGEN_CCACHE' in os.environ:
  cmake_args += [f'-DUSE_CCACHE=ON']

packages = ['netgen', 'pyngcore']

if 'darwin' in sys.platform:
    cmake_args += [
        '-DNG_INSTALL_DIR_LIB=netgen',
        '-DNG_INSTALL_DIR_PYTHON=.',
        '-DNG_INSTALL_DIR_BIN=bin',
        '-DNG_INSTALL_DIR_CMAKE=netgen/cmake',
        '-DNG_INSTALL_DIR_INCLUDE=netgen/include',
        '-DNG_INSTALL_DIR_RES=share',
    ]
elif 'win' in sys.platform:
    cmake_args += [
        '-A Win64',
        '-DNG_INSTALL_DIR_BIN=netgen',
        '-DNG_INSTALL_DIR_PYTHON=.',
        '-DNG_INSTALL_DIR_LIB=netgen/lib',
        '-DNG_INSTALL_DIR_CMAKE=netgen/cmake',
        '-DNG_INSTALL_DIR_INCLUDE=netgen/include',
    ]
elif 'linux' in sys.platform:
    name_dir = name.replace('-','_')
    cmake_args += [
        f'-DNG_INSTALL_DIR_LIB={py_install_dir}/{name_dir}.libs',
        '-DNG_INSTALL_DIR_BIN=bin',
        '-DNG_INSTALL_DIR_INCLUDE=include/netgen',
        '-DTCL_INCLUDE_PATH=/usr/include',
        '-DTK_INCLUDE_PATH=/usr/include',
    ]
    packages = []

cmake_args += [
        '-DUSE_SUPERBUILD:BOOL=ON',
        '-DUSE_CCACHE:BOOL=ON',
        '-DUSE_GUI=ON',
        '-DUSE_NATIVE_ARCH=OFF',
        '-DBUILD_ZLIB=ON',
        '-DBUILD_OCC=ON',
        '-DUSE_OCC=ON',
        '-DBUILD_FOR_CONDA=ON',
        f'-DNETGEN_PYTHON_PACKAGE_NAME={name}',
        '-DBUILD_STUB_FILES=OFF',
]

pyprefix = pathlib.Path(sys.prefix).as_posix()
cmake_args += [f'-DCMAKE_PREFIX_PATH={pyprefix}', f'-DPython3_ROOT_DIR={pyprefix}']

setup(
    name=name,
    version=version,
    description="Netgen",
    author='The Netgen team',
    license="LGPL2.1",
    packages=packages,
    #package_dir={'netgen': 'python'},
    tests_require=['pytest'],
    #include_package_data=True,
    cmake_process_manifest_hook=install_filter,
    cmake_args = cmake_args,
    setup_requires=setup_requires,
    entry_points={
    'console_scripts': [
        'netgen = netgen.__main__:main',
    ],
},
)
