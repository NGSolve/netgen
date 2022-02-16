import glob
import os
import sys

from skbuild import setup
import skbuild.cmaker
from subprocess import check_output
from distutils.sysconfig import get_python_lib;

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
version = '.dev'.join(version)

py_install_dir = get_python_lib(1,0,'').replace('\\','/')

name = "netgen-mesher"
arch = None
cmake_args = [
        f'-DNETGEN_VERSION_GIT={git_version}',
        f'-DNETGEN_VERSION_PYTHON={version}',
    ]

if 'NETGEN_ARCH' in os.environ:
  arch = os.environ['NETGEN_ARCH']
  name += "-"+arch
  flag = '/'+arch if 'win' in sys.platform else '-march=core-avx2'
  cmake_args += [f'-DNG_COMPILE_FLAGS={flag}']

if 'NETGEN_CCACHE' in os.environ:
  cmake_args += [f'-DUSE_CCACHE=ON']

packages = ['netgen', 'pyngcore']

if 'darwin' in sys.platform:
    cmake_args += [
        f'-DNG_INSTALL_DIR_LIB=netgen',
        f'-DNG_INSTALL_DIR_PYTHON=.',
        f'-DNG_INSTALL_DIR_CMAKE=lib/cmake',
        f'-DNG_INSTALL_DIR_BIN=bin',
    ]
elif 'win' in sys.platform:
    cmake_args += [
        '-A Win64',
        f'-DNG_INSTALL_DIR_BIN=netgen',
        f'-DNG_INSTALL_DIR_PYTHON=.',
        f'-DNG_INSTALL_DIR_LIB=Library/lib',
    ]
elif 'linux' in sys.platform:
    name_dir = name.replace('-','_')
    cmake_args += [
        f'-DNG_INSTALL_DIR_LIB={py_install_dir}/{name_dir}.libs',
        f'-DNG_INSTALL_DIR_BIN=bin',
    ]
    packages = []

cmake_args += [
        '-DUSE_SUPERBUILD:BOOL=ON',
        '-DUSE_CCACHE:BOOL=ON',
        '-DUSE_GUI=ON',
        '-DUSE_NATIVE_ARCH=OFF',
        '-DNG_INSTALL_DIR_INCLUDE=include/netgen',
        '-DBUILD_ZLIB=ON',
        '-DBUILD_OCC=ON',
        '-DUSE_OCC=ON',
        '-DBUILD_FOR_CONDA=ON',
        f'-DNETGEN_PYTHON_PACKAGE_NAME={name}',
        '-DBUILD_STUB_FILES=OFF',
]

if 'PYDIR' in os.environ:
    cmake_args += [f'-DCMAKE_PREFIX_PATH={os.environ["PYDIR"]}']

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
