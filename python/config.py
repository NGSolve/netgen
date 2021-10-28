def _cmake_to_bool(s):
    return s.upper() not in ['', '0','FALSE','OFF','N','NO','IGNORE','NOTFOUND']

is_python_package    = _cmake_to_bool("@SKBUILD@")

BUILD_FOR_CONDA     = _cmake_to_bool("@BUILD_FOR_CONDA@")
BUILD_STUB_FILES    = _cmake_to_bool("@BUILD_STUB_FILES@")
CHECK_RANGE         = _cmake_to_bool("@CHECK_RANGE@")
DEBUG_LOG           = _cmake_to_bool("@DEBUG_LOG@")
ENABLE_CPP_CORE_GUIDELINES_CHECK  = _cmake_to_bool("@ENABLE_CPP_CORE_GUIDELINES_CHECK@")
ENABLE_UNIT_TESTS   = _cmake_to_bool("@ENABLE_UNIT_TESTS@")
INSTALL_PROFILES    = _cmake_to_bool("@INSTALL_PROFILES@")
INTEL_MIC           = _cmake_to_bool("@INTEL_MIC@")
TRACE_MEMORY        = _cmake_to_bool("@TRACE_MEMORY@")
USE_CCACHE          = _cmake_to_bool("@USE_CCACHE@")
USE_CGNS            = _cmake_to_bool("@USE_CGNS@")
USE_GUI             = _cmake_to_bool("@USE_GUI@")
USE_INTERNAL_TCL    = _cmake_to_bool("@USE_INTERNAL_TCL@")
USE_JPEG            = _cmake_to_bool("@USE_JPEG@")
USE_MPEG            = _cmake_to_bool("@USE_MPEG@")
USE_MPI             = _cmake_to_bool("@USE_MPI@")
USE_MPI4PY          = _cmake_to_bool("@USE_MPI4PY@")
USE_NATIVE_ARCH     = _cmake_to_bool("@USE_NATIVE_ARCH@")
USE_NUMA            = _cmake_to_bool("@USE_NUMA@")
USE_OCC             = _cmake_to_bool("@USE_OCC@")
USE_PYTHON          = _cmake_to_bool("@USE_PYTHON@")
USE_SPDLOG          = _cmake_to_bool("@USE_SPDLOG@")

CMAKE_INSTALL_PREFIX  = "@CMAKE_INSTALL_PREFIX@"
NG_INSTALL_DIR_PYTHON   = "@NG_INSTALL_DIR_PYTHON_DEFAULT@"
NG_INSTALL_DIR_BIN      = "@NG_INSTALL_DIR_BIN_DEFAULT@"
NG_INSTALL_DIR_LIB      = "@NG_INSTALL_DIR_LIB_DEFAULT@"
NG_INSTALL_DIR_INCLUDE  = "@NG_INSTALL_DIR_INCLUDE_DEFAULT@"
NG_INSTALL_DIR_CMAKE    = "@NG_INSTALL_DIR_CMAKE_DEFAULT@"
NG_INSTALL_DIR_RES      = "@NG_INSTALL_DIR_RES_DEFAULT@"

NETGEN_PYTHON_RPATH_BIN = "@NETGEN_PYTHON_RPATH_BIN@"
NETGEN_PYTHON_RPATH     = "@NETGEN_PYTHON_RPATH@"

NG_COMPILE_FLAGS           = "@NG_COMPILE_FLAGS@"
ngcore_compile_options     = "@ngcore_compile_options@"
ngcore_compile_definitions = "@ngcore_compile_definitions@"

NETGEN_VERSION = "@NETGEN_VERSION@"
NETGEN_VERSION_GIT = "@git_version_string@"
NETGEN_VERSION_PYTHON = "@NETGEN_VERSION_PYTHON@"

NETGEN_VERSION_MAJOR = "@NETGEN_VERSION_MAJOR@"
NETGEN_VERSION_MINOR = "@NETGEN_VERSION_MINOR@"
NETGEN_VERSION_TWEAK = "@NETGEN_VERSION_TWEAK@"
NETGEN_VERSION_PATCH = "@NETGEN_VERSION_PATCH@"
NETGEN_VERSION_HASH = "@NETGEN_VERSION_HASH@"

version = NETGEN_VERSION_GIT
