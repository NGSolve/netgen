import os
import sys

from . import config
_netgen_bin_dir=os.path.realpath(os.path.join(os.path.dirname(__file__),'..',config.NETGEN_PYTHON_RPATH_BIN))
_netgen_lib_dir=os.path.realpath(os.path.join(os.path.dirname(__file__),'..',config.NETGEN_PYTHON_RPATH))

__diagnostics_template = """
Netgen diagnostics:
    sys.platform:          {sys.platform}
    sys.executable:        {sys.executable}
    sys.version:           {sys.version}
    Netgen python version: {config.PYTHON_VERSION}
    Netgen path            {__file__}
    Netgen config          {config.__file__}
    Netgen version         {config.NETGEN_VERSION}
    sys.path: {sys.path}
"""

def _get_diagnostics():
    return __diagnostics_template.format(sys=sys, config=config, __file__=__file__)

# compare compile-time and run-time python version
def _check_python_version():
    sys_version = f"{sys.version_info.major}.{sys.version_info.minor}"
    compile_version = f"{config.PYTHON_VERSION_MAJOR}.{config.PYTHON_VERSION_MINOR}"

    if sys_version != compile_version:
        print(_get_diagnostics(), file=sys.stderr)
        raise RuntimeError(f"Python version mismatch: compile-time version is {compile_version}, run-time version is {sys_version}")

_check_python_version()

if sys.platform.startswith('win'):
    v = sys.version_info
    if v.major == 3 and v.minor >= 8:
        os.add_dll_directory(_netgen_bin_dir)
    else:
        os.environ['PATH'] += ';'+_netgen_bin_dir

del sys
del os

from . import libngpy

from netgen.libngpy._meshing import _Redraw

def Redraw(*args, **kwargs):
    return _Redraw(*args, **kwargs)
