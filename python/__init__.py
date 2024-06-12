import os
import sys
import importlib.metadata
from pathlib import Path

from . import config
_netgen_bin_dir=os.path.realpath(os.path.join(os.path.dirname(__file__),'..',config.NETGEN_PYTHON_RPATH_BIN))
_netgen_lib_dir=os.path.realpath(os.path.join(os.path.dirname(__file__),'..',config.NETGEN_PYTHON_RPATH))

def _add_shared_lib_path(p: Path):
    print("Adding shared search path", f.locate().parent)
    if sys.platform.startswith('win'):
        os.add_dll_directory(p)
    elif 'linux' in sys.platform:
        if str(p) not in os.environ['LD_LIBRARY_PATH']:
            os.environ['LD_LIBRARY_PATH'] = str(p) + ':' + os.environ['LD_LIBRARY_PATH']
    elif 'darwin' in sys.platform:
        if str(p) not in os.environ['DYLD_LIBRARY_PATH']:
            os.environ['DYLD_LIBRARY_PATH'] = str(p) + ':' + os.environ['DYLD_LIBRARY_PATH']
try:
    importlib.metadata.metadata('netgen-occt')
    for f in importlib.metadata.files('netgen-occt'):
        if f.match('*libTKernel*'):
            print("Adding shared search path", f.locate().parent)
            _add_shared_lib_path(f.locate().parent)
except importlib.metadata.PackageNotFoundError:
    pass

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
    os.environ['PATH'] += ';'+_netgen_bin_dir

del sys
del os

from pyngcore import Timer
from . import libngpy

from netgen.libngpy._meshing import _Redraw

def Redraw(*args, **kwargs):
    return _Redraw(*args, **kwargs)

def TimeFunction(func, name=None):
    name = name or func.__qualname__
    timer = Timer(name)
    def retfunc(*args,**kwargs):
        with timer:
            ret = func(*args, **kwargs)
        return ret
    return retfunc

