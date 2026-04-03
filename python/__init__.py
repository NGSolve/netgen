import os
import sys
from pathlib import Path

from . import config
_netgen_bin_dir=os.path.realpath(os.path.join(os.path.dirname(__file__),'..',config.NETGEN_PYTHON_RPATH_BIN))
_netgen_lib_dir=os.path.realpath(os.path.join(os.path.dirname(__file__),'..',config.NETGEN_PYTHON_RPATH))

def load_occ_libs():
    try:
        try:
            import importlib.metadata as metadata
        except ImportError:
            import importlib_metadata as metadata
        import ctypes
        metadata.metadata('netgen-occt')
        lib_names = [
            "TKOffset",
            "TKFillet",
            "TKDEIGES",
            "TKBool",
            "TKDESTEP",
            "TKXSBase",
            "TKDESTL",
            "TKXCAF",
            "TKVCAF",
            "TKCAF",
            "TKBO",
            "TKPrim",
            "TKLCAF",
            "TKCDF",
            "TKV3d",
            "TKHLR",
            "TKMesh",
            "TKService",
            "TKShHealing",
            "TKTopAlgo",
            "TKGeomAlgo",
            "TKBRep",
            "TKGeomBase",
            "TKG3d",
            "TKG2d",
            "TKMath",
            "TKDE",
            "TKernel",
        ]
        lib_names.reverse()
        lib_paths = {}
        for f in metadata.files('netgen-occt'):
            if f.match('*libTK*') or f.match("*.dll"):
                p = f.locate()
                name = p.name.split('.')[0].lower().replace("lib","")
                lib_paths[name] = str(p)
        for lib_name in lib_names:
            p = lib_paths[lib_name.lower()]
            ctypes.CDLL(p, mode=ctypes.RTLD_GLOBAL)

    except metadata.PackageNotFoundError:
        pass

load_occ_libs()

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

