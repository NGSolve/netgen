import os
import sys

from . import config
_netgen_bin_dir=os.path.realpath(os.path.join(os.path.dirname(__file__),'..',config.NETGEN_PYTHON_RPATH_BIN))
_netgen_lib_dir=os.path.realpath(os.path.join(os.path.dirname(__file__),'..',config.NETGEN_PYTHON_RPATH))

if sys.platform.startswith('win'):
    if sys.version >= '3.8':
        os.add_dll_directory(_netgen_bin_dir)
    else:
        os.environ['PATH'] += ';'+_netgen_bin_dir

del sys
del os

from . import libngpy

from netgen.libngpy._meshing import _Redraw

def Redraw(*args, **kwargs):
    return _Redraw(*args, **kwargs)
