import os
import sys

_netgen_bin_dir=os.path.realpath(os.path.join(os.path.dirname(__file__),'..','@NETGEN_PYTHON_RPATH_BIN@'))
_netgen_lib_dir=os.path.realpath(os.path.join(os.path.dirname(__file__),'..','@NETGEN_PYTHON_RPATH@'))

if sys.platform.startswith('win'):
    if sys.version >= '3.8':
        os.add_dll_directory(_netgen_bin_dir)
    else:
        os.environ['PATH'] += ';'+_netgen_bin_dir

del sys
del os

from . import libngpy

def Redraw(*args, **kwargs):
    try:
        if libngpy.meshvis._Redraw(*args, **kwargs):
            import netgen
            import tkinter
            cnt = 0
            while(netgen.gui.win.tk.dooneevent(tkinter._tkinter.DONT_WAIT) and cnt < 100):
                cnt += 1
    except:
        pass


