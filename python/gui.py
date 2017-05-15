import netgen

def StartGUI():
    import os
    from tkinter import Tk

    from . import _netgen_lib_dir
    from . import _netgen_bin_dir

    global win
    win = Tk()
    dir_path = os.path.dirname(os.path.realpath(__file__))
    win.tk.eval("source "+os.path.realpath(os.path.join(_netgen_bin_dir, 'ng.tcl')).replace('\\','/'))

if not netgen.libngpy._meshing._netgen_executable_started:
    StartGUI()
