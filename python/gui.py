import netgen

from . import libngguipy
from . import libngpy

def StartGUI():
    from tkinter import Tk
    from netgen import config
    import sys, os

    try:
        # the GUI tries to load ngsolve.tcl (which loads ngsolve shared libraries)
        # BUT might fail to load dependencies (like mkl), these are handled by the
        # ngsolve __init__.py script, so import ngsolve from python already here
        import ngsolve
    except:
        pass

    global win
    win = Tk()
    win.tk.eval('lappend ::auto_path ' + netgen._netgen_lib_dir)
    win.tk.eval('lappend ::auto_path ' + netgen._netgen_bin_dir)
    # load with absolute path to avoid issues on MacOS
    win.tk.eval('load "'+netgen._netgen_lib_dir.replace('\\','/')+'/libnggui[info sharedlibextension]" gui')

    if config.is_python_package and 'darwin' in sys.platform:
        # libngsolve and other libraries are installed into netgen python dir to keep relative installation paths, but tcl won't find them there automatically
        netgen_dir = os.path.abspath(os.path.dirname(netgen.__file__))
        win.tk.eval(f'set netgen_library_dir {netgen_dir}')

    win.tk.eval( netgen.libngpy._meshing._ngscript)

    try:
        from IPython import get_ipython
        ipython = get_ipython()
        ipython.magic('gui tk')
    except:
        pass

    def _Redraw(*args, **kwargs):
        if libngpy._meshing._Redraw(*args, **kwargs):
            import netgen
            import tkinter
            cnt = 0
            while(win.tk.dooneevent(tkinter._tkinter.DONT_WAIT) and cnt < 100):
                cnt += 1

    netgen._Redraw = _Redraw
    _Redraw(blocking=True)

    

def Snapshot(w,h, filename=None):
    netgen.Redraw(blocking=True)
    import numpy
    image = netgen.libngguipy.Snapshot(w, h)
    image = numpy.array(image, dtype=numpy.uint8).reshape(h, w, 3)
    image = image[::-1,:,:]
    if filename:
        import PIL.Image
        im = PIL.Image.fromarray(image)
        im.save(filename)
    return image

def run_pyfile(filename):
    with open(filename) as f:
        exec(f.read(), {})

if __name__ == "__main__":
    import sys, threading
    StartGUI()
    if(len(sys.argv) > 1):
        if(sys.argv[1].endswith(".py")):
            t = threading.Thread(target=run_pyfile, args=(sys.argv[1],),
                                 daemon=True)
            t.start()
    win.mainloop()
else:
    if not netgen.libngpy._meshing._netgen_executable_started:
        import os
        if not "NETGEN_DOCUMENTATION_RST_FORMAT" in os.environ:
            StartGUI()
