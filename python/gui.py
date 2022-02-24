import netgen

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
    win.tk.eval('load "'+netgen._netgen_lib_dir.replace('\\','/')+'/libgui[info sharedlibextension]" gui')

    if config.is_python_package and 'darwin' in sys.platform:
        # libngsolve and other libraries are installed into netgen python dir to keep relative installation paths, but tcl won't find them there automatically
        netgen_dir = os.path.abspath(os.path.dirname(netgen.__file__))
        win.tk.eval(f'set netgen_library_dir {netgen_dir}')

    win.tk.eval( netgen.libngpy._meshing._ngscript)

    
if not netgen.libngpy._meshing._netgen_executable_started:
    import os
    if not "NETGEN_DOCUMENTATION_RST_FORMAT" in os.environ:
        StartGUI()

def Snapshot(w,h, filename=None):
    netgen.Redraw(blocking=True)
    import numpy
    image = netgen.libngpy.Snapshot(w, h)
    image = numpy.array(image, dtype=numpy.uint8).reshape(h, w, 3)
    image = image[::-1,:,:]
    if filename:
        import PIL.Image
        im = PIL.Image.fromarray(image)
        im.save(filename)
    return image
