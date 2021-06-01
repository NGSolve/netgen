import netgen

def StartGUI():
    from tkinter import Tk

    global win
    win = Tk()
    win.tk.eval('lappend ::auto_path ' + netgen._netgen_lib_dir)
    win.tk.eval('lappend ::auto_path ' + netgen._netgen_bin_dir)
    # load with absolute path to avoid issues on MacOS
    win.tk.eval('load "'+netgen._netgen_lib_dir.replace('\\','/')+'/libgui[info sharedlibextension]" gui')
    win.tk.eval( netgen.libngpy._meshing._ngscript)

    try:
        from IPython import get_ipython
        ipython = get_ipython()
        ipython.magic('gui tk')
    except:
        pass
    
if not netgen.libngpy._meshing._netgen_executable_started:
    import os
    if not "NETGEN_DOCUMENTATION_RST_FORMAT" in os.environ:
        try:
            StartGUI()
        except:
            pass

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
