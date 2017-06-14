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

if not netgen.libngpy._meshing._netgen_executable_started:
    StartGUI()
