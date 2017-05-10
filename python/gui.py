from tkinter import Tk

from . import _netgen_lib_dir
from . import _netgen_bin_dir

win = Tk()
import os
dir_path = os.path.dirname(os.path.realpath(__file__))
win.tk.eval("source "+os.path.realpath(os.path.join(_netgen_bin_dir, 'ng.tcl')).replace('\\','/'))

# %gui tk
