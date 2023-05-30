import importlib.util, threading, sys, os

def _py_handler(f):
    spec = importlib.util.spec_from_file_location(os.path.basename(f), os.path.abspath(f))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def _geo_handler(f):
    from netgen.csg import CSGeometry
    print("load", f)
    geo = CSGeometry(f)
    geo.Draw()

def _step_handler(f):
    from netgen.occ import OCCGeometry
    print("load", f)
    geo = OCCGeometry(f)
    geo.Draw()

def _stl_handler(f):
    from netgen.stl import STLGeometry
    print("load", f)
    geo = STLGeometry(f)
    geo.Draw()

_file_handler = {}
_file_handler['.py'] = _py_handler
_file_handler['.geo'] = _geo_handler
_file_handler['.step'] = _step_handler
_file_handler['.stl'] = _stl_handler

def handle_arguments():
    import __main__
    import os.path
    argv = sys.argv
    if len(argv)>1:
      _, ext = os.path.splitext(argv[1])
      if ext in _file_handler:
          _file_handler[ext](argv[1])

def main():
    import netgen
    # Use Redraw without event handling
    netgen.Redraw = netgen._Redraw

    from .gui import win
    th = threading.Thread(target=handle_arguments)
    th.start()
    win.tk.mainloop()

if __name__ == "__main__":
    sys.exit(main())
