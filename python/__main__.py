import imp, threading, sys

def handle_arguments():
    import __main__
    argv = sys.argv
    if len(argv)>1 and argv[1].endswith(".py"):
      with open(argv[1]) as pyfile:
          imp.load_module('__main__', pyfile, argv[1], (".py", "r", imp.PY_SOURCE))

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
