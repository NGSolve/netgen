import netgen
import pytest

def test_gui():
    try:
        from tkinter import Tk
        win = Tk()
    except:
        pytest.skip("can't create a window")
    import netgen.gui

