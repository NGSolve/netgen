import math
import numpy as np
from time import time
import os

import webgui_jupyter_widgets
from webgui_jupyter_widgets import BaseWebGuiScene, encodeData, WebGuiDocuWidget
import webgui_jupyter_widgets.widget as wg

from packaging.version import parse

if parse(webgui_jupyter_widgets.__version__) >= parse("0.2.18"):
    _default_width = None
    _default_height = None
else:
    _default_width = "100%"
    _default_height = "50vh"


_registered_draw_types = {}


def register_draw_type(*types):
    def inner(func):
        for typ in types:
            _registered_draw_types[typ] = func

    return inner


class WebGLScene(BaseWebGuiScene):
    def __init__(self, obj, args=[], kwargs={}):
        self.obj = obj
        self.args = args
        self.kwargs = kwargs
        self.encoding = "b64"

    def GetData(self, set_minmax=True):
        self.kwargs["encoding"] = self.encoding
        typ = type(self.obj)
        d = None
        if type(self.obj) in _registered_draw_types:
            d = _registered_draw_types[typ](self.obj, self.args, self.kwargs)
        else:
            import inspect

            for t in inspect.getmro(typ):
                if t in _registered_draw_types:
                    print("check type", t)
                    d = _registered_draw_types[t](self.obj, self.args, self.kwargs)
                    break
        if d is None and hasattr(self.obj, "_webgui_data"):
            d = self.obj._webgui_data()
            bp = d["Bezier_trig_points"]
            for i in range(len(bp)):
                bp[i] = encodeData(np.array(bp[i], dtype=np.float32))

            ep = d["edges"]
            for i in range(len(ep)):
                ep[i] = encodeData(np.array(ep[i], dtype=np.float32))

        if d is None:
            raise RuntimeError(f"Cannot draw object of type {typ}")

        args = self.args
        kwargs = self.kwargs
        if "clipping" in kwargs:
            clipping = kwargs["clipping"]
            d["clipping"] = True
            if isinstance(clipping, dict):
                allowed_args = ("x", "y", "z", "dist", "function", "pnt", "vec")
                if "vec" in clipping:
                    vec = clipping["vec"]
                    clipping["x"] = vec[0]
                    clipping["y"] = vec[1]
                    clipping["z"] = vec[2]
                if "pnt" in clipping:
                    d["mesh_center"] = list(clipping["pnt"])
                for name, val in clipping.items():
                    if not (name in allowed_args):
                        raise Exception(
                            "Only {} allowed as arguments for clipping!".format(
                                ", ".join(allowed_args)
                            )
                        )
                    d["clipping_" + name] = val

        if "on_init" in kwargs:
            d["on_init"] = kwargs["on_init"]

        if set_minmax:
            if "min" in kwargs:
                d["funcmin"] = kwargs["min"]
            if "max" in kwargs:
                d["funcmax"] = kwargs["max"]
            d["autoscale"] = kwargs["autoscale"]

        if "vectors" in kwargs:
            d["vectors"] = True
            if isinstance(kwargs["vectors"], dict):
                for name, val in kwargs["vectors"].items():
                    if not (name in ("grid_size", "offset")):
                        raise Exception(
                            'Only "grid_size" and "offset" allowed as arguments for vectors!'
                        )
                    d["vectors_" + name] = val

        if "eval_function" in kwargs:
            d["user_eval_function"] = kwargs["eval_function"]

        # see shaders/utils.h for value explanation (function_mode)
        if "eval_" in kwargs:
            eval_ = kwargs["eval"]
            if isinstance(eval_, int):
                d["eval"] = eval_
            elif eval_ == "norm":
                d["eval"] = 3
            elif eval_ == "real":
                d["eval"] = 5
            elif eval_ == "imag":
                d["eval"] = 6

        return d


bezier_trig_trafos = {}  # cache trafos for different orders

# def Draw(shape, clipping=None, js_code=None, filename=""):
#     # todo: also handle occ geometry, list of shapes, etc.

#     scene = WebGLScene(shape, clipping=clipping, on_init=js_code)

#     if wg._IN_IPYTHON:
#         if wg._IN_GOOGLE_COLAB:
#             from IPython.display import display, HTML
#             html = scene.GenerateHTML()
#             display(HTML(html))
#         else:
#             # render scene using widgets.DOMWidget
#             scene.Draw()
#             return scene
#     else:
#         if filename:
#             scene.GenerateHTML(filename=filename)
#         return scene


def _get_draw_default_args():
    return dict(
        name="function",
        order=2,
        draw_vol=True,
        draw_surf=True,
        autoscale=True,
        deformation=False,
        interpolate_multidim=False,
        animate=False,
        objects=[],
        nodal_p1=False,
        settings={},
        width=_default_width,
        height=_default_height,
    )


def Draw(obj, *args, **kwargs):
    kwargs_with_defaults = _get_draw_default_args()
    kwargs_with_defaults.update(kwargs)

    scene = WebGLScene(obj, args, kwargs_with_defaults)
    if wg._IN_IPYTHON:
        if wg._IN_GOOGLE_COLAB:
            from IPython.display import display, HTML

            html = scene.GenerateHTML()
            display(HTML(html))
        else:
            import webgui_jupyter_widgets as wjw
            from packaging.version import parse

            # render scene using widgets.DOMWidget
            if parse(wjw.__version__) < parse("0.2.15"):
                scene.Draw()
            else:
                scene.Draw(
                    kwargs_with_defaults["width"], kwargs_with_defaults["height"]
                )
            return scene
    else:
        if "filename" in kwargs_with_defaults:
            scene.GenerateHTML(filename=kwargs_with_defaults["filename"])
        return scene


def _DrawDocu(obj, *args, **kwargs):
    kwargs_with_defaults = _get_draw_default_args()
    kwargs_with_defaults.update(kwargs)
    scene = WebGLScene(obj, args, kwargs)
    import json

    docu_path = os.environ["NETGEN_DOCUMENTATION_OUT_DIR"]
    src_path = os.environ["NETGEN_DOCUMENTATION_SRC_DIR"]
    cwd_path = os.path.abspath(".")
    rel_path = os.path.relpath(".", src_path)
    path = os.path.join(docu_path, rel_path)

    if not os.path.exists(path):
        os.makedirs(path)
    counter_file = os.path.join(docu_path, ".counter")
    if os.path.exists(counter_file):
        file_counter = int(open(counter_file, "r").read()) + 1
    else:
        file_counter = 0

    open(counter_file, "w").write(str(file_counter))

    data_file = "render_data_{}.json".format(file_counter)
    data_file_abs = os.path.join(path, data_file)
    preview_file = "preview_{}.png".format(file_counter)
    preview_file_abs = os.path.join(path, preview_file)

    widget = WebGuiDocuWidget()
    widget.value = {"render_data": data_file, "preview": preview_file}
    scene.widget = widget
    data = scene.GetData()
    json.dump(data, open(data_file_abs, "w"))
    scene.MakeScreenshot(preview_file_abs, 1200, 600)
    scene.Redraw = lambda: None
    from IPython.display import display, HTML

    display(widget)
    return scene


if "NETGEN_DOCUMENTATION_SRC_DIR" in os.environ:
    # we are buiding the documentation, some things are handled differently:
    # 1) Draw() is generating a .png (using headless chromium via selenium) and a render_data.json
    #    to show a preview image and load the render_data only when requested by user
    # 2) return a NGSDocuWebGuiWidget instead of NGSWebGuiWidget implementing the preview/load on demand of webgui

    _Draw = Draw
    Draw = _DrawDocu
