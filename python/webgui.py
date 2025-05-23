import math
import numpy as np
from time import time
import os

try:
    import webgui_jupyter_widgets
    from webgui_jupyter_widgets import BaseWebGuiScene, WebGuiDocuWidget
    import webgui_jupyter_widgets.widget as wg
except ImportError:
    class BaseWebGuiScene:
        pass

    wg = None

def encodeData( data, dtype=None, encoding='b64' ):
    import numpy as np
    from base64 import b64encode
    dtype = dtype or data.dtype
    values = np.array(data.flatten(), dtype=dtype)
    if encoding=='b64':
        return b64encode(values).decode("ascii")
    elif encoding=='binary':
        return values.tobytes()
    else:
        raise RuntimeError("unknown encoding" + str(encoding))

from packaging.version import parse

import netgen.meshing as ng

if wg is not None and parse(webgui_jupyter_widgets.__version__) >= parse("0.2.18"):
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


_bernstein_cache = {}


def GetIBernsteinBasis(etype, order):
    if (etype, order) in _bernstein_cache:
        return _bernstein_cache[(etype, order)]
    bvals = None

    if etype == "segment":

        def Binomial(n, i):
            return math.factorial(n) / math.factorial(i) / math.factorial(n - i)

        def Bernstein(x, i, n):
            return Binomial(n, i) * x**i * (1 - x) ** (n - i)

        bvals = np.zeros(
            (order + 1, order + 1), dtype=float
        )  # .Matrix(order+1,order+1)
        for i in range(order + 1):
            for j in range(order + 1):
                bvals[i, j] = Bernstein(i / order, j, order)

    if etype == "trig":

        def BernsteinTrig(x, y, i, j, n):
            return (
                math.factorial(n)
                / math.factorial(i)
                / math.factorial(j)
                / math.factorial(n - i - j)
                * x**i
                * y**j
                * (1 - x - y) ** (n - i - j)
            )

        og = order
        ndtrig = int((og + 1) * (og + 2) / 2)
        bvals = np.zeros((ndtrig, ndtrig))
        ii = 0
        for ix in range(og + 1):
            for iy in range(og + 1 - ix):
                jj = 0
                for jx in range(og + 1):
                    for jy in range(og + 1 - jx):
                        bvals[ii, jj] = BernsteinTrig(ix / og, iy / og, jx, jy, og)
                        jj += 1
                ii += 1

    if bvals is None:
        raise RuntimeError(f"Unkown element type {etype}")

    ibvals = _bernstein_cache[(etype, order)] = np.linalg.inv(bvals)
    return ibvals


def GetWireframePoints(etype, order):
    n = order
    if etype == "trig":
        return np.array(
            [(i / n, 0) for i in range(n + 1)]
            + [(0, i / n) for i in range(n + 1)]
            + [(i / n, 1.0 - i / n) for i in range(n + 1)]
        )
    if etype == "quad":
        return np.array(
            [(i / n, 0) for i in range(n + 1)]
            + [(0, i / n) for i in range(n + 1)]
            + [(i / n, 1.0) for i in range(n + 1)]
            + [(1.0, i / n) for i in range(n + 1)]
        )

    raise RuntimeError(f"Unknown element type {etype}")


def GetElementPoints(etype, order):
    n = order
    if etype == "trig":
        return np.array(
            [(i / n, j / n) for j in range(n + 1) for i in range(n + 1 - j)]
        )
    if etype == "quad":
        return np.array(
            [(i / n, j / n) for j in range(n + 1) for i in range(n + 1 - j)]
            + [(1 - i / n, 1 - j / n) for j in range(n + 1) for i in range(n + 1 - j)]
        )

    raise RuntimeError(f"Unknown element type {etype}")


def MapBernstein(pnts, etype, order):
    """
    Maps function values at equidistant control points to the Bernstein basis function.
    Parameters:
       pnts (numpy.ndarray): The input control points with shape (number_of_elements, points_per_element, function_dimension)
            point_per_element must be a multiple of the basis size
       etype (str): Element type (currently ignored and trig assumed)
       order (int): Polynomial order

    Returns:
        numpy.ndarray: The mapped points with the shape (points_per_element, number_of_elements, function_dimension)
    """
    ibvals = GetIBernsteinBasis(etype, order)
    # for wireframe or subdivided elements, we have multiple point sets per element
    # so do a reshape to simulate more elements with correct number of control points per element instead
    if pnts.shape[1] != ibvals.shape[0]:
        pnts = pnts.reshape((-1, ibvals.shape[0], pnts.shape[2]))

    points = np.zeros(pnts.shape, dtype=np.float32).transpose(1, 0, 2)
    for i in range(points.shape[2]):
        points[:, :, i] = np.tensordot(ibvals, pnts[:, :, i], axes=(1, 1))
    return points


@register_draw_type(ng.Mesh)
def GetData(mesh, args, kwargs):
    d = {}
    d["gui_settings"] = kwargs["settings"]
    d["mesh_dim"] = mesh.dim

    pmin, pmax = mesh.bounding_box
    diag = pmax - pmin
    pmid = pmin + 0.5 * diag
    d["mesh_center"] = [pmid[i] for i in range(3)]
    d["mesh_radius"] = diag.Norm()

    d["funcdim"] = 0
    d["show_mesh"] = True
    d["draw_surf"] = True
    d["funcmin"] = 0.0
    d["funcmax"] = 1.0

    # Generate surface element data
    # webgui code assumes 4 scalar fields (x,y,z, mesh_index)
    # TODO: other element types than trigs
    order = kwargs["order"]
    refpts = GetElementPoints("trig", order)
    pnts = np.ndarray((len(mesh.Elements2D()), refpts.shape[0], 4))
    mesh.CalcElementMapping(refpts, pnts)

    # set mesh_index
    for i, el in enumerate(mesh.Elements2D()):
        pnts[i, :, 3] = el.index - 1
    fds = mesh.FaceDescriptors()
    d["colors"] = [fd.color +(fd.transparency,) for fd in fds]
    d["mesh_regions_2d"] = len(fds)
    d["names"] = [fd.bcname for fd in fds]

    d["Bezier_trig_points"] = MapBernstein(pnts, "trig", order)
    d["order2d"] = order

    # Generate wireframe data
    refpts = GetWireframePoints("trig", order)
    pnts = np.ndarray((len(mesh.Elements2D()), refpts.shape[0], 4))
    mesh.CalcElementMapping(refpts, pnts)
    d["Bezier_points"] = MapBernstein(pnts, "segment", order)
    d["show_wireframe"] = True

    # TODO: Generate edge data
    d["edges"] = []

    # encode data as b64
    for name in ["Bezier_trig_points", "edges", "Bezier_points"]:
        pnew = []
        for plist in d[name]:
            pnew.append(encodeData(np.array(plist, dtype=np.float32)))
        d[name] = pnew
    return d

class WebGLScene(BaseWebGuiScene):
    class Widget:
        def __init__(self):
            self.value = {}

    def __init__(self, obj, args=[], kwargs={}):
        self.widget = self.Widget()
        self.obj = obj
        self.args = args
        self.kwargs = kwargs
        self.encoding = "b64"

    def Redraw(self, *args, **kwargs):
        if args or kwargs:
            if 'show' not in kwargs:
                kwargs['show'] = False

            new_scene = Draw(*args, **kwargs)
            self.obj = new_scene.obj
            self.args = new_scene.args
            self.kwargs = new_scene.kwargs
        super().Redraw()

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

        if "js_code" in kwargs:
            d["on_init"] = kwargs["js_code"]

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

        if "fullscreen" in kwargs:
            d["fullscreen"] = kwargs["fullscreen"]
        if "gui_settings" not in d:
            d["gui_settings"] = self.kwargs["settings"]

        if "euler_angles" in kwargs:
            camera = d["gui_settings"].get("camera", {})
            camera["euler_angles"] = kwargs["euler_angles"]
            d["gui_settings"]['camera'] = camera

        d["objects"] = []
        for obj in kwargs["objects"]:
            if isinstance(obj, dict):
                d["objects"].append(obj)
            else:
                d["objects"].append(obj._GetWebguiData())

        if 'center' in kwargs:
            center = list(kwargs['center'])
            if len(center) == 2:
                center.append(0.)
            d["mesh_center"] = center

        if 'radius' in kwargs:
            d["mesh_radius"] = kwargs['radius']

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
        fullscreen=False,
        scale=1.0,
        width=_default_width,
        height=_default_height,
    )


def Draw(obj, *args, show=True, **kwargs):
    kwargs_with_defaults = _get_draw_default_args()
    kwargs_with_defaults.update(kwargs)

    scene = WebGLScene(obj, args, kwargs_with_defaults)
    if show and wg is not None and wg._IN_IPYTHON:
        if wg._IN_GOOGLE_COLAB:
            from IPython.display import display, HTML

            html = scene.GenerateHTML()
            display(HTML(html))
            return
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
    if "filename" in kwargs_with_defaults:
        scene.GenerateHTML(filename=kwargs_with_defaults["filename"])
    return scene

async def _MakeScreenshot(data, png_file, width=1200, height=600):
    """Uses playwright to make a screenshot of the given html file."""
    # pylint: disable=import-outside-toplevel
    from playwright.async_api import async_playwright
    from webgui_jupyter_widgets.html import GenerateHTML, getScreenshotHTML

    html_file = png_file + ".html"
    GenerateHTML(data, filename=html_file, template=getScreenshotHTML())

    async with async_playwright() as play:
        browser = await play.chromium.launch()
        page = await browser.new_page(viewport={"width": width, "height": height})
        await page.goto(f"file://{os.path.abspath(html_file)}")
        await page.screenshot(path=png_file)
        await browser.close()

def _DrawDocu(obj, *args, **kwargs):
    import json
    import asyncio

    kwargs_with_defaults = _get_draw_default_args()
    kwargs_with_defaults.update(kwargs)
    scene = WebGLScene(obj, args, kwargs_with_defaults)

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
    asyncio.run(_MakeScreenshot(data, preview_file_abs, 1200, 600))
    scene.Redraw = lambda: None
    from IPython.display import display, HTML

    display(widget)
    return scene


if "NETGEN_DOCUMENTATION_SRC_DIR" in os.environ:
    # use nest_asyncio to allow reentrant asyncio when executing jupyter notebooks
    import nest_asyncio
    nest_asyncio.apply()

    # we are buiding the documentation, some things are handled differently:
    # 1) Draw() is generating a .png (using headless chromium via selenium) and a render_data.json
    #    to show a preview image and load the render_data only when requested by user
    # 2) return a NGSDocuWebGuiWidget instead of NGSWebGuiWidget implementing the preview/load on demand of webgui

    _Draw = Draw
    Draw = _DrawDocu
