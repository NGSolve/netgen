import math
import numpy as np
from time import time
import os

import anywidget
import traitlets

try:
    from IPython import get_ipython
    _IN_IPYTHON = get_ipython() is not None
except ImportError:
    _IN_IPYTHON = False


class BaseWebGuiScene:
    pass

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

import netgen.meshing as ng

_default_width = "100%"
_default_height = "20rem"

_registered_draw_types = {}


def register_draw_type(*types):
    """Decorator to register a drawing handler for one or more types.

    The decorated function is called as ``func(obj, args, kwargs)`` where
    *obj* is the object passed to :func:`Draw`, *args* are its positional
    arguments, and *kwargs* is a dict of keyword arguments merged with the
    default Draw options.  It must return a dict of render data.

    Parameters
    ----------
    *types
        One or more types for which the handler should be invoked.

    Example
    -------
    >>> @register_draw_type(MyMesh)
    ... def draw_my_mesh(obj, args, kwargs):
    ...     return {"mesh_dim": 3, ...}
    """
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

class WebGuiWidget(anywidget.AnyWidget):
    _esm = """
    function loadWebgui() {
        // webgui's dist is a Vite IIFE build that assigns a global `webgui`
        // var; it has no module exports, so it must be loaded as a classic
        // <script> and read off of window rather than import()-ed.
        if (window.webgui) return Promise.resolve(window.webgui);
        if (!window._webgui_loading) {
            window._webgui_loading = new Promise((resolve, reject) => {
                const s = document.createElement("script");
                s.src = "https://cdn.jsdelivr.net/npm/webgui@0.2.39/dist/webgui.js";
                s.onload = () => resolve(window.webgui);
                s.onerror = reject;
                document.head.appendChild(s);
            });
        }
        return window._webgui_loading;
    }

    // dat.gui (the control panel webgui adds, top-right) injects its stylesheet
    // into the main document <head>. When this widget is hosted by MyST it lives
    // inside a shadow root, and head styles do not cross the shadow boundary, so
    // the panel renders unstyled. Copy the dat.gui stylesheet into our shadow
    // root so it applies. No-op in a normal Jupyter page (no shadow root).
    function adoptDatGuiStyles(el) {
        const root = el.getRootNode();
        if (typeof ShadowRoot === 'undefined' || !(root instanceof ShadowRoot)) return;
        if (root.querySelector('style[data-webgui-datgui]')) return;
        let copied = false;
        for (const style of document.querySelectorAll('head style')) {
            const css = style.textContent || '';
            if (/\\.dg[\\s.,{>:]/.test(css)) {
                const clone = document.createElement('style');
                clone.setAttribute('data-webgui-datgui', '');
                clone.textContent = css;
                root.appendChild(clone);
                copied = true;
            }
        }
        return copied;
    }

    // The render data is normally the `value` trait itself. When `value` is a
    // string it is instead a URL to a JSON file holding the data: the static
    // MyST build offloads the (large) render data to a sidecar file next to the
    // notebook to keep pages small (see scripts/widgets_to_directives.py). Fetch
    // it in that case; otherwise use the value as-is.
    function resolveRenderData(value) {
        if (typeof value === 'string') {
            return fetch(value).then((resp) => resp.json());
        }
        return Promise.resolve(value);
    }

    // Interactive scene view (port of WebguiView in the legacy widget.ts).
    function renderScene(webgui, model, el) {
        el.classList.add('webgui-widget');
        // size the widget as requested from Python; the canvas container fills it

        const scene = new webgui.Scene();
        const container = document.createElement('div');
        container.style.width = model.get("width");
        container.style.height = model.get("height");
        el.appendChild(container);

        // defer to the next tick so the container is laid out before webgui
        // measures it (otherwise offsetHeight is read as ~0 -> tiny canvas)
        setTimeout(() => {
            resolveRenderData(model.get("value")).then((render_data) => {
                container.style.width = model.get("width");
                container.style.height = model.get("height");
                if(Object.keys(render_data).length === 0)
                    return;
                scene.init(container, render_data);
                scene.container.style.position = 'relative';
                scene.render();
                adoptDatGuiStyles(el);
                // dat.gui may inject its <style> a tick after construction; retry once.
                setTimeout(() => adoptDatGuiStyles(el), 100);
            });
        }, 0);

        // redraw: Python sets widget.value -> push new data into the scene
        model.on('change:value', () => {
            resolveRenderData(model.get("value")).then((render_data) => {
                scene.updateRenderData(render_data);
            });
        });
    }

    // Documentation view (port of WebguiDocuView): show a preview image and
    // only load the interactive scene + render data on click.
    function renderDocu(webgui, model, el) {
        const files = model.get("value");
        const container = document.createElement('div');
        container.className = 'webgui_container';
        container.style.width = '100%';
        container.innerHTML = `
            <img src="${files['preview']}" class="image">
            <div class="webgui_overlay webgui_tooltip">
                <span class="webgui_tooltiptext"> Click to load interactive WebGUI </span>
            </div>`;
        const div = document.createElement('div');
        div.appendChild(container);
        el.appendChild(div);

        container.addEventListener('click', () => {
            document.body.style.cursor = 'wait';
            fetch(files['render_data'])
                .then((resp) => resp.json())
                .then((render_data) => {
                    document.body.style.cursor = '';
                    const style = `width: ${el.clientWidth}px; height: ${el.clientHeight}px;`;
                    container.remove();
                    const pel = el.children[0];
                    pel.innerHTML = '';
                    pel.setAttribute('style', style);
                    const scene = new webgui.Scene();
                    scene.init(pel, render_data);
                    scene.render();
                    adoptDatGuiStyles(el);
                    setTimeout(() => adoptDatGuiStyles(el), 100);
                });
        });
    }

    async function render({ model, el }) {
        const webgui = await loadWebgui();
        const value = model.get("value");
        if (value && value.render_data !== undefined && value.preview !== undefined) {
            renderDocu(webgui, model, el);
        } else {
            renderScene(webgui, model, el);
        }
    }
    export default { render };
    """
    _css = """
    .webgui-widget {
        padding: 0px 2px;
    }
    .webgui_container { position: relative; }
    .webgui_container .image { width: 100%; display: block; }
    .webgui_overlay {
        position: absolute; top: 0; left: 0; width: 100%; height: 100%;
        cursor: pointer;
    }
    .webgui_tooltip .webgui_tooltiptext {
        visibility: hidden;
        background-color: rgba(85, 85, 85, 0.9); color: #fff;
        text-align: center; border-radius: 6px; padding: 5px 10px;
        position: absolute; top: 50%; left: 50%;
        transform: translate(-50%, -50%); white-space: nowrap;
    }
    .webgui_tooltip:hover .webgui_tooltiptext { visibility: visible; }

    .label3d {
      -moz-user-select: none;
      -ms-user-select:none;
      -o-user-select:none;
      -webkit-user-select: none;
      display:block;
      position: absolute;
      transform: translate(-50%, -50%);
      user-select:none;
      z-index: 1;
    }

    .tooltiptext {
      width: 120px;
      background-color: #555;
      color: #fff;
      text-align: center;
      border-radius: 6px;
      padding: 5px 0;
      position: absolute;
      z-index: 1;
      left: 50%;
      margin-left: -60px;
      opacity: 1;
      transition: opacity 0.3s;
    }


    """
    value = traitlets.Dict({}).tag(sync=True)
    width = traitlets.Unicode(_default_height).tag(sync=True)
    height = traitlets.Unicode(_default_height).tag(sync=True)


class WebGLScene(BaseWebGuiScene):
    def __init__(self, obj, args=[], kwargs={}):
        self.widget = WebGuiWidget()
        self.obj = obj
        self.args = args
        self.kwargs = kwargs
        self.encoding = "b64"
        width = kwargs.get("width", _default_width)
        height = kwargs.get("height", _default_height)
        if isinstance(width, (int,float)):
            width = str(width) + "px"
        if isinstance(height, (int,float)):
            height = str(height) + "px"
        self.widget.width = str(width)
        self.widget.height = str(height)
        self.widget.value = self.GetData()

    def __repr__(self):
        return "WebGLScene"

    def GenerateHTML(self, filename=None, template=None):
        self.encoding = 'b64'
        return GenerateHTML(self.GetData(), filename, template)

    def MakeScreenshot(self, filename, width=1200, height=600):
        self.encoding = 'b64'
        return _MakeScreenshot(self.GetData(), filename, width, height)

    def Draw(self, width = None, height = None):
        from IPython.display import display
        if width is not None:
            self.widget.width = str(width)
        if height is not None:
            self.widget.height = str(height)
        display(self.widget)

    def Redraw(self, *args, **kwargs):
        if args or kwargs:
            if 'show' not in kwargs:
                kwargs['show'] = False

            new_scene = Draw(*args, **kwargs)
            self.obj = new_scene.obj
            self.args = new_scene.args
            self.kwargs = new_scene.kwargs
        if self.widget is not None:
            self.widget.value = self.GetData()

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

        id_ = kwargs.get("id", None)
        js_code = kwargs.get("js_code", "")
        if id_ is not None:
            js_code += f"""
                if(window._webgui_scenes === undefined)
                    window._webgui_scenes = {{}};
                window._webgui_scenes['{id_}'] = scene;
            """
        if js_code:
            d["on_init"] = js_code

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

    def DownloadScreenshot(self, filename="image.jpg"):
        from IPython.display import Javascript, display

        if "id" not in self.kwargs:
            raise Exception(
                "To make a screenshot, please provide an id='some_unique_id' argument to Draw()"
            )
        format = filename.split(".")[-1].lower()
        display(
            Javascript(
                f"""
        {{
            var scene = window._webgui_scenes['{self.kwargs['id']}'];
            console.log("screenshot of scene", scene);
            const toimage = () => {{
                var link = document.createElement('a');
                link.href = scene.renderer.domElement.toDataURL('image/{format}');
                link.download = '{filename}';
                link.click();
                scene.event_handlers['afterrender'].pop(toimage);
            }};
            scene.on('afterrender', toimage);
            scene.render();
        }}
        """
            )
        )


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
    """Visualise a mesh or field in the webgui (Jupyter or standalone).

    The object type determines how it is rendered.  Netgen meshes
    (``netgen.meshing.Mesh``) are supported out of the box.  When
    *ngsolve* is imported, ``ngsolve.Mesh``, ``CoefficientFunction``
    and ``GridFunction`` are registered as well via
    :func:`register_draw_type`.

    Parameters
    ----------
    obj : mesh, function, or any registered draw type
        The object to visualise.  For an ``ngsolve.CoefficientFunction``
        pass the mesh or region as the first positional argument.
    show : bool
        Display the widget immediately in Jupyter (default ``True``).
    order : int
        Polynomial order used for visualisation (default 2).
    draw_vol : bool
        Draw volume elements (default ``True``).
    draw_surf : bool
        Draw surface elements (default ``True``).
    deformation : bool or GridFunction
        Apply deformation to the mesh (default ``False``).
    scale : float
        Scaling factor for the deformation (default 1.0).
    autoscale : bool
        Automatically determine color range (default ``True``).
    min, max : float
        Explicit color-range bounds (override *autoscale*).
    clipping : dict
        Clipping-plane specification, e.g.
        ``{"x": 1, "y": 0, "z": 0, "dist": 0}`` or with
        ``"vec"``/``"pnt"`` keys.
    vectors : bool or dict
        Show vector arrows.  Pass a dict with ``"grid_size"`` and/or
        ``"offset"`` for fine control.
    eval_function : str
        Custom JavaScript evaluation expression for the color map.
    objects : list
        Additional overlay objects to render.
    settings : dict
        GUI settings forwarded to the viewer.  Supported keys:

        - **Colormap** (*dict*) — ``autoscale`` (bool), ``ncolors`` (int),
          ``min`` / ``max`` (float).
        - **Clipping** (*dict*) — ``enable`` (bool), ``function`` (bool),
          ``x``, ``y``, ``z``, ``dist`` (float).
        - **Light** (*dict*) — ``ambient``, ``diffuse`` (float, 0–1),
          ``shininess`` (float, 0–100), ``specularity`` (float, 0–1).
        - **Vectors** (*dict*) — ``show`` (bool), ``grid_size`` (int),
          ``offset`` (float).
        - **Misc** (*dict*) — ``subdivision`` (int, 1–20),
          ``line_thickness`` (int, 1–20), ``fast_draw`` (bool).
        - **Complex** (*dict*) — ``phase`` (float), ``speed`` (float),
          ``animate`` (bool).
        - **Multidim** (*dict*) — ``t`` (float), ``speed`` (float),
          ``animate`` (bool).
        - **Objects** (*dict*) — visibility toggles by name, e.g.
          ``{"Edges": True, "Wireframe": False, "Surface": True}``.
        - **camera** (*dict*) — ``euler_angles`` ([x, y, z] in degrees),
          ``transformations`` (list of rotation/move dicts).
        - **axes_labels** (*list*) — axis labels, default
          ``["X", "Y", "Z"]``.

        Flat shortcuts are also accepted: ``autoscale``,
        ``colormap_min``, ``colormap_max``, ``colormap_ncolors``,
        ``subdivision``, ``line_thickness``, ``edges``, ``mesh``,
        ``elements``.
    euler_angles : tuple
        Initial camera orientation as ``[x, y, z]`` in degrees.
    center : list
        Override the centre of the scene.
    radius : float
        Override the bounding-sphere radius.
    interpolate_multidim : bool
        Interpolate between multi-dim components (default ``False``).
    animate : bool
        Animate multi-dim data (default ``False``).
    animate_complex : bool
        Animate the complex phase (default ``False``).
    colors : list
        Custom color table.
    fullscreen : bool
        Open in fullscreen mode (default ``False``).
    width, height : str
        Widget size, e.g. ``"100%"`` / ``"50vh"``.
    filename : str
        If given, export the scene to an HTML file.
    nodal_p1 : bool
        Use nodal P1 interpolation (default ``False``).
    js_code : str
        JavaScript snippet executed after scene initialisation.
    id : str
        Register the scene under this name in
        ``window._webgui_scenes`` for JS access.

    Returns
    -------
    WebGLScene
        The scene object.  Call ``scene.Redraw()`` to update after
        changes.

    See Also
    --------
    register_draw_type : Register a handler for additional types.
    """
    kwargs_with_defaults = _get_draw_default_args()
    kwargs_with_defaults.update(kwargs)

    scene = WebGLScene(obj, args, kwargs_with_defaults)
    if show and _IN_IPYTHON:
        from IPython.display import display
        display(scene.widget)
    if "filename" in kwargs_with_defaults:
        scene.GenerateHTML(filename=kwargs_with_defaults["filename"])
    return scene

_html_template = """
<!DOCTYPE html>
<html>
    <head>
        <title>NGSolve WebGUI</title>
        <meta name='viewport' content='width=device-width, user-scalable=no'/>
        <style>
            body{
                margin:0;
                overflow:hidden;
            }
            canvas{
                cursor:grab;
                cursor:-webkit-grab;
                cursor:-moz-grab;
            }
            canvas:active{
                cursor:grabbing;
                cursor:-webkit-grabbing;
                cursor:-moz-grabbing;
            }
        </style>
    </head>
    <body>
          <script src="https://cdn.jsdelivr.net/npm/webgui@0.2.39/dist/webgui.js"></script>
          </script>
          <script>
            {render}
            const scene = new webgui.Scene();
            scene.init(document.body, render_data, {preserveDrawingBuffer: false});
          </script>
    </body>
</html>
"""

def getScreenshotHTML():
    return _html_template.replace("preserveDrawingBuffer: false", "preserveDrawingBuffer: true")

def GenerateHTML(data, filename=None, template=None):
    if template is None:
        template = _html_template
    import json
    data = json.dumps(data)

    jscode = "var render_data = {}\n".format(data)
    html = _html_template.replace('{render}', jscode )

    if filename is not None:
        open(filename,'w').write( html )
    return html




async def _MakeScreenshot(data, png_file, width=1200, height=600):
    """Uses playwright to make a screenshot of the given html file."""
    # pylint: disable=import-outside-toplevel
    from playwright.async_api import async_playwright

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

    widget = WebGuiWidget()
    widget.value = {"render_data": data_file, "preview": preview_file}
    scene.widget = widget
    data = scene.GetData()
    json.dump(data, open(data_file_abs, "w"))
    asyncio.run(_MakeScreenshot(data, preview_file_abs, 1200, 600))
    scene.Redraw = lambda: None
    from IPython.display import display

    if scene.widget is not None:
        display(scene.widget)
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
