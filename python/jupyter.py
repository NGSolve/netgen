import json
import os

import requests
from IPython.display import HTML, DisplayHandle, Javascript, display

_id = 0  # counter for unique id for each generated div element

# VSCode has a strange meaning of "vh" in jupyter notebooks, use pixels instead
_default_height = "500px" if "VSCODE_PID" in os.environ else "50vh"

_webgui_loaded = False


def load_webgui(url: str = "", embed_code: bool = False):
    global _webgui_loaded
    if _webgui_loaded:
        return
    _webgui_loaded = True

    url = url or "https://cdn.jsdelivr.net/npm/webgui@0.2.38/dist/webgui.js"

    if embed_code:
        if url.startswith("http"):
            webgui_code = requests.get(url).text
        else:
            with open(url, "r") as f:
                webgui_code = f.read()
        display(HTML(f"<script>{webgui_code}</script>"))
    else:
        display(HTML(f"<script src={url}> </script>"))


def render(
    data, width="100%", height=_default_height, handle: DisplayHandle | None = None
):
    load_webgui()
    if handle is None:
        global _id
        _id += 1
        id = _id
    else:
        id = handle.display_id.split("_")[-1]

    el_id = f"_webgui_root_{id}"

    js_code = _render_js_template.replace("{{data}}", json.dumps(data))
    js_code = js_code.replace("{{id}}", str(id))
    js_code = js_code.replace("{{el_id}}", el_id)
    if handle:
        # this is a redraw call
        handle.update(Javascript(js_code))
    else:
        # this is a first draw call -> create the root div element for the webgui
        display(
            HTML(
                f'<div class="webgui-widget" id={el_id} style="position: relative; width: {width}; height:{height};" />'
            )
        )
        return display(Javascript(js_code), display_id=f"webgui_render_{id}")


_render_js_template = """
{
    if(window._webgui_scenes === undefined) {
        window._webgui_scenes = {};
    }
    const waitForWebgui = (callback, counter) => {
        counter = counter || 0;
        if (window.webgui) {
            callback();
        } else {
            if (counter < 20) {
                setTimeout(() => {
                    waitForWebgui(callback, counter + 1);
                }, 200);
            }
            else {
                console.log("Error: Webgui not loaded");
            }
        }
    };

    const draw = () => {
        const data = JSON.parse(`{{data}}`);
        let scene = window._webgui_scenes[{{id}}];
        // console.log("have data, scene = ", scene);
        if(scene === undefined) {
            console.log("init scene");
            const root = document.getElementById("{{el_id}}");
            console.log("root element", root);
            scene = new webgui.Scene();
            scene.init(root, data);
            window._webgui_scenes[{{id}}] = scene;
        }
        else {
            scene.updateRenderData(data);
        }
        // console.log("scene", scene);
    }
    waitForWebgui(draw);
}
"""
