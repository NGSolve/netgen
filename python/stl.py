from netgen.libngpy._stl import *
from netgen.libngpy._meshing import MeshingParameters


def stl_meshing_func (geom, **args):
    if "mp" in args:
        return GenerateMesh (geom, args["mp"])
    else:
        return GenerateMesh (geom, MeshingParameters (**args))
#     return GenerateMesh (geom, MeshingParameters (**args))

STLGeometry.GenerateMesh = stl_meshing_func
