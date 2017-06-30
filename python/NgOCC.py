from netgen.libngpy._NgOCC import *
from netgen.libngpy._meshing import MeshingParameters

def NgOCC_meshing_func (geom, **args):
    if "mp" in args:
        return GenerateMesh (geom, args["mp"])
    else:
        return GenerateMesh (geom, MeshingParameters (**args))

OCCGeometry.GenerateMesh = NgOCC_meshing_func
