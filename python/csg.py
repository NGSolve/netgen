from netgen.libngpy._csg import *
from netgen.libngpy._meshing import MeshingParameters
from netgen.libngpy._meshing import Pnt
from netgen.libngpy._meshing import Vec


try:
    import libngpy.csgvis as csgvis
    from libngpy.csgvis import MouseMove
    CSGeometry.VS = csgvis.VS
    SetBackGroundColor = csgvis.SetBackGroundColor
    del csgvis

    def VS (obj):
        return obj.VS()

except:
    pass


def csg_meshing_func (geom, **args):
    if "mp" in args:
        return GenerateMesh (geom, args["mp"])
    else:
        return GenerateMesh (geom, MeshingParameters (**args))
#     return GenerateMesh (geom, MeshingParameters (**args))

CSGeometry.GenerateMesh = csg_meshing_func


unit_cube = CSGeometry()
unit_cube.Add (OrthoBrick(Pnt(0,0,0), Pnt(1,1,1)))

