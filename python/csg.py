from .libngpy._csg import *
from .libngpy._meshing import Pnt, Vec, Trafo
from .meshing import meshsize

try:
    from . import csgvis
    from .csgvis import MouseMove
    CSGeometry.VS = csgvis.VS
    SetBackGroundColor = csgvis.SetBackGroundColor
    del csgvis

    def VS (obj):
        return obj.VS()
except:
    pass

unit_cube = CSGeometry()
p1 = Plane(Pnt(0,0,0),Vec(-1,0,0)).bc("back")
p2 = Plane(Pnt(1,1,1),Vec(1,0,0)).bc("front")
p3 = Plane(Pnt(0,0,0),Vec(0,-1,0)).bc("left")
p4 = Plane(Pnt(1,1,1),Vec(0,1,0)).bc("right")
p5 = Plane(Pnt(0,0,0),Vec(0,0,-1)).bc("bottom")
p6 = Plane(Pnt(1,1,1),Vec(0,0,1)).bc("top")
unit_cube.Add (p1*p2*p3*p4*p5*p6, col=(0,0,1))

