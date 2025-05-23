""" Netgen OCC documentation 

This module provides the Netgen OCCGeometry class, as well as Python wrappers 
for OpenCascade functionality. This allows the construction from scratch, as well as
editing geometries imported via STEP format.

Most of the classes are directly py-wrapped from the OCCT classes. 
For more detailed documentation consult the OCCT docs, a good starting point is
https://dev.opencascade.org/doc/refman/html/class_b_rep_builder_a_p_i___make_shape.html
"""

from .config import USE_OCC
if not USE_OCC:
    raise ImportError("Netgen was not built with Opencascade support")

from .libngpy._NgOCC import *
from .meshing import meshsize, IdentificationType


gp_Ax3 = Axes
gp_Ax1 = Axis

Translation = gp_Trsf.Translation
Rotation = gp_Trsf.Rotation
Mirror = gp_Trsf.Mirror


def Rectangle(l,w): return WorkPlane().Rectangle(l,w)
def MoveTo(x,y): return WorkPlane().MoveTo(x,y)
def LineTo(x,y): return WorkPlane().LineTo(x,y)
def Line(l): return WorkPlane().Line(l)    


unit_square_shape = WorkPlane().Line(1, name="bottom").Rotate(90) \
  .Line(1, name="right").Rotate(90) \
  .Line(1, name="top").Rotate(90) \
  .Line(1, name="left").Rotate(90).Face()

unit_square = OCCGeometry(unit_square_shape, dim=2)



uc_shape = Box((0,0,0),(1,1,1))
uc_shape.faces.Max(X).name = "front"
uc_shape.faces.Min(X).name = "back"
uc_shape.faces.Max(Y).name = "right"
uc_shape.faces.Min(Y).name = "left"
uc_shape.faces.Max(Z).name = "top"
uc_shape.faces.Min(Z).name = "bottom"
    
unit_cube = OCCGeometry(uc_shape)

