""" Netgen OCC documentation 

This module provides the Netgen OCCGeometry class, as well as Python wrappers 
for OpenCascade functionality. This allows the construction from scratch, as well as
editing geometries imported via STEP format.

Most of the classes are directly py-wrapped from the OCCT classes. 
For more detailed documentation consult the OCCT docs, a good starting point is
https://dev.opencascade.org/doc/refman/html/class_b_rep_builder_a_p_i___make_shape.html
"""

from .libngpy._NgOCC import *
from .meshing import meshsize


gp_Ax3 = Axes
gp_Ax1 = Axis

Translation = gp_Trsf.Translation
Rotation = gp_Trsf.Rotation
Mirror = gp_Trsf.Mirror
