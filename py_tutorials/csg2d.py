from ngsolve import *

from random import random, seed
ngsglobals.msg_level = 0

import netgen
from pyngcore import *
from netgen.geom2d import *

seed(4)

def GenerateMesh():

    g = CSG2d()
    g1 = CSG2d()
    outer = Rectangle(0, 1, 0, 1,"outer","outer")
    inner = Solid2d()

    for i in range(30):
        cx = random()
        cy = random()
        r  = 0.03+0.05*random()
        print("Add Circle", i, cx, cy, r, flush = True)
        circle = Circle(cx, cy, r, "circle"+str(i), "circle"+str(i))
        g1.Add(circle)
        inner += circle
        outer -= circle


    g.Add(inner)
    g.Add(outer)
    geo = g.GenerateSplineGeometry()
    Draw(geo)

    # draw this geometry for checking ff the final mesh/geometry is correct
    # g1.Add(outer)
    # geo1 = g1.GenerateSplineGeometry()
    # Draw(geo1)

    print('generate mesh')
    m = geo.GenerateMesh(maxh=0.1)
    mesh = Mesh(m)
    mesh.Curve(3)
    Draw(mesh)

    return mesh

from ngsolve import Draw
with PajeTrace():
    mesh = GenerateMesh()
