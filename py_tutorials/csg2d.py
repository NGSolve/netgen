from random import random, seed
from ngsolve import Draw, Mesh

import netgen
from pyngcore import *
from netgen.geom2d import *

seed(4)

g = CSG2d()
outer = Rectangle((0, 0), (1, 1), "outer","outer")
inner = Solid2d()

for i in range(30):
    cx = random()
    cy = random()
    r  = 0.03+0.05*random()
    print("Add Circle", i, cx, cy, r, flush = True)
    circle = Circle((cx, cy), r, "circle"+str(i), "circle"+str(i))
    inner += circle
    outer -= circle


g.Add(inner)
g.Add(outer)
geo = g.GenerateSplineGeometry()

m = geo.GenerateMesh(maxh=0.1)

try:
    from ngsolve import Draw, Mesh
    Draw(geo)
    mesh = Mesh(m)
    mesh.Curve(3)
    Draw(mesh)
except:
    pass
