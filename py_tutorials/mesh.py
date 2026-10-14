from netgen.meshing import *
from netgen.csg import *

geo = CSGeometry("shaft.geo")

param = MeshingParameters()
param.maxh = 10
print (param)

m1 = GenerateMesh (geo, param)
m1.SecondOrder()

m1.Export('shaft.mesh', 'Neutral Format')



Save (m1, "mesh.vol", geo)

