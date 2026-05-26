from netgen.occ import *
from netgen.meshing import NetgenGeometry, Mesh as NGMesh
import numpy as np


class UnitSphereGeometry(NetgenGeometry):
    def midpoint(self, newp, p1, p2, secpoint=0.5):
        p1 = np.array([p1[0], p1[1], p1[2]])
        p2 = np.array([p2[0], p2[1], p2[2]])
        p = p1 + secpoint * (p2 - p1)
        self.project(newp, p)

    def project(self, newp, p):
        pt = np.array([p[0], p[1], p[2]])
        pt /= np.linalg.norm(pt)
        newp[0] = pt[0]
        newp[1] = pt[1]
        newp[2] = pt[2]

    def PointBetweenEdge(self, p1, p2, secpoint, surfi1, surfi2, ep1, ep2, newp, newgi):
        self.midpoint(newp, p1, p2, secpoint)

    def PointBetween(self, p1, p2, secpoint, surfi, gi1, gi2, newp, newgi):
        self.midpoint(newp, p1, p2, secpoint)

    def ProjectPointGI(self, surfind, p, gi):
        self.project(p, p)
        return True


m = OCCGeometry(Sphere(Pnt(0, 0, 0), 1)).GenerateMesh()
my_geo = UnitSphereGeometry()
m.SetGeometry(my_geo)
m.Refine()
m.Curve(3)

from ngsolve import *

mesh = Mesh(m)
Draw(mesh)
