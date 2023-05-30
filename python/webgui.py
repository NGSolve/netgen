import math
import numpy as np
from time import time
import os

from webgui_jupyter_widgets import BaseWebGuiScene, encodeData, WebGuiDocuWidget
import webgui_jupyter_widgets.widget as wg

class WebGLScene(BaseWebGuiScene):
    def __init__(self, mesh, clipping, on_init):
        from IPython.display import display, Javascript
        import threading
        self.mesh = mesh
        self.clipping = clipping
        self.on_init = on_init

    def GetData(self, set_minmax=True):
        import json
        #  d = BuildRenderData(self.mesh, self.cf, self.order, draw_surf=self.draw_surf, draw_vol=self.draw_vol, deformation=self.deformation, region=self.region)
        d = self.mesh._webgui_data()
        bp = d['Bezier_trig_points']
        for i in range(len(bp)):
            bp[i] = encodeData(np.array(bp[i], dtype=np.float32))

        ep = d['edges']
        for i in range(len(ep)):
            ep[i] = encodeData(np.array(ep[i], dtype=np.float32))

        if self.clipping is not None:
            d['clipping'] = True
            if isinstance(self.clipping, dict):
                allowed_args = ("x", "y", "z", "dist", "function", "pnt", "vec")
                if "vec" in self.clipping:
                    vec = self.clipping["vec"]
                    self.clipping["x"] = vec[0]
                    self.clipping["y"] = vec[1]
                    self.clipping["z"] = vec[2]
                if "pnt" in self.clipping:
                    d['mesh_center'] = list(self.clipping["pnt"])
                for name, val in self.clipping.items():
                    if not (name in allowed_args):
                        raise Exception('Only {} allowed as arguments for clipping!'.format(", ".join(allowed_args)))
                    d['clipping_' + name] = val

        if self.on_init:
            d['on_init'] = self.on_init

        return d


bezier_trig_trafos = { }  # cache trafos for different orders

def BuildRenderData(mesh, func, order=2, draw_surf=True, draw_vol=True, deformation=None, region=True):
    d = {}
    d['ngsolve_version'] = "Netgen"
    d['mesh_dim'] = 3 # mesh.dim TODO

    d['order2d'] = 1
    d['order3d'] = 0

    d['draw_vol'] = False
    d['draw_surf'] = True
    d['funcdim'] = 1

    func2 = None
    if func and func.is_complex:
        d['is_complex'] = True
        func1 = func[0].real
        func2 = ngs.CoefficientFunction( (func[0].imag, 0.0) )
    elif func and func.dim>1:
        func1 = func[0]
        func2 = ngs.CoefficientFunction( tuple(func[i] if i<func.dim else 0.0 for i in range(1,3)) ) # max 3-dimensional functions
        d['funcdim'] = func.dim
    elif func:
        func1 = func
        d['funcdim'] = 1
    else:
        # no function at all -> we are just drawing a mesh, eval mesh element index instead
        mats = mesh.GetMaterials()
        bnds = mesh.GetBoundaries()
        nmats = len(mesh.GetMaterials())
        nbnds = len(mesh.GetBoundaries())
        n = max(nmats, nbnds)
        func1 = ngs.CoefficientFunction(list(range(n)))
        n_regions = [0, 0, nmats, nbnds]
        d['mesh_regions_2d'] = n_regions[mesh.dim]
        d['mesh_regions_3d'] = nmats if mesh.dim==3 else 0
        d['funcdim'] = 0
    func1 = ngs.CoefficientFunction( (ngs.x, ngs.y, ngs.z, func1 ) )
    func0 = ngs.CoefficientFunction( (ngs.x, ngs.y, ngs.z, 0.0 ) )
    if deformation is not None:
        func1 += ngs.CoefficientFunction((deformation, 0.0))
        func0 += ngs.CoefficientFunction((deformation, 0.0))

    d['show_wireframe'] = False
    d['show_mesh'] = True
    if order2d>0:
        og = order2d
        d['show_wireframe'] = True
        d['show_mesh'] = True
        timer2.Start()

        timer3Bvals.Start()

        # transform point-values to Bernsteinbasis
        def Binomial(n,i): return math.factorial(n) / math.factorial(i) / math.factorial(n-i)
        def Bernstein(x, i, n): return Binomial(n,i) * x**i*(1-x)**(n-i)
        Bvals = ngs.Matrix(og+1,og+1)
        for i in range(og+1):
            for j in range(og+1):
                Bvals[i,j] = Bernstein(i/og, j, og)
        iBvals = Bvals.I
        timer3Bvals.Stop()        
        # print (Bvals)
        # print (iBvals)

                
        Bezier_points = [] 

        ipts = [(i/og,0) for i in range(og+1)] + [(0, i/og) for i in range(og+1)] + [(i/og,1.0-i/og) for i in range(og+1)]
        ir_trig = ngs.IntegrationRule(ipts, [0,]*len(ipts))
        ipts = [(i/og,0) for i in range(og+1)] + [(0, i/og) for i in range(og+1)] + [(i/og,1.0) for i in range(og+1)] + [(1.0, i/og) for i in range(og+1)]
        ir_quad = ngs.IntegrationRule(ipts, [0,]*len(ipts))

        vb = [ngs.VOL, ngs.BND][mesh.dim-2]
        if region and region.VB() == vb:
            vb = region
        cf = func1 if draw_surf else func0
        timer2map.Start()
        pts = mesh.MapToAllElements({ngs.ET.TRIG: ir_trig, ngs.ET.QUAD: ir_quad}, vb)
        timer2map.Stop()
        pmat = cf(pts)

        timermult.Start()
        pmat = pmat.reshape(-1, og+1, 4)
        if False:
            BezierPnts = np.tensordot(iBvals.NumPy(), pmat, axes=(1,1))
        else:
            BezierPnts = np.zeros( (og+1, pmat.shape[0], 4) )
            for i in range(4):
                ngsmat = ngs.Matrix(pmat[:,:,i].transpose()) 
                BezierPnts[:,:,i] = iBvals * ngsmat
        timermult.Stop()
        
        timer2list.Start()        
        for i in range(og+1):
            Bezier_points.append(encodeData(BezierPnts[i], dtype=np.float32))
        timer2list.Stop()        

        d['Bezier_points'] = Bezier_points

        ipts = [(i/og,0) for i in range(og+1)]
        ir_seg = ngs.IntegrationRule(ipts, [0,]*len(ipts))
        vb = [ngs.VOL, ngs.BND, ngs.BBND][mesh.dim-1]
        if region and region.VB() == vb:
            vb = region
        pts = mesh.MapToAllElements(ir_seg, vb)
        pmat = func0(pts)
        pmat = pmat.reshape(-1, og+1, 4)
        edge_data = np.tensordot(iBvals.NumPy(), pmat, axes=(1,1))
        edges = []
        for i in range(og+1):
            edges.append(encodeData(edge_data[i], dtype=np.float32))
        d['edges'] = edges

        ndtrig = int((og+1)*(og+2)/2)
        
        if og in bezier_trig_trafos.keys():
            iBvals_trig = bezier_trig_trafos[og]
        else:
            def BernsteinTrig(x, y, i, j, n):
                return math.factorial(n)/math.factorial(i)/math.factorial(j)/math.factorial(n-i-j) \
                  * x**i*y**j*(1-x-y)**(n-i-j)
            Bvals = ngs.Matrix(ndtrig, ndtrig)
            ii = 0
            for ix in range(og+1):
                for iy in range(og+1-ix):
                    jj = 0
                    for jx in range(og+1):
                        for jy in range(og+1-jx):
                            Bvals[ii,jj] = BernsteinTrig(ix/og, iy/og, jx, jy, og)
                            jj += 1
                    ii += 1
            iBvals_trig = Bvals.I
            bezier_trig_trafos[og] = iBvals_trig

            
        # Bezier_points = [ [] for i in range(ndtrig) ]
        Bezier_points = []
        
        ipts = [(i/og,j/og) for j in range(og+1) for i in range(og+1-j)]
        ir_trig = ngs.IntegrationRule(ipts, [0,]*len(ipts))
        ipts = ([(i/og,j/og) for j in range(og+1) for i in range(og+1-j)] + 
                [(1-i/og,1-j/og) for j in range(og+1) for i in range(og+1-j)])
        ir_quad = ngs.IntegrationRule(ipts, [0,]*len(ipts))
        
        vb = [ngs.VOL, ngs.BND][mesh.dim-2]
        if region and region.VB() == vb:
            vb = region
        pts = mesh.MapToAllElements({ngs.ET.TRIG: ir_trig, ngs.ET.QUAD: ir_quad}, vb)

        pmat = ngs.CoefficientFunction( func1 if draw_surf else func0 ) (pts)

        
        funcmin = np.min(pmat[:,3])
        funcmax = np.max(pmat[:,3])
        pmin = np.min(pmat[:,0:3], axis=0)
        pmax = np.max(pmat[:,0:3], axis=0)

        mesh_center = 0.5*(pmin+pmax)
        mesh_radius = np.linalg.norm(pmax-pmin)/2
        
        pmat = pmat.reshape(-1, len(ir_trig), 4)

        BezierPnts = np.tensordot(iBvals_trig.NumPy(), pmat, axes=(1,1))

        for i in range(ndtrig):
            Bezier_points.append(encodeData(BezierPnts[i], dtype=np.float32))


        d['Bezier_trig_points'] = Bezier_points    
        d['mesh_center'] = list(mesh_center)
        d['mesh_radius'] = mesh_radius



    if func:
        d['funcmin'] = funcmin
        d['funcmax'] = funcmax
    return d

def Draw(shape, clipping=None, js_code=None, filename=""):
    # todo: also handle occ geometry, list of shapes, etc.

    scene = WebGLScene(shape, clipping=clipping, on_init=js_code)

    if wg._IN_IPYTHON:
        if wg._IN_GOOGLE_COLAB:
            from IPython.display import display, HTML
            html = scene.GenerateHTML()
            display(HTML(html))
        else:
            # render scene using widgets.DOMWidget
            scene.Draw()
            return scene
    else:
        if filename:
            scene.GenerateHTML(filename=filename)
        return scene

