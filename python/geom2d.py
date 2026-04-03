from .libngpy._geom2d import SplineGeometry, Solid2d, CSG2d, Rectangle, Circle, EdgeInfo, PointInfo
from .meshing import meshsize
import math as math

unit_square = SplineGeometry()
_pnts = [ (0,0), (1,0), (1,1), (0,1) ]
_lines = [ (0,1,1,"bottom"), (1,2,2,"right"), (2,3,3,"top"), (3,0,4,"left") ]
_pnums = [unit_square.AppendPoint(*p) for p in _pnts]
for l1,l2,bc,bcname in _lines:
    unit_square.Append( ["line", _pnums[l1], _pnums[l2]], bc=bcname)


def MakeRectangle (geo, p1, p2, bc=None, bcs=None, **args):
    p1x, p1y = p1
    p2x, p2y = p2
    p1x,p2x = min(p1x,p2x), max(p1x, p2x)
    p1y,p2y = min(p1y,p2y), max(p1y, p2y)
    if not bcs: bcs=4*[bc]
    pts = [geo.AppendPoint(*p) for p in [(p1x,p1y), (p2x, p1y), (p2x, p2y), (p1x, p2y)]]
    for p1,p2,bc in [(0,1,bcs[0]), (1, 2, bcs[1]), (2, 3, bcs[2]), (3, 0, bcs[3])]:
        geo.Append( ["line", pts[p1], pts[p2]], bc=bc, **args)

def MakeCircle (geo, c, r, **args):
    cx,cy = c
    pts = [geo.AppendPoint(*p) for p in [(cx,cy-r), (cx+r,cy-r), (cx+r,cy), (cx+r,cy+r), \
                                         (cx,cy+r), (cx-r,cy+r), (cx-r,cy), (cx-r,cy-r)]]
    for p1,p2,p3 in [(0,1,2), (2,3,4), (4, 5, 6), (6, 7, 0)]:
        geo.Append( ["spline3", pts[p1], pts[p2], pts[p3]], **args)

    

def CreatePML(geo, pml_size, tol=1e-12):
    """Create a pml layer around the geometry. This function works only on convex geometries and
the highest existing domain number must be named by using the function geo.SetMaterial(domnr, name).
Points in the geometry are assumed to be the same if (pt1 - pt2).Norm() < tol.
Returned is a dict with information to create the pml layer:
  normals: A dict from the names of the linear pml domains to the normal vectors pointing inside the pml."""

    def Start(spline):
        if spline.rightdom == 0:
            return spline.StartPoint()
        return spline.EndPoint()
    def End(spline):
        if spline.rightdom == 0:
            return spline.EndPoint()
        return spline.StartPoint()

    splines = []
    for i in range(geo.GetNSplines()):
        splines.append(geo.GetSpline(i))
    border = []
    is_closed = False
    current_endpoint = None
    while not is_closed:
        for spline in splines:
            if spline.leftdom == 0 or spline.rightdom == 0:
                if current_endpoint is not None:
                    if (Start(spline)-current_endpoint).Norm() < tol:
                        border.append(spline)
                        current_endpoint = End(spline)
                        if (current_endpoint - startpoint).Norm() < tol:
                            is_closed = True
                        break
                else:
                    startpoint = Start(spline)
                    current_endpoint = End(spline)
                    border.append(spline)
                    break
        else:
            raise Exception("Couldn't find closed spline around domain")
    endpointindex_map = []
    for spline in border:
        pnt = End(spline)
        for i in range(geo.GetNPoints()):
            if (pnt - geo.GetPoint(i)).Norm() < tol:
                endpointindex_map.append(i)
                break
        else:
            raise Exception("Couldn't find endpoint of spline in geometry")
    start_ndoms = ndoms = geo.GetNDomains() + 1
    new_spline_domains = []
    normals = {}
    duplicate_cnt = 0

    for i, spline in enumerate(border):
        if i == 0:
            global_start = Start(spline) + pml_size * spline.GetNormal(0)
            global_start_pnt = current_start = geo.AppendPoint(global_start[0], global_start[1])
        next_spline = border[(i+1)%len(border)]
        new_end =  End(spline) + pml_size * spline.GetNormal(1)
        spline_name = geo.GetBCName(spline.bc)

        if "pml_" + spline_name in normals \
        and normals["pml_" + spline_name] != spline.GetNormal(0):
            duplicate_cnt += 1
            spline_name = spline_name + "_duplicate_" + str(duplicate_cnt)

        if (new_end - global_start).Norm() < tol:
            new_spline_domains.append(ndoms)
            geo.Append(["line", current_start, global_start_pnt], bc="outer_" + spline_name, leftdomain = ndoms)
            geo.Append(["line", global_start_pnt, endpointindex_map[i]], leftdomain=ndoms, rightdomain=start_ndoms)
            geo.SetMaterial(ndoms, "pml_" + spline_name)
            normals["pml_" + spline_name] = spline.GetNormal(0)
            ndoms += 1
            break
        end = geo.AppendPoint(new_end[0], new_end[1])
        new_spline_domains.append(ndoms)
        geo.Append(["line", current_start, end], bc="outer_" + spline_name, leftdomain = ndoms)
        geo.Append(["line", end, endpointindex_map[i]], leftdomain=ndoms, rightdomain=ndoms+1)
        geo.SetMaterial(ndoms, "pml_" + spline_name)
        normals["pml_" + spline_name] = spline.GetNormal(0)
        ndoms += 1
        new_start = Start(next_spline) + pml_size * next_spline.GetNormal(0)
        if (new_start - global_start).Norm() < tol:
            geo.Append(["line", end, global_start_pnt], bc="outer", leftdomain = ndoms)
            geo.Append(["line", global_start_pnt, endpointindex_map[i]], leftdomain=ndoms, rightdomain=start_ndoms)
            geo.SetMaterial(ndoms, "pml_corner")
            ndoms += 1
            break
        if (new_end - new_start).Norm() < tol:
            current_start = end
        else:
            current_start = geo.AppendPoint(new_start[0], new_start[1])
            geo.Append(["line", end, current_start], bc="outer", leftdomain = ndoms)
            geo.Append(["line", current_start, endpointindex_map[i]], leftdomain=ndoms, rightdomain=ndoms+1)
            geo.SetMaterial(ndoms, "pml_corner")
            ndoms += 1
    for spline, domnr in zip(border, new_spline_domains):
        if spline.leftdom == 0:
            spline.leftdom = domnr
        else:
            spline.rightdom = domnr
    return {"normals" : normals}

SplineGeometry.AddCircle = lambda geo, c, r, **args : MakeCircle(geo, c, r, **args)
SplineGeometry.AddRectangle = lambda geo, p1, p2, **args : MakeRectangle(geo, p1, p2, **args)
SplineGeometry.AddSegment = lambda *args, **kwargs : SplineGeometry.Append(*args, **kwargs)
SplineGeometry.AddPoint = lambda *args, **kwargs : SplineGeometry.AppendPoint(*args, **kwargs)
SplineGeometry.CreatePML = CreatePML

bc = lambda s : EdgeInfo(bc=s)
maxh = lambda h : EdgeInfo(maxh=h)
def cp(p_or_px, py_or_none = None):
    if py_or_none is None:
        return EdgeInfo(control_point=p)
    else:
        return EdgeInfo(control_point=(p_or_px,py_or_none))


def Ellipse(center, a, b, bc="ellipse", mat="ellipse"):
    """Creates ellipse centered at point center with principle axis a and b.

    Parameters
    ---------
    center : Vec2
      center of ellipse
    a : Vec2
      first principle axis, needs to be perpendicular to b
    b : Vec2
      second principle axis, needs to be perpendicular to a
    bc : string
      boundary name
    mat : string
      material name
    """
    if abs(a[0]*b[0] + a[1]*b[1]) > 1e-12:
        raise Exception("In Ellipse: principle axis a and b are not perpendicular")
    
    ellipse = Circle( center=(0,0), radius=1.0, mat=mat, bc=bc )
    
    alpha = math.pi/2-math.atan2(a[0],a[1])
    r_a = math.sqrt(a[0]**2+a[1]**2)
    r_b = math.sqrt(b[0]**2+b[1]**2)
    ellipse.Scale( (r_a,r_b) )
    ellipse.Rotate( alpha/math.pi*180, center=(0,0) )
    ellipse.Move( center )
    
    return ellipse
