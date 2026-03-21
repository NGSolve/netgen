Netgen mesh generator (ksugahar fork)

This fork adds APIs for externally imported mesh curving, enabling
`mesh.Curve(order)` with any CAD kernel (e.g. Coreform Cubit / ACIS)
without requiring OCC geometry or STEP files.

## Added APIs

### 1. `Element2D.SetGeomInfo(vertex_index, u, v)`

Set UV parametric coordinates on surface mesh elements.
Required to provide initial UV hints for `mesh.Curve(order)`.

```python
el = ngmesh.Elements2D()[0]
el.SetGeomInfo(0, u0, v0)  # vertex 0
el.SetGeomInfo(1, u1, v1)  # vertex 1
el.SetGeomInfo(2, u2, v2)  # vertex 2
```

PR: https://github.com/NGSolve/netgen/pull/232

### 2. `CallbackGeometry(project, normal, num_surfaces, tangent=None)`

Geometry backend with Python callbacks for surface projection.
Replaces OCC geometry for `mesh.Curve(order)`.

```python
from netgen.meshing import CallbackGeometry

def project(surfnr, x, y, z, u_hint, v_hint):
    # Project point onto surface, return (x, y, z, u, v)
    uv = cubit.surface(surfnr).u_v_from_position([x, y, z])
    pos = cubit.surface(surfnr).position_from_u_v(uv[0], uv[1])
    return pos[0], pos[1], pos[2], uv[0], uv[1]

def normal(surfnr, x, y, z):
    n = cubit.surface(surfnr).normal_at([x, y, z])
    return n[0], n[1], n[2]

geo = CallbackGeometry(project, normal, cubit.get_surface_count())
ngmesh.SetGeometry(geo)
mesh = Mesh(ngmesh)
mesh.Curve(3)  # High-order curving via ACIS — no OCC, no STEP
```

### 3. `Mesh.CalcSurfacesOfNode()`

Rebuild surface element lookup tables after manually adding elements.
Required for `HDivSurface` (BEM) on externally imported meshes.

### 4. `Mesh.RebuildSurfaceElementLists()`

Rebuild internal linked lists of surface elements per FaceDescriptor.
Required for correct BEM assembly after manual mesh construction.

## Install

```bash
# 1. Install official ngsolve
pip install ngsolve

# 2. Replace netgen-mesher with this fork (from GitHub Releases)
pip install <wheel-url> --force-reinstall
```

Wheels: https://github.com/ksugahar/netgen/releases

## Upstream

Based on Netgen 6.2.2602. Find the original at https://ngsolve.org
