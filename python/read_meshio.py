from netgen.meshing import *

def ReadViaMeshIO(filename):
    import meshio
    import numpy as np
    
    # print ("reading via meshio:", filename)

    m = meshio.read(filename)
    pts = m.points
    
    mesh = Mesh(dim=pts.shape[1])
    # mesh.AddPoints ( np.asarray(pts, dtype=np.float64) )  # needed for correct little/big endian
    mesh.AddPoints ( pts )  

    fd = mesh.Add (FaceDescriptor(bc=1))
    for cb in m.cells:
        # mesh.AddElements(dim=cb.dim, index=1, data=np.asarray(cb.data, dtype=np.int32), base=0)
        mesh.AddElements(dim=cb.dim, index=1, data=cb.data, base=0)
            
    return mesh

