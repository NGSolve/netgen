from .libngpy._meshing import *

class _MeshsizeObject:
    pass

meshsize = _MeshsizeObject()

meshsize.very_coarse = MeshingParameters(curvaturesafety=1,
                                         segmentsperedge=0.3,
                                         grading=0.7)
meshsize.coarse = MeshingParameters(curvaturesafety=1.5,
                                    segmentsperedge=0.5,
                                    grading=0.5)
meshsize.moderate = MeshingParameters(curvaturesafety=2,
                                      segmentsperedge=1,
                                      grading=0.3)
meshsize.fine = MeshingParameters(curvaturesafety=3,
                                  segmentsperedge=2,
                                  grading=0.2)
meshsize.very_fine = MeshingParameters(curvaturesafety=5,
                                       segmentsperedge=3,
                                       grading=0.1)
