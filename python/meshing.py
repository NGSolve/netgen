from .libngpy._meshing import *

class _MeshsizeObject:
    @property
    def very_coarse(self):
        return MeshingParameters(curvaturesafety=1,
                                 segmentsperedge=0.3,
                                 grading=0.7,
                                 resthsurfcurvfac=0.25,
                                 resthchartdistfac=0.8,
                                 resthlinelengthfac=0.2,
                                 resthcloseedgefac=0.5,
                                 resthminedgelen=0.002,
                                 resthedgeanglefac=0.25,
                                 resthsurfmeshcurvfac=1.)
    @property
    def coarse(self):
        return MeshingParameters(curvaturesafety=1.5,
                                 segmentsperedge=0.5,
                                 grading=0.5,
                                 resthsurfcurvfac=0.5,
                                 resthchartdistfac=1,
                                 resthlinelengthfac=0.35,
                                 resthcloseedgefac=1,
                                 resthminedgelen=0.02,
                                 resthedgeanglefac=0.5,
                                 resthsurfmeshcurvfac=1.5)
    @property
    def moderate(self):
        return MeshingParameters(curvaturesafety=2,
                                 segmentsperedge=1,
                                 grading=0.3,
                                 resthsurfcurvfac=1.,
                                 resthchartdistfac=1.5,
                                 resthlinelengthfac=0.5,
                                 resthcloseedgefac=2,
                                 resthminedgelen=0.2,
                                 resthedgeanglefac=1,
                                 resthsurfmeshcurvfac=2.)
    @property
    def fine(self):
        return MeshingParameters(curvaturesafety=3,
                                 segmentsperedge=2,
                                 grading=0.2,
                                 resthsurfcurvfac=1.5,
                                 resthchartdistfac=2,
                                 resthlinelengthfac=1.5,
                                 resthcloseedgefac=3.5,
                                 resthminedgelen=1.,
                                 resthedgeanglefac=1.5,
                                 resthsurfmeshcurvfac=3.)

    @property
    def very_fine(self):
        return  MeshingParameters(curvaturesafety=5,
                                  segmentsperedge=3,
                                  grading=0.1,
                                  resthsurfcurvfac=3,
                                  resthchartdistfac=5,
                                  resthlinelengthfac=3,
                                  resthcloseedgefac=5,
                                  resthminedgelen=2.,
                                  resthedgeanglefac=3.,
                                  resthsurfmeshcurvfac=5.)

meshsize = _MeshsizeObject()
