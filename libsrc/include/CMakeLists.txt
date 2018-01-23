install(FILES nginterface.h nginterface_v2.hpp mydefs.hpp DESTINATION ${NG_INSTALL_DIR_INCLUDE} COMPONENT netgen_devel)

install(FILES
  acisgeom.hpp csg.hpp geometry2d.hpp gprim.hpp incopengl.hpp
  inctcl.hpp incvis.hpp linalg.hpp meshing.hpp myadt.hpp mydefs.hpp
  mystdlib.h nginterface_v2_impl.hpp occgeom.hpp ngsimd.hpp
  opti.hpp parallel.hpp parallelinterface.hpp stlgeom.hpp visual.hpp
  DESTINATION ${NG_INSTALL_DIR_INCLUDE}/include COMPONENT netgen_devel
  )

