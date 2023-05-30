target_sources(nglib PRIVATE
        Partition_Inter2d.cxx Partition_Inter3d.cxx 
        Partition_Loop.cxx Partition_Loop2d.cxx Partition_Loop3d.cxx Partition_Spliter.cxx 
        occgenmesh.cpp occgeom.cpp occmeshsurf.cpp python_occ.cpp
        python_occ_basic.cpp python_occ_shapes.cpp
        occ_face.cpp occ_edge.cpp occ_vertex.cpp occ_utils.cpp
)

if(USE_GUI)
    target_sources(nggui PRIVATE vsocc.cpp occpkg.cpp)
endif(USE_GUI)

install(FILES
  occgeom.hpp occmeshsurf.hpp vsocc.hpp occ_utils.hpp
  occ_vertex.hpp occ_edge.hpp occ_face.hpp occ_solid.hpp
  DESTINATION ${NG_INSTALL_DIR_INCLUDE}/occ COMPONENT netgen_devel
)
