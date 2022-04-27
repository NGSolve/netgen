target_sources(nglib PRIVATE
        nginterface.cpp nginterface_v2.cpp 
        read_fnf_mesh.cpp readtetmesh.cpp readuser.cpp writeabaqus.cpp writediffpack.cpp 
        writedolfin.cpp writeelmer.cpp writefeap.cpp writefluent.cpp writegmsh.cpp writejcm.cpp 
        writepermas.cpp writetecplot.cpp writetet.cpp writetochnog.cpp writeuser.cpp
        wuchemnitz.cpp writegmsh2.cpp writeOpenFOAM15x.cpp rw_cgns.cpp
)

install(FILES
  writeuser.hpp
  DESTINATION ${NG_INSTALL_DIR_INCLUDE}/interface COMPONENT netgen_devel
)
