add_definitions(-DNGINTERFACE_EXPORTS)
if(NOT WIN32)
    set(mesh_object_libs 
        $<TARGET_OBJECTS:la>
        $<TARGET_OBJECTS:gprim>
        $<TARGET_OBJECTS:gen>
        )
endif(NOT WIN32)

add_library(mesh ${NG_LIB_TYPE}
        adfront2.cpp adfront3.cpp bisect.cpp boundarylayer.cpp 
        clusters.cpp curvedelems.cpp delaunay.cpp delaunay2d.cpp	    
        geomsearch.cpp global.cpp hprefinement.cpp improve2.cpp		    
        improve2gen.cpp improve3.cpp localh.cpp meshclass.cpp		    
        meshfunc.cpp meshfunc2d.cpp meshing2.cpp meshing3.cpp		    
        meshtool.cpp meshtype.cpp msghandler.cpp netrule2.cpp		    
        netrule3.cpp parser2.cpp parser3.cpp prism2rls.cpp		    
        pyramid2rls.cpp pyramidrls.cpp quadrls.cpp refine.cpp		    
        ruler2.cpp ruler3.cpp secondorder.cpp smoothing2.5.cpp		    
        smoothing2.cpp smoothing3.cpp specials.cpp tetrarls.cpp		    
        topology.cpp triarls.cpp validate.cpp bcfunctions.cpp   
        parallelmesh.cpp  paralleltop.cpp  paralleltop.hpp basegeom.cpp     
        python_mesh.cpp hexarls.cpp
        ../../ng/onetcl.cpp
        ${mesh_object_libs}
        )

if(APPLE)
    set_target_properties( mesh PROPERTIES SUFFIX ".so")
endif(APPLE)

if(NOT WIN32)
    target_link_libraries( mesh ${ZLIB_LIBRARIES} ${MPI_CXX_LIBRARIES} ${PYTHON_LIBRARIES} ${METIS_LIBRARY})
    install( TARGETS mesh ${NG_INSTALL_DIR})
endif(NOT WIN32)

install(FILES
  adfront2.hpp adfront3.hpp basegeom.hpp bcfunctions.hpp bisect.hpp
  boundarylayer.hpp classifyhpel.hpp clusters.hpp curvedelems.hpp
  findip2.hpp findip.hpp geomsearch.hpp global.hpp hpref_hex.hpp
  hprefinement.hpp hpref_prism.hpp hpref_pyramid.hpp hpref_quad.hpp
  hpref_segm.hpp hpref_tet.hpp hpref_trig.hpp improve2.hpp improve3.hpp
  localh.hpp meshclass.hpp meshfunc.hpp meshing2.hpp meshing3.hpp
  meshing.hpp meshtool.hpp meshtype.hpp msghandler.hpp paralleltop.hpp
  ruler2.hpp ruler3.hpp specials.hpp topology.hpp validate.hpp
  DESTINATION ${NG_INSTALL_DIR_INCLUDE}/meshing COMPONENT netgen_devel
)
