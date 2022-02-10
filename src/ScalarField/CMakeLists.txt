add_library(ScalarField STATIC 
    ScalarField.cpp 
    ScalarField.h
)

target_include_directories(ScalarField PRIVATE 
    ${SRC}
	${OPENMESH_INCLUDE_DIR}
)

set(SCALARFIELD_LIBS glog::glog ${OPENMESH_LIBRARIES})
#if(CImg_FOUND)
#    set(SCALARFIELD_LIBS ${SCALARFIELD_LIBS} CImg::CImg)
#    target_compile_definitions(ScalarField PRIVATE SCALARFIELD_USING_CIMG)
#endif()

#if(NOT MSVC)
#	set(SCALARFIELD_LIBS ${SCALARFIELD_LIBS} X11)
#endif()

target_link_libraries(ScalarField ${SCALARFIELD_LIBS})

target_compile_definitions(ScalarField PUBLIC -D_SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING)

add_library(ScalarField::ScalarField ALIAS ScalarField)
