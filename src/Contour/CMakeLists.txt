add_library(Contour STATIC 
	ContourType.h
	Edge.h
	Point.h
    Shape.h
	Shape.cpp
	Types.h
)

find_package(CUDA  REQUIRED)
target_include_directories(Contour PUBLIC 
	${SRC}
	### externals ###
    ${OPENMESH_INCLUDE_DIR}
	${CUDA_INCLUDE_DIRS} 
)

target_link_libraries(Contour glog::glog ${OPENMESH_LIBRARIES})

target_compile_definitions(Contour PUBLIC -D_SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING)

add_library(Contour::Contour ALIAS Contour)
