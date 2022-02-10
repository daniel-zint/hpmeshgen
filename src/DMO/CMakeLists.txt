add_library(DMO SHARED 
    DmoMesh.h
    DmoParams.h
    gpuErrchk.h
    kernel.cu
    Metrics.h
    Solver.h
    UniformGrid.h
    Vertex.h
    DMO_PUBLIC.h
)

target_include_directories(DMO PUBLIC 
	${SRC}
	### externals ###
    ${OPENMESH_INCLUDE_DIR}
    ${CUDA_INCLUDE_DIRS}
)

target_link_libraries(DMO 
	glog::glog 
	${OPENMESH_LIBRARIES}
	ScalarField::ScalarField
)

target_compile_options(DMO PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-std=c++17>)

set_target_properties(DMO PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

target_compile_definitions(DMO PUBLIC -D_SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING)
target_compile_definitions(DMO PRIVATE -DDMO_COMPILING)

add_library(DMO::DMO ALIAS DMO)