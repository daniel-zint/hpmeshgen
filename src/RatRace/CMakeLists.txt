add_library(RatRace STATIC 
    RatRace.cpp 
    RatRace.h
    DualHybridGraph.cpp 
    DualHybridGraph.h
)
find_package(CUDA  REQUIRED)

target_include_directories(RatRace PRIVATE 
	${SRC}
    ${OPENMESH_INCLUDE_DIR}
    ${CUDA_INCLUDE_DIRS} 
)

target_link_libraries(RatRace ${LIBS})
target_compile_definitions(RatRace PUBLIC -D_SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING)

add_library(RatRace::RatRace ALIAS RatRace)
