file(GLOB GRIDEVALUATION_SRC_FILES 
    *.cpp 
    *.h
)
find_package(CUDA REQUIRED) 
add_library(GridEvaluation STATIC
    ${GRIDEVALUATION_SRC_FILES}
)

target_link_libraries(GridEvaluation PRIVATE glog::glog)

target_include_directories(GridEvaluation PRIVATE 
    ### externals ###
    ${OPENMESH_INCLUDE_DIR}
    ${CUDA_INCLUDE_DIRS} 
)

target_include_directories(GridEvaluation PUBLIC
    ..
)

target_compile_definitions(GridEvaluation PUBLIC -D_SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING)

add_library(GridEvaluation::GridEvaluation ALIAS GridEvaluation)
