cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

if( CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR )
    message( FATAL_ERROR "Please select another Build Directory ! (and give it a clever name, like bin_Visual2012_64bits/)" )
endif()
if( CMAKE_SOURCE_DIR MATCHES " " )
	message( "Your Source Directory contains spaces. If you experience problems when compiling, this can be the cause." )
endif()
if( CMAKE_BINARY_DIR MATCHES " " )
	message( "Your Build Directory contains spaces. If you experience problems when compiling, this can be the cause." )
endif()

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

#project name
project(HPMeshGen2 LANGUAGES CUDA CXX)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED True)
set (CMAKE_POSITION_INDEPENDENT_CODE ON)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

if (CMAKE_CUDA_COMPILER)
    MESSAGE(STATUS "Generate with CUDA")
    enable_language(CUDA)
    #set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -gencode arch=compute_35,code=sm_35")
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-extended-lambda")
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-relaxed-constexpr")
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --use_fast_math")
	set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)
	set(CMAKE_CUDA_ARCHITECTURES 61)
else()
    MESSAGE(FATAL_ERROR "Cuda not found")
endif()

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH} )

########## required libraries ##########

#openmp
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
	set(LIBS ${LIBS} OpenMP::OpenMP_CXX)
endif()

#glog
find_package(glog REQUIRED)
set(LIBS ${LIBS} glog::glog)

#OpenMesh
find_package(OpenMesh REQUIRED)
set(LIBS ${LIBS} ${OPENMESH_LIBRARIES})
add_definitions(-D_USE_MATH_DEFINES)

#Eigen
#find_package(Eigen3 REQUIRED)
#set(LIBS ${LIBS} Eigen3::Eigen)

#nanoflann
#find_package(nanoflann CONFIG REQUIRED)
#set(LIBS ${LIBS} nanoflann::nanoflann)

#CImg
#find_package(CImg CONFIG)
#if(CImg_FOUND)
#	set(LIBS ${LIBS} CImg::CImg)
#endif()

if(NOT MSVC)
	#X11
	#MESSAGE(STATUS "Add X11 for non windows systems")
	#find_package(X11 REQUIRED)
	#set(LIBS ${LIBS} X11)
    MESSAGE(STATUS "Add stdc++fs for non windows systems")
	set(LIBS ${LIBS} stdc++fs)
endif()

if(MSVC)
    # enable parallel build in Visual Studio
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")
endif()


#if (MSVC)
#    add_compile_options(/W4 /WX)
#else()
#    add_compile_options(-Wall -Wextra -pedantic -Werror)
#endif()


########## executables ##########

set(SRC ${PROJECT_SOURCE_DIR}/src)

add_subdirectory(src)

########## unit tests ##########

option(BUILD_TESTS "Build unit tests" OFF)
if(${BUILD_TESTS})
    message(STATUS "Build unit tests")
    enable_testing()
    include(GoogleTest)
    add_subdirectory(tests)
endif()
