# Group Stuff
source_group("Config" FILES 
    Config.cpp 
    Config.h
)

#collect everything within the src folder
file(GLOB PROJECT_SRC_FILES 
    *.cpp 
    *.h
    *.cu
    */*.cpp 
    */*.h
    */*.cu
)

set(BASE_DIRECTORY "${CMAKE_SOURCE_DIR}")
configure_file(baseDirectory.h.in baseDirectory.h)

#add project executable
add_executable(${PROJECT_NAME}
    ${PROJECT_SRC_FILES}
)

target_include_directories(${PROJECT_NAME} PRIVATE 
    ${SRC}
    ${CMAKE_CURRENT_BINARY_DIR}
    ### externals ###
    ${OPENMESH_INCLUDE_DIR}
)

target_link_libraries(${PROJECT_NAME} PRIVATE 
    ${LIBS} 
    InputFileReader::InputFileReader
    Stopwatch::Stopwatch
	RatRace::RatRace
	Blossom5::Blossom5
	DMO::DMO
    Contour::Contour
	ScalarField::ScalarField
	GridEvaluation::GridEvaluation
)
target_compile_definitions(${PROJECT_NAME} PUBLIC -D_SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING)

if(MSVC)
    set_target_properties(${PROJECT_NAME} PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}")
endif()