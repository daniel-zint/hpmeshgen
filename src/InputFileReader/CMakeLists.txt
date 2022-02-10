add_library(InputFileReader STATIC
    InputFileReader.cpp 
    InputFileReader.h
)

target_link_libraries(InputFileReader PRIVATE glog::glog)

target_compile_definitions(InputFileReader PUBLIC -D_SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING)

add_library(InputFileReader::InputFileReader ALIAS InputFileReader)