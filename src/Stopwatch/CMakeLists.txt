add_library(Stopwatch STATIC
    Stopwatch.cpp 
    Stopwatch.h
)

target_compile_definitions(Stopwatch PUBLIC -D_SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING)

add_library(Stopwatch::Stopwatch ALIAS Stopwatch)