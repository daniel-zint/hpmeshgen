set(BLOSSOM_INCLUDE ${SRC}/Blossom5)
file(GLOB BLOSSOM_FILES
	*.cpp
	*.h
	*/*.cpp
	*/*.h
)
add_library(Blossom5 STATIC 
    ${BLOSSOM_FILES}
)
target_compile_definitions(Blossom5 PUBLIC -D_SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING)

add_library(Blossom5::Blossom5 ALIAS Blossom5)