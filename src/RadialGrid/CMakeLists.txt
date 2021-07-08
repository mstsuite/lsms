set(SRC
    RadialGrid.cpp
)

# Libraries --------------------------------------------------------------------
add_library(RadialGrid STATIC ${SRC})

target_include_directories(RadialGrid PRIVATE 
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/src
)