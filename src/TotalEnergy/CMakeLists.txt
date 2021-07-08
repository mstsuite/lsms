# Sources ----------------------------------------------------------------------
set(SRC
    calculateTotalEnergy.cpp
    janake_c.f
    zeropt_c.f
)

# Libraries --------------------------------------------------------------------
add_library(TotalEnergy STATIC ${SRC})

# Includes ---------------------------------------------------------------------
target_include_directories(TotalEnergy PRIVATE 
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/src
)