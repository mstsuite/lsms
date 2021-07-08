# Sources ----------------------------------------------------------------------
set(C_SRC
    calculateMadelungMatrices.cpp
)

set(F_SRC
    cal_madelung_matrix.f
    getkncut.f
    getrscut.f
    getstruc.f
    interf.f
    interfsmr.f
    lattice.f
    madewd.f
    madewdj.f
    madsum.f
    ord3v.f
    pqintg_c.f
    lmfacts.f
    bessj.f
)

# Libraries --------------------------------------------------------------------
add_library(Madelung STATIC ${C_SRC} ${F_SRC})

# Includes ---------------------------------------------------------------------
target_include_directories(Madelung PRIVATE 
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/src
)