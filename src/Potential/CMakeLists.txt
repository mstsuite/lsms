# Sources ----------------------------------------------------------------------
set(C_SRC
    calculateChargesPotential.cpp
    getvmt.cpp
    getXCName.cpp
    interpolatePotential.cpp
    libxcInterface.cpp
)

set(F_SRC
    getqm_mt.f
    newexchg.f
    alpha2_c.f
    newpot_c.f
    epcorr.f
    rs.f
    calculateASAro3.f
)

# Libraries --------------------------------------------------------------------
add_library(Potential STATIC ${C_SRC} ${F_SRC})

# Includes ---------------------------------------------------------------------
target_include_directories(Potential PRIVATE 
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/src
)