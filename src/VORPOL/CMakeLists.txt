# Sources ----------------------------------------------------------------------
set(SRC
    setup_vorpol.f
    setup_boundary.f
    filter_edge.f
    chkedge.f
    invm3.f
    chkbnd.f
    celbnd.f
    polyhedron.f
    sort.f
    sortidx.f
    rcritpts.f
    volvor.f
    stepyll.f
    caltnode.f
    intpl0.f
    calsig.f
    intphi.f
    inter_dip.f
    inter.f
    interstitial.f
    sigma.f
    setupVorpol.cpp
)

# Libraries --------------------------------------------------------------------
add_library(VORPOL STATIC ${SRC})

# Includes ---------------------------------------------------------------------
target_include_directories(VORPOL PRIVATE 
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/src
)