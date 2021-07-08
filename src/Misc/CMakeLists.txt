# Sources ----------------------------------------------------------------------
set(SRC
    zeroout.f
    fstop.f
    newint.f
    fit.f
    interp.f
    trltog.f
    spin_trafo.f
    fitpot.f
    rwave.f
    initwave.f
    dfv.f
    dfv_new.f
    rzextr.f
    zsphbesjh.f
    zsphbesj.f
    cinterp.f
    cgaunt_c.f
    gaunt.f
    quadrature.cpp
    ifacts_c.f
    clock_time.c
    cmtruni.f
    constraint.f
    v_plus_minus.f
    congauss_c.f
    u_sigma_u.f
    fnpi.f
    wrtmtx.f
    mbeqa.f
    cnewint.f
    ylag.f
    ricbes.f
    newder.f
    clebsch.f
    rsimp.f
    matrot1.f
    matr.f
    rotmat.f
    bulirsch_stoer.f
    mod_midpoint.f
    bulirschStoerIntegrator.cpp
    calculateGauntCoeficients.cpp
    associatedLegendreFunction.cpp
    readLastLine.cpp
    stop_with_backtrace.cpp
)

if(ESSL_WORKAROUND)
    list(APPEND SRC essl_workaround.f)
endif()

# Target -----------------------------------------------------------------------
add_library(Misc STATIC ${SRC})

# Includes ---------------------------------------------------------------------
target_include_directories(Misc PRIVATE 
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/src
)
