


target_sources(lsmscore
        PUBLIC
        associatedLegendreFunction.cpp
        associatedLegendreFunction.hpp
        bulirsch_stoer.f
        bulirschStoerIntegrator.cpp
        bulirschStoerIntegrator.hpp
        calculateGauntCoeficients.cpp
        calculateGaussLegendrePoints.hpp
        cgaunt_c.f
        clock_time.c
        cgc.h
        cinterp.f
        clebsch.f
        cmtruni.f
        cnewint.f
        Coeficients.cpp
        Coeficients.hpp
        congauss_c.f
        constraint.f
        dfv.f
        dfv_new.f
        fit.f
        fitpot.f
        fnpi.f
        fstop.f
        gaunt.f
        Gaunt_gmp.hpp
        Gaunt.hpp
        getclm.f
        ifacts_c.f
        Indices.hpp
        initwave.f
        integrateOneDim.hpp
        interp.f
        matr.f
        matrot1.f
        mbeqa.f
        mod_midpoint.f
        newder.f
        newint.f
        plglmax.f
        quadrature.cpp
        rationalFit.hpp
        readLastLine.cpp
        readLastLine.hpp
        ricbes.f
        rotmat.f
        rsimp.f
        rwave.f
        rzextr.f
        spin_trafo.f
        stop_with_backtrace.hpp
        trltog.f
        u_sigma_u.f
        v_plus_minus.f
        wrtmtx.f
        ylag.f
        zeroout.f
        zsphbes.f
        zsphbesj.f
        zsphbesjh.f
        io.hpp

        utils.hpp
        integrator.hpp
        poisson.hpp
        )

if(USE_ESSL)
  target_sources(lsmscore
        PUBLIC
        essl_workaround.f
        )
  message(STATUS "Adding essl_workaround.f")
endif()

target_include_directories(
        lsmscore
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)


