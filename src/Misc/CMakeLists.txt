


target_sources(lsmscore
        PUBLIC

        SphericalHarmonics.cpp
        SphericalHarmonics.hpp

        associatedLegendreFunction.cpp
        associatedLegendreFunction.hpp
        bulirsch_stoer.f
        bulirschStoerIntegrator.cpp
        bulirschStoerIntegrator.hpp
        calculateGauntCoeficients.cpp
        calculateGaussLegendrePoints.hpp

        cgaunt_c.f90

        clock_time.c
        ClebschGordan.hpp
        ClebschGordan.cpp
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
        Indices.cpp

        spherical_harmonics.hpp
        spherical_harmonics.f90

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
        diff.hpp

        PotentialType.hpp
        num_digits.hpp

        ChargePlotter.cpp
        ChargePlotter.hpp

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


