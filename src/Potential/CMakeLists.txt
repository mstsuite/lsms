
target_sources(lsmscore
        PUBLIC

        xc_lda_potential.hpp

        alpha2_c.f
        calculateASAro3.f
        calculateChargesPotential.cpp
        calculateChargesPotential.hpp
        epcorr.f
        getqm_mt.f
        getvmt.cpp
        getvmt.hpp
        getXCName.cpp
        getXCName.hpp
        interpolatePotential.cpp
        interpolatePotential.hpp

        newexchg.f
        newpot_c.f
        PotentialShifter.hpp
        PotentialShifter.cpp
        rs.f
        lsf_functional.hpp
        lsf_functional.cpp
        newFunctionalInterface.cpp
        newFunctionalInterface.hpp
        common.hpp

        newPotential.cpp
        newPotential.hpp

        XCLibxc.cpp
        XCLibxc.hpp
        XCBase.cpp
        XCBase.hpp
        XCLDA.cpp
        XCLDA.hpp

        )

target_include_directories(
        lsmscore
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)


