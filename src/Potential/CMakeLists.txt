
target_sources(lsmscore
        PUBLIC
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
        libxcInterface.cpp
        libxcInterface.hpp
        newexchg.f
        newpot_c.f
        PotentialShifter.hpp
        PotentialShifter.cpp
        rs.f
        )

target_include_directories(
        lsmscore
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)


