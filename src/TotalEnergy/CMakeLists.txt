

target_sources(lsmscore
        PUBLIC
        calculateTotalEnergy.cpp
        calculateTotalEnergy.hpp
        localTotalEnergy.cpp
        localTotalEnergy.hpp
        zeropt_c.f
        janake_c.f
        DFTEnergy.cpp
        DFTEnergy.hpp
        )

target_include_directories(
        lsmscore
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)