

target_sources(lsmscore PUBLIC
        calculateMultipoleMadelung.hpp
        calculateMultipoleMadelung.cpp
        MultipoleMadelung.cpp
        MultipoleMadelung.hpp
        lattice_utils.cpp
        lattice_utils.hpp
        madelung.cpp
        madelung.hpp
        madelung_term.hpp
        integer_factors.hpp
        integer_factors.cpp
        utils.hpp debug.hpp
        common.hpp
        spherical_harmonics.hpp
        spherical_harmonics.f90
        gaunt_factor.hpp
        gaunt_factor.cpp
        gaunt_factor.f90
        integer_factors.f90
        legendre.f90
        clm.f90
        GauntFactor.cpp
        GauntFactor.hpp
        gaussq.f90)

target_include_directories(lsmscore PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
