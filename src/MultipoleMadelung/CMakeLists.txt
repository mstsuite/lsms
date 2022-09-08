

target_sources(lsmscore PUBLIC
        calculateMultipoleMadelung.hpp
        calculateMultipoleMadelung.cpp
        MultipoleMadelung.cpp
        MultipoleMadelung.hpp
        lattice_utils.cpp
        lattice_utils.hpp
        monopole_madelung.cpp
        monopole_madelung.hpp
        madelung_term.hpp
        madelung_term.cpp
        integer_factors.hpp
        integer_factors.cpp
        spherical_harmonics.f90
        spherical_harmonics.hpp
        utils.hpp
        debug.hpp
        common.hpp)




target_include_directories(lsmscore PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
