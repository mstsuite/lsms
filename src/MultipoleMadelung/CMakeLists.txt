

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

        utils.hpp
        debug.hpp
        common.hpp)




target_include_directories(lsmscore PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
