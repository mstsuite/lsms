

target_sources(lsmscore
        PUBLIC
        calculateCoreStates.hpp
        calculateCoreStates.cpp
        coreSolver.cpp
        coreSolver.hpp
        corslv_c.f
        deepst_c.f
        getcor_c.f
        invals_c.f
        inwhnk_c.f
        inws_c.f
        outws_c.f
        richnk_c.f
        semcst_c.f
        )

target_include_directories(lsmscore PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
