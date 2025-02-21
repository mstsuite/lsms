
#
# Add GTest library
#

include(FetchContent)
set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
set(BUILD_GTEST ON CACHE BOOL "" FORCE)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG release-1.12.1
)

set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)
include(GoogleTest)

##################################################################
# 1. Test: Lua Libraries
add_executable(test_lua_library test_lua_library.cpp)
target_sources(test_lua_library PUBLIC accel_common.hpp)
target_include_directories(test_lua_library PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(test_lua_library PUBLIC Lua::Lua)
target_link_libraries(test_lua_library PUBLIC lsmscore)
target_link_libraries(test_lua_library PUBLIC gtest_main)

foreach (ref_file "input.lua")
    add_custom_command(
            TARGET test_lua_library
            POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy
            "${CMAKE_CURRENT_SOURCE_DIR}/${ref_file}"
            "${CMAKE_CURRENT_BINARY_DIR}/${ref_file}")
endforeach ()

gtest_discover_tests(test_lua_library)
##################################################################

##################################################################
# 2. Test: XC
if (USE_LIBXC)

    add_executable(test_xc test_xc.cpp)
    target_link_libraries(test_xc PUBLIC gtest_main)
    target_sources(test_xc PUBLIC
            ${CMAKE_SOURCE_DIR}/src/Potential/alpha2_c.f
            xc.f90
            xc.hpp
            helpers.hpp
            )

    target_link_libraries(test_xc PRIVATE lsmscore)
    target_include_directories(test_xc PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
    target_link_libraries(test_xc PUBLIC libxc::libxc)
    target_sources(test_xc PUBLIC accel_common.hpp)
    target_include_directories(test_xc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

    foreach (ref_file "result_pbe_libxc.out" "result_pbe.out")
        add_custom_command(
                TARGET test_xc
                POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy
                "${CMAKE_CURRENT_SOURCE_DIR}/${ref_file}"
                "${CMAKE_CURRENT_BINARY_DIR}/${ref_file}")
    endforeach ()

    gtest_discover_tests(test_xc)

endif ()
##################################################################

##################################################################
# 3. Test: Integration Routines
add_executable(test_integrator test_integrator.cpp)
target_sources(test_integrator PUBLIC accel_common.hpp)
target_include_directories(test_integrator PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(test_integrator PUBLIC lsmscore)
target_link_libraries(test_integrator PUBLIC gtest_main)

gtest_discover_tests(test_integrator)
##################################################################

###################################################################
# 4. Test: Poisson solver routines
add_executable(test_poisson test_poisson.cpp)
target_sources(test_poisson PUBLIC accel_common.hpp)
target_include_directories(test_poisson PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(test_poisson PUBLIC lsmscore)
target_link_libraries(test_poisson PUBLIC gtest_main)

gtest_discover_tests(test_poisson)
###################################################################

##################################################################
# 5. Test: Integration of local energy
add_executable(test_integrate_energy test_integrate_energy.cpp)
target_sources(test_integrate_energy PUBLIC accel_common.hpp)
target_include_directories(test_integrate_energy PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(test_integrate_energy PUBLIC gtest_main)
target_link_libraries(test_integrate_energy PRIVATE lsmscore)

foreach (ref_file "ezrho1.out" "ezrho2.out" "ezrho3.out")
    add_custom_command(
            TARGET test_integrate_energy
            POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy
            "${CMAKE_CURRENT_SOURCE_DIR}/${ref_file}"
            "${CMAKE_CURRENT_BINARY_DIR}/${ref_file}")
endforeach ()

gtest_discover_tests(test_integrate_energy)
##################################################################

##################################################################
# 6. Test: Differentiation
add_executable(test_diff test_diff.cpp)
target_sources(test_diff PUBLIC accel_common.hpp)
target_include_directories(test_diff PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(test_diff PUBLIC gtest_main)
target_link_libraries(test_diff PRIVATE lsmscore)

gtest_discover_tests(test_diff)
##################################################################

##################################################################
# 7. Test: Gaunt coefficients
add_executable(test_gaunt test_gaunt.cpp)
target_sources(test_gaunt PUBLIC accel_common.hpp)
target_include_directories(test_gaunt PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(test_gaunt PUBLIC gtest_main)
target_link_libraries(test_gaunt PRIVATE lsmscore)

gtest_discover_tests(test_gaunt)
##################################################################

##################################################################
# 8. Test: Spherical harmonics
add_executable(test_spherical_harmonics test_spherical_harmonics.cpp)
target_sources(test_spherical_harmonics PUBLIC accel_common.hpp)
target_include_directories(test_spherical_harmonics PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(test_spherical_harmonics PUBLIC gtest_main)
target_link_libraries(test_spherical_harmonics PRIVATE lsmscore)

gtest_discover_tests(test_spherical_harmonics)
##################################################################

##################################################################
# 9. Test: Madelung calculation
add_executable(test_madelung test_madelung.cpp)
target_link_libraries(test_madelung PUBLIC gtest_main)
target_link_libraries(test_madelung PRIVATE lsmscore)
target_sources(test_madelung PUBLIC accel_common.hpp)
target_include_directories(test_madelung PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

gtest_discover_tests(test_madelung)
##################################################################

##################################################################
# 10. Test: Multipole Madelung calculation
add_executable(test_multipole_madelung test_multipole_madelung.cpp)
target_link_libraries(test_multipole_madelung PUBLIC gtest_main)
target_link_libraries(test_multipole_madelung PRIVATE lsmscore)
target_sources(test_multipole_madelung PUBLIC accel_common.hpp)
target_include_directories(test_multipole_madelung PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

gtest_discover_tests(test_multipole_madelung)
##################################################################

##################################################################
# 11. Test: Eigensolver routines
add_executable(test_eigensolver test_eigensolver.cpp)
target_link_libraries(test_eigensolver PUBLIC gtest_main)
target_link_libraries(test_eigensolver PRIVATE lsmscore)
target_sources(test_eigensolver PUBLIC accel_common.hpp)
target_include_directories(test_eigensolver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

gtest_discover_tests(test_eigensolver)
##################################################################

##################################################################
# 12. Test: Mixer routines
add_executable(test_mixer test_mixer.cpp vector_fun.f90)
target_link_libraries(test_mixer PUBLIC gtest_main)
target_link_libraries(test_mixer PRIVATE lsmscore)
target_sources(test_mixer PUBLIC accel_common.hpp)
target_include_directories(test_mixer PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

gtest_discover_tests(test_mixer)
##################################################################

##################################################################
# 13. Test: Mixer routines
add_executable(test_atomicdft test_atomicdft.cpp)
target_link_libraries(test_atomicdft PUBLIC gtest_main)
target_link_libraries(test_atomicdft PRIVATE lsmscore)
target_sources(test_atomicdft PUBLIC accel_common.hpp)
target_include_directories(test_atomicdft PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

gtest_discover_tests(test_atomicdft)
##################################################################

##################################################################
# 13. Test: Single site routines
add_executable(test_ss test_single_site.cpp)
target_link_libraries(test_ss PUBLIC gtest_main)
target_link_libraries(test_ss PRIVATE lsmscore)
target_sources(test_ss PUBLIC accel_common.hpp)
target_include_directories(test_ss PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

FetchContent_Declare(matplotplusplus
        GIT_REPOSITORY https://github.com/alandefreitas/matplotplusplus
        GIT_TAG origin/master) # or whatever tag you want


FetchContent_GetProperties(matplotplusplus)
if(NOT matplotplusplus_POPULATED)
    FetchContent_Populate(matplotplusplus)
    add_subdirectory(${matplotplusplus_SOURCE_DIR} ${matplotplusplus_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

target_link_libraries(test_ss PUBLIC matplot)


gtest_discover_tests(test_ss)
##################################################################