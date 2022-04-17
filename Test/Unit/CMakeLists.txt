
# Add gtest
include(FetchContent)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG release-1.11.0
)

set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)
include(GoogleTest)

include(FortranCInterface)
FortranCInterface_HEADER(FC.h MACRO_NAMESPACE "FC_")
FortranCInterface_VERIFY(CXX)

##################################################################
# 1. Test Lua Libraries
add_executable(test_lua_library test_lua_library.cpp)
target_link_libraries(test_lua_library PUBLIC Lua::Lua)
add_test(
        NAME test_lua_library
        COMMAND $<TARGET_FILE:test_lua_library>
)
##################################################################

##################################################################
# 2. Test XC
add_executable(test_xc test_xc.cpp)
target_link_libraries(test_xc PUBLIC gtest_main)
target_sources(test_xc PUBLIC
        ${CMAKE_SOURCE_DIR}/src/Potential/alpha2_c.f
        ${CMAKE_CURRENT_BINARY_DIR}/FC.h
        xc.f90
        xc.hpp
        helpers.hpp
        )

target_link_libraries(test_xc PRIVATE lsmscore)

target_include_directories(test_xc PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(test_xc PUBLIC libxc::libxc)

gtest_discover_tests(test_xc)
##################################################################

##################################################################
# 3. Test Integration Routines
add_executable(test_integrator test_integrator.cpp)
target_link_libraries(test_integrator PUBLIC lsmscore)
add_test(
        NAME test_integrator
        COMMAND $<TARGET_FILE:test_integrator>
)
##################################################################

###################################################################
# 4. Test Poisson Solver Routines
add_executable(test_poisson test_poisson.cpp)
target_link_libraries(test_poisson PUBLIC lsmscore)
add_test(
        NAME test_poisson
        COMMAND $<TARGET_FILE:test_poisson>
)
###################################################################