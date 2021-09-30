

###################################################################
# 1. Test Lua Libraries
add_executable(test_lua_library test_lua_library.cpp)
target_link_libraries(test_lua_library PUBLIC Lua::Lua)
add_test(
        NAME test_lua_library
        COMMAND $<TARGET_FILE:test_lua_library>
)
###################################################################
###################################################################
# 2. Test Integration Routines
add_executable(test_integrator test_integrator.cpp)
target_link_libraries(test_integrator PUBLIC lsmscore)
add_test(
        NAME test_integrator
        COMMAND $<TARGET_FILE:test_integrator>
)
###################################################################
###################################################################
# 3. Test Poisson Solver Routines
add_executable(test_poisson test_poisson.cpp)
target_link_libraries(test_poisson PUBLIC lsmscore)
add_test(
        NAME test_poisson
        COMMAND $<TARGET_FILE:test_poisson>
)
###################################################################