

# 1. Test Lua Libraries
add_executable(test_lua_library test_lua_library.cpp)
target_link_libraries(test_lua_library PUBLIC Lua::Lua)
target_link_libraries(test_lua_library PUBLIC gtest_main)

# Discover all tests
include(GoogleTest)
gtest_discover_tests(test_lua_library)