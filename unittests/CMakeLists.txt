

#function(add_test test_name library_target)

    #add_executable(${testname} ${testname})

    #target_link_libraries(test_math PUBLIC gtest_main)

    #if (TARGET ${library_target})

    #else ()
    #    message(FATAL_ERROR "Target doesn't exists")
    #endif ()

    #get_target_property(LIB_INCLUDES ${library_target} INCLUDE_DIRECTORIES)
    #if (_LIB_INCLUDES)
    #    foreach (dir ${LIB_INCLUDES})
    #        target_include_directories(${exename} ${dir})
    #    endforeach ()
    #endif ()


    #target_include_directories(${exename} PRIVATE $<TARGET_PROPERTY:${library_target},Fortran_MODULE_DIRECTORY>)
    #target_link_directories(${exename} PRIVATE $<TARGET_PROPERTY:${library_target},LIBRARY_OUTPUT_DIRECTORY>)
    #target_link_libraries(${exename} PRIVATE -l${library_target})
    #add_dependencies(${exename} ${library_target})

    #get_target_property(_compile_defs ${library_target} INTERFACE_COMPILE_DEFINITIONS)
    # if (_compile_defs)
    #     target_compile_definitions(${exename} PUBLIC ${_compile_defs})
    # endif ()

    #get_target_property(_compile_opts ${library_target} INTERFACE_COMPILE_OPTIONS)
    #if (_compile_opts)
    #target_compile_options(${exename} PUBLIC ${_compile_opts})
    #endif ()

#endmacro()


# 1. Test Lua Libraries
add_executable(test_lua_library test_lua_library.cpp)
target_link_libraries(test_lua_library Lua::Lua)
target_link_libraries(test_lua_library test_main)

# 2. Test
add_executable(test_math test_math.cpp)
target_link_libraries(test_math lsmscore)
target_link_libraries(test_math gtest_main)


# Discover all tests
include(GoogleTest)
gtest_discover_tests(test_lua_library)
gtest_discover_tests(test_math)
