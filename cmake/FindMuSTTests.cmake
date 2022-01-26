#
# Created by Franco P. Moitzi
#

macro(add_must_test library_target testname)


    set(exename test_${testname})

    add_executable(${exename} ${testname}.cpp ${TEST_EXTRA_SOURCE})

    set_target_properties(${exename} PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests)


    foreach (objectlib ${TEST_EXTRA_OBJECTS})
        add_dependencies(${exename} ${objectlib})
        target_sources(${exename} PUBLIC $<TARGET_OBJECTS:${objectlib}>)
        target_include_directories(${exename}
                PRIVATE
                $<TARGET_PROPERTY:${objectlib},Fortran_MODULE_DIRECTORY>)
    endforeach ()

    get_target_property(LIB_INCLUDES lsmscore INCLUDE_DIRECTORIES)
    foreach (dir ${LIB_INCLUDES})
        target_include_directories(${exename} PUBLIC ${dir})
    endforeach ()

    target_include_directories(${exename} PRIVATE $<TARGET_PROPERTY:${library_target},Fortran_MODULE_DIRECTORY>)
    target_link_directories(${exename} PRIVATE $<TARGET_PROPERTY:${library_target},LIBRARY_OUTPUT_DIRECTORY>)
    target_link_libraries(${exename} PRIVATE -l${library_target})
    add_dependencies(${exename} ${library_target})

    get_target_property(_compile_defs ${library_target} INTERFACE_COMPILE_DEFINITIONS)
    if (_compile_defs)
        target_compile_definitions(${exename} PUBLIC ${_compile_defs})
    endif ()

    get_target_property(_compile_opts ${library_target} INTERFACE_COMPILE_OPTIONS)
    if (_compile_opts)
        target_compile_options(${exename} PUBLIC ${_compile_opts})
    endif ()

    target_compile_options(${exename} PUBLIC ${TEST_FLAGS})

    if (WITH_MPI)
        target_compile_definitions(${exename} PUBLIC MPI)
        target_link_libraries(${exename} PRIVATE MPI::MPI_Fortran)
    endif ()

    add_test(NAME ${testname}
            COMMAND $<TARGET_FILE:${exename}>
            WORKING_DIRECTORY
            ${CMAKE_BINARY_DIR}/tests)

    # File with input data for the test
    set(init_file "${testname}.init")
    if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${init_file}")
        add_custom_command(
                TARGET ${exename} POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy
                "${CMAKE_CURRENT_SOURCE_DIR}/${init_file}"
                "${CMAKE_BINARY_DIR}/tests/${init_file}")
    endif ()

    # File with reference results
    set(ref_file "${testname}.ref")
    if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${ref_file}")
        add_custom_command(
                TARGET ${exename} POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy
                "${CMAKE_CURRENT_SOURCE_DIR}/${ref_file}"
                "${CMAKE_BINARY_DIR}/tests/${ref_file}")
    endif ()


endmacro()


macro(add_fortran_object_library lib main_lib sources)
    add_library(${lib} OBJECT ${sources})
    add_dependencies(${lib} ${main_lib})
    target_include_directories(${lib} PRIVATE $<TARGET_PROPERTY:${main_lib},INTERFACE_INCLUDE_DIRECTORIES>)
    target_include_directories(${lib} PRIVATE $<TARGET_PROPERTY:${main_lib},Fortran_MODULE_DIRECTORY>)
    set_target_properties(${lib} PROPERTIES
            Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/lib/test_modules
            )
endmacro()
