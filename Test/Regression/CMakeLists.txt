
#add_executable(lsms_comparer compare.cpp)
#set_target_properties(lsms_comparer
#        PROPERTIES
#        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
#        )

if (NOT USE_ACCELERATOR_CUDA_C AND NOT USE_ACCELERATOR_HIP)

    message(STATUS "Regression tests are used")

    find_package(UnixCommands)
    find_program(DIFF diff)

    list(APPEND TESTS 01_al_fcc_asa 01_al_fcc_asa_restart)

#    foreach(_tests IN LISTS TESTS)
#
#        file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${_tests}
#                DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
#                PATTERN test.sh.in
#                EXCLUDE
#                )
#
#        configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${_tests}/test.sh.in
#                ${CMAKE_CURRENT_BINARY_DIR}/${_tests}/test.sh.in)
#
#        file(GENERATE
#                OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${_tests}/test.sh
#                INPUT ${CMAKE_CURRENT_BINARY_DIR}/${_tests}/test.sh.in)
#
#        add_test(NAME regtest_${_tests}
#                COMMAND
#                ${BASH} test.sh
#                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${_tests}
#                )
#
#    endforeach()


else ()

    message(STATUS "Regression tests are not used")

endif ()