
if (NOT USE_ACCELERATOR_CUDA_C AND NOT USE_ACCELERATOR_HIP)

    message(STATUS "Regression tests are used")

    find_package(UnixCommands)

    list(APPEND TESTS
            01_al_fcc_asa
            02_al_fcc_conv_asa
            03_fe_bcc_asa
            04_fe_bcc_conv_asa
            05_cuzn_b2_asa)

    configure_file(
            ${CMAKE_CURRENT_SOURCE_DIR}/test_regression.py.in
            ${CMAKE_CURRENT_BINARY_DIR}/test_regression.py)

    foreach (_tests IN LISTS TESTS)

        file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${_tests}
                DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
                )

    endforeach ()

    find_package(Python COMPONENTS Interpreter REQUIRED)

    add_test(NAME RegressionTests
            COMMAND "${Python_EXECUTABLE} -m pytest -s"
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            )

else ()

    message(STATUS "Regression tests are not used")

endif ()