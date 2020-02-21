if (NOT "$ENV{OLCF_ESSL_ROOT}" STREQUAL "")
    set(ESSL_PATH $ENV{OLCF_ESSL_ROOT})
endif()

if (NOT "$ENV{ESSL_ROOT}" STREQUAL "")
    set(ESSL_PATH $ENV{ESSL_ROOT})
endif()

if (ESSL_PATH)
    find_path(ESSL_INCLUDE_DIR
        NAMES essl.h
        PATHS ${ESSL_PATH}/include
    )

    find_library(ESSL_LIBRARY
      NAMES essl
      PATHS ${ESSL_PATH}/lib64
    )

    find_package_handle_standard_args(ESSL
        FOUND_VAR ESSL_FOUND
        REQUIRED_VARS
            ESSL_LIBRARY
            ESSL_INCLUDE_DIR
    )

    if(ESSL_FOUND AND NOT TARGET ESSL)
        add_library(ESSL UNKNOWN IMPORTED)
        set_target_properties(ESSL PROPERTIES
            IMPORTED_LOCATION ${ESSL_LIBRARY}
            INTERFACE_INCLUDE_DIRECTORIES ${ESSL_INCLUDE_DIR}
        )
    endif()

    mark_as_advanced(
        ESSL_INCLUDE_DIR
        ESSL_LIBRARY
    )
else()
    message(FATAL_ERROR "Could not find ESSL, please set ESSL_ROOT.")
endif()

