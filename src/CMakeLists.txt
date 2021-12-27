
#
#
#
add_executable(lsms_main)
set_target_properties(lsms_main PROPERTIES OUTPUT_NAME "lsms")

set_target_properties(lsms_main
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
        )

add_library(lsmscore STATIC)

target_link_libraries(lsms_main
        PUBLIC lsmscore)

add_executable(lsms_wl)
set_target_properties(lsms_wl PROPERTIES OUTPUT_NAME "wl-lsms")

set_target_properties(lsms_wl
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
        )

target_link_libraries(lsms_main
        PUBLIC lsmscore)
#
#
set_target_properties(lsmscore PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
        ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
        Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)


target_include_directories(
        lsmscore
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)


target_link_libraries(
        lsmscore
        PUBLIC
        gfortran
)

add_subdirectory(Accelerator)

add_subdirectory(Communication)

add_subdirectory(Core)

add_subdirectory(LuaInterface)

add_subdirectory(Main)

add_subdirectory(Madelung)

add_subdirectory(Misc)

add_subdirectory(MultipleScattering)

add_subdirectory(Potential)

add_subdirectory(RadialGrid)

add_subdirectory(SingleSite)

add_subdirectory(TotalEnergy)

add_subdirectory(VORPOL)


target_link_libraries(lsms_main PUBLIC lsmscore)
