
add_library(mixing_ref SHARED)

target_sources(mixing_ref
        PUBLIC
        nocaseCompare.F90
        MixingModule.F90
        KindParamModule.F90
        MathParamModule.F90
        PublicParamDefinitionsModule.F90
        PublicTypeDefinitionsModule.F90
        )

set_target_properties(mixing_ref PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mods)

target_include_directories(mixing_ref PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/mods)

target_compile_options(mixing_ref PUBLIC -fcheck=all -finit-local-zero)


add_executable(test_mixing test_mixing.f90)

target_link_libraries(test_mixing PUBLIC mixing_ref)
target_link_libraries(test_mixing PUBLIC LAPACK::LAPACK)