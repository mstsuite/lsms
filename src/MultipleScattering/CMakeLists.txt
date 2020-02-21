set(SRC
    calculateTauMatrix.cpp
    makegij_c.f
    setgij.f
    block_inverse_fortran.f
    zblock_lu.F
    wasinv.f
    zmar1.f
    wasinv_p.f
    zuqmx.f
    zutfx.f
    zucpx.f
    zaxpby.f
    zrandn.f
    tau_inv_postproc.f
    trgtol.f
    green_function.f
    gf_local.f
    int_zz_zj.f
    mdosms_c.f
    mgreen_c.f
    green_function_rel.f
    write_kkrmat.f
    relmtrx.f
    gfill.f
    gafill.f
    magnet.f
    magnetic_dens.f
    new_dens.f
    block_inverse.cpp
    zblock_lu_cpp.cpp
    zblock_lu_cublas.cpp
)

# Target -----------------------------------------------------------------------
add_library(MultipleScattering STATIC ${SRC})

# Libraries --------------------------------------------------------------------
target_link_libraries(MultipleScattering
    PUBLIC
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
)

if(ENABLE_OPENMP AND OpenMP_CXX_FOUND)
    target_link_libraries(MultipleScattering PUBLIC OpenMP::OpenMP_CXX)
endif()

target_include_directories(MultipleScattering PRIVATE 
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/src
)
