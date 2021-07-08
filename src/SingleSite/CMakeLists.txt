# Sources ----------------------------------------------------------------------
set(C_SRC
    SingleSiteScattering.cpp
    readSingleAtomData_hdf5.cpp
    readSingleAtomData_bigcell.cpp
    writeSingleAtomData_hdf5.cpp
    writeSingleAtomData_bigcell.cpp
    checkAntiFerromagneticStatus.cpp
)

set(F_SRC
    F_readSingleAtomData_bigcell.f90
    single_scatterer_nonrel.f
    semrel_c.f
    scalar_m.f
    single_scatterer_rel.f
    spzwafu.f
    csbf.f
    matops.f
    gjinv.f
    dirmag1-op.f
    dirmag2-op.f
    brmat.f
    F_writeSingleAtomData_bigcell.f90
)

# Libraries --------------------------------------------------------------------
add_library(SingleSite STATIC ${C_SRC} ${F_SRC})

# Includes ---------------------------------------------------------------------
target_include_directories(SingleSite PRIVATE 
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/src
)
