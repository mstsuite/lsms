#
# Created by Franco P. Moitzi
#


target_sources(
        lsmscore
        PUBLIC
        accelerator_initialize.F
        accelerator_finalize.F
        accelerator_interface.h
        Accelerator.hpp
        Accelerator.cpp
        zblock_lu_cuda_c.cu
        fortran.c
        fortran.h
        inverse.cu
        DeviceStorage.cu
        DeviceStorage.hpp
)


target_include_directories(
        lsmscore
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)