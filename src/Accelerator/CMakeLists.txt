#
# Created by Franco P. Moitzi
#
if (USE_ACCELERATOR_CUDA_C)
    target_sources(
            lsmscore
            PUBLIC

            inverse.cu
            inverse.h

            operations.h

            DeviceStorage.cu
            DeviceStorage.hpp

            zblock_lu_cuda_c.cu

            cudaDoubleComplex.hpp
            deviceCheckError.hpp

            )
endif()

if (USE_ACCELERATOR_HIP)
    target_sources(
            lsmscore
            PUBLIC
            DeviceStorage_HIP.cpp
            DeviceStorage.hpp
            deviceCheckError.hpp
    )
endif ()

target_sources(
        lsmscore
        PUBLIC
        accelerator_initialize.F
        accelerator_finalize.F
        Accelerator.cpp
        Accelerator.hpp
        fortran.c
        fortran.h
)

if (USE_ACCELERATOR_HIP)
    target_sources(
            lsmscore
            PUBLIC
    )
endif ()

target_include_directories(
        lsmscore
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)
