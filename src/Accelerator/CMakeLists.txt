# Sources ----------------------------------------------------------------------
set(SRC
    accelerator_initialize.F
    accelerator_finalize.F
    accelerator_interface.c
    Accelerator.cpp
    fortran.c
)

if(CUDA_CXX)
    list(APPEND SRC
        buildKKRMatrix_gpu.cu
        buildKKRMatrix_kernels.cu
        zblock_lu_cuda_c.cu
        DeviceStorage.cu
        inverse.cu
    )
else()
    list(APPEND SRC
        no_gpu_buildKKRMatrix.cpp
        no_gpu_kernel_buildKKRMatrix.cpp
    )
endif()

# Target -----------------------------------------------------------------------
add_library(Accelerator STATIC ${SRC})

# Libraries --------------------------------------------------------------------
#target_link_libraries(Accelerator
#    PUBLIC
#        ${BLAS_LIBRARIES}
#        ${LAPACK_LIBRARIES}
#)

# Includes ---------------------------------------------------------------------
target_include_directories(Accelerator PRIVATE 
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/src
)