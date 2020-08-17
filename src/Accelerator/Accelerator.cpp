/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */

#include "Accelerator.hpp"
#include <iostream>
#if defined(ACCELERATOR_CUDA_C)
#include <cuda_runtime.h>
#endif
#if defined(ACCELERATOR_HIP)
#include <hip_runtime.h>
#endif

void acceleratorInitialize(int sz, int nthreads)
{
  accelerator_initialize_(&sz);
}

void acceleratorFinalize(void)
{
  accelerator_finalize_();
#if defined(ACCELERATOR_CUDA_C)
  cudaDeviceReset();
#endif
#if defined(ACCELERATOR_HIP)
  hipDeviceReset();
#endif
}

void acceleratorPrint(void)
{
#if defined(ACCELERATOR_CUDA_C)
  cudaError_t cuda_error;
  int cuda_count;
  cudaDeviceProp cuda_prop;
  cuda_error = cudaGetDeviceCount(&cuda_count);
  std::cout << "Found " << cuda_count << " CUDA GPUs." << std::endl;
  for (int i = 0; i < cuda_count; i++)
  {
    cuda_error = cudaGetDeviceProperties(&cuda_prop, i);
    std::cout << "Device " << i << ": " << cuda_prop.name << std::endl;
  }
#endif
#if defined(ACCELERATOR_HIP)
  hipError_t hipError;
  int deviceCount;
  hipDeviceProp_t hipProp;
  hipError = hipGetDeviceCount(&deviceCount);
  std::cout << "Found " << deviceCount << " HIP GPUs." << std::endl;
  for(int i=0; i<deviceCount; i++)
  {
    hipError = hipGetDeviceProperties(&hipProp, i);
    std::cout << "Device " << i << ": " << hipProp.name << std::endl;
  }
#else
  ;
#endif
}
