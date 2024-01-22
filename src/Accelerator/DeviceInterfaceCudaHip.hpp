/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */
#ifndef LSMS_DEVICEINTERFACECUDAHIP_HPP
#define LSMS_DEVICEINTERFACECUDAHIP_HPP

#if defined(ACCELERATOR_CUDA_C) && defined(ACCELERATOR_HIP)
#error Both ACCELERATOR_CUDA_C and ACCELERATOR_HIP defined. Only ONE of these allowed!
#endif

#if defined(ACCELERATOR_CUDA_C)

#include <cuda.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>

typedef cuDoubleComplex deviceDoubleComplex;

const auto deviceMemcpyHostToDevice = cudaMemcpyHostToDevice;
const auto deviceMemcpyDeviceToHost = cudaMemcpyDeviceToHost;

typedef cudaStream_t deviceStream_t;
typedef cudaError_t deviceError_t;
typedef cudaEvent_t deviceEvent_t;

const auto deviceSuccess = cudaSuccess;

inline deviceError_t deviceMalloc (void** devPtr, size_t size )
  {return cudaMalloc (devPtr, size);}

inline deviceError_t deviceFree (void* devPtr)
{return cudaFree(devPtr);}

inline deviceError_t deviceMallocHost (void** devPtr, size_t size )
  {return cudaMallocHost (devPtr, size);}

inline deviceError_t deviceFreeHost (void* devPtr)
{return cudaFreeHost(devPtr);}

inline deviceError_t deviceMemcpy (void* dst, const void* src, size_t count, cudaMemcpyKind kind)
  {return cudaMemcpy(dst, src, count, kind);}

inline deviceError_t deviceMemcpyAsync (void* dst, const void* src, size_t count, cudaMemcpyKind kind, cudaStream_t stream = 0)
  {return cudaMemcpyAsync (dst, src, count, kind, stream);}

const auto deviceEventDefault = cudaEventDefault;
const auto deviceEventBlockingSync = cudaEventBlockingSync;
const auto deviceEventDisableTiming = cudaEventDisableTiming;
const auto deviceEventInterprocess = cudaEventInterprocess;

inline deviceError_t deviceEventCreateWithFlags(cudaEvent_t *e, unsigned f)
  {return cudaEventCreateWithFlags(e, f);}

inline deviceError_t deviceEventDestroy(cudaEvent_t e)
  {return cudaEventDestroy(e);}

inline deviceError_t deviceStreamCreate(cudaStream_t *stream)
  {return cudaStreamCreate(stream);}

inline deviceError_t deviceStreamDestroy(cudaStream_t stream)
  {return cudaStreamDestroy(stream);}

#endif
#if defined(ACCELERATOR_HIP)

#include <hip/hip_runtime.h>
#include <hip/hip_complex.h>
#include <hipblas.h>

typedef hipDoubleComplex deviceDoubleComplex;

const auto deviceMemcpyHostToDevice = hipMemcpyHostToDevice;
const auto deviceMemcpyDeviceToHost = hipMemcpyDeviceToHost;

typedef hipStream_t deviceStream_t;
typedef hipError_t deviceError_t;
typedef hipEvent_t deviceEvent_t;

const auto deviceSuccess = hipSuccess;

inline deviceError_t deviceMalloc (void** devPtr, size_t size )
  {return hipMalloc (devPtr, size);}

inline deviceError_t deviceMallocHost (void** devPtr, size_t size )
  { return hipHostMalloc(devPtr, size, hipHostMallocNumaUser); }
  // {return hipMallocHost (devPtr, size); }
  // /* {return hipHostMalloc(devPtr, size, 0); /* hipMallocHost (devPtr, size);  } */

inline deviceError_t deviceFree (void* devPtr)
  {return hipFree(devPtr);}

inline deviceError_t deviceFreeHost (void* devPtr)
  { return hipHostFree(devPtr); }
  // {return hipFreeHost(devPtr); }
  // /* {return hipHostFree(devPtr); /* hipFreeHost(devPtr);  } */

inline deviceError_t deviceMemcpy (void* dst, const void* src, size_t count, hipMemcpyKind kind)
  {return hipMemcpy(dst, src, count, kind);}

inline deviceError_t deviceMemcpyAsync (void* dst, const void* src, size_t count, hipMemcpyKind kind, hipStream_t stream = 0)
  {return hipMemcpyAsync (dst, src, count, kind, stream);}

const auto deviceEventDefault = hipEventDefault;
const auto deviceEventBlockingSync = hipEventBlockingSync;
const auto deviceEventDisableTiming = hipEventDisableTiming;
const auto deviceEventInterprocess = hipEventInterprocess;

inline deviceError_t deviceEventCreateWithFlags(hipEvent_t *e, unsigned f)
  {return hipEventCreateWithFlags(e, f);}

inline deviceError_t deviceEventDestroy(hipEvent_t e)
  {return hipEventDestroy(e);}

inline deviceError_t deviceStreamCreate(hipStream_t *stream)
  {return hipStreamCreate(stream);}	

inline deviceError_t deviceStreamDestroy(hipStream_t stream)
  {return hipStreamDestroy(stream);}	

#endif

#endif
