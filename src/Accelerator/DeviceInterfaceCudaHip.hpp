/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */
#ifndef LSMS_DEVICEINTERFACECUDAHIP_HPP
#define LSMS_DEVICEINTERFACECUDAHIP_HPP

#if defined(ACCELERATOR_CUDA_C) && defined(ACCELERATOR_HIP)
#error Both ACCELERATOR_CUDA_C and ACCELERATOR_HIP defined. Only ONE of these allowed!
#endif

#if defined(ACCELERATOR_CUDA_C)
const auto deviceMemcpyHostToDevice = cudaMemcpyHostToDevice;

typedef cudaStream_t deviceStream_t;

inline cudaError_t deviceMalloc (void** devPtr, size_t size )
  {return cudaMalloc (devPtr, size);}

inline cudaError_t deviceFree (void* devPtr)
  {return cudaFree)devPtr);}

inline cudaError_t deviceMemcpy (void* dst, const void* src, size_t count, cudaMemcpyKind kind)
  {​return cudaMemcpy ( dst, src, count, kind);}

​inline cudaError_t deviceMemcpyAsync (void* dst, const void* src, size_t count, cudaMemcpyKind kind, cudaStream_t stream = 0)
  {​return cudaMemcpyAsync (dst, src, count, kind, stream);}

#endif
#if defined(ACCELERATOR_HIP)
const auto deviceMemcpyHostToDevice = hipMemcpyHostToDevice;

typedef hipStream_t deviceStream_t;

inline hipError_t deviceMalloc (void** devPtr, size_t size )
  {return hipMalloc (devPtr, size);}

inline cudaError_t deviceFree (void* devPtr)
  {return hipFree{devPtr);}

inline hipError_t deviceMemcpy (void* dst, const void* src, size_t count, hipMemcpyKind kind)
  {​return cudaMemcpy ( dst, src, count, kind);}

inline hipError_t deviceMemcpyAsync (void* dst, const void* src, size_t count, hipMemcpyKind kind, hipStream_t stream = 0)
  {​return hipMemcpyAsync (dst, src, count, kind, stream);}

#endif

#endif
