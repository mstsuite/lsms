#ifndef ME_CUDACHECKERROR_H
#define ME_CUDACHECKERROR_H

#if defined(ACCELERATOR_CUDA_C) && defined(ACCELERATOR_HIP)
#error Both ACCELERATOR_CUDA_C and ACCELERATOR_HIP defined. Only ONE of these allowed!
#endif

// inline int cudaCheckError(void) {return 0;}
#include <cstdio>

#if defined(ACCELERATOR_CUDA_C)
#define deviceCheckError() {                                          \
 cudaError_t e=cudaGetLastError();                                 \
 if(e!=cudaSuccess) {                                              \
   printf("Cuda failure %s:%d: '%s'\n",__FILE__,__LINE__,cudaGetErrorString(e));           \
   exit(0); \
 }                                                                 \
}

#define cublasCheckError(e) {                                          \
 if(e!=CUBLAS_STATUS_SUCCESS) {                                              \
   printf("CUBLAS failure %s:%d:\n",__FILE__,__LINE__);           \
   exit(0); \
 }                                                                 \
}
#endif
#if defined(ACCELERATOR_HIP)
#define deviceCheckError() {                                          \
 hipError_t e=hipGetLastError();                                 \
 if(e!=hipSuccess) {                                              \
   printf("HIP failure %s:%d: '%s'\n",__FILE__,__LINE__,hipGetErrorString(e));           \
   exit(0); \
 }                                                                 \
}

#define cublasCheckError(e) {}
//								       \
  // if(e!=CUBLAS_STATUS_SUCCESS) {				       \
  // printf("CUBLAS failure %s:%d:\n",__FILE__,__LINE__);           \
  // exit(0); \
//  }                                                                 \
}
#endif

#endif
