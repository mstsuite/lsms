/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */

#ifndef LSMS_DEVICE_STORAGE_HPP
#define LSMS_DEVICE_STORAGE_HPP

#include "Real.hpp"
#include "Complex.hpp"

#if defined(ACCELERATOR_CUDA_C)
#include <cublas_v2.h>
#include <cusolverDn.h>
#endif

#if defined(ACCELERATOR_HIP)
#include <hip/hip_runtime_api.h> 
#include <rocsolver.h>
#endif

#include "SingleSite/AtomData.hpp"

#ifdef _OPENMP
#include <omp.h>
#else
#ifndef LSMS_DUMMY_OPENMP
#define LSMS_DUMMY_OPENMP
inline int omp_get_max_threads() {return 1;}
inline int omp_get_num_threads() {return 1;}
inline int omp_get_thread_num() {return 0;}
#endif
#endif

// #include "DeviceMatrix.hpp"

template <class T> class DeviceMatrix;

#if defined(ACCELERATOR_CUDA_C)
extern "C" Complex* get_dev_m_();
extern "C" Complex* get_dev_bgij_();
extern "C" Complex* get_dev_tmat_n_();
extern "C" int* get_dev_ipvt_();
extern "C" cudaStream_t get_stream_(const int &id);
extern "C" cublasHandle_t get_cublas_handle_();
extern "C" cudaEvent_t get_cuda_event_();
extern "C" Complex* get_host_m_(const int &max_nrmat_ns);
#endif

static const int MAX_THREADS=16;
class DeviceStorage {
private:
  static int nThreads;
  static Complex *dev_m[MAX_THREADS], *dev_bgij[MAX_THREADS], *dev_tmat_n[MAX_THREADS];
  static Complex *dev_tau[MAX_THREADS], *dev_tau00[MAX_THREADS], *dev_t0[MAX_THREADS], *dev_t[MAX_THREADS];
  static int *dev_ipvt[MAX_THREADS];
#if defined(ACCELERATOR_CUDA_C)
  static cublasHandle_t cublas_h[MAX_THREADS];
  static cusolverDnHandle_t cusolverDnHandle[MAX_THREADS];  
  static cudaEvent_t event[MAX_THREADS];
  static cudaStream_t stream[MAX_THREADS][2];
#endif
#if defined (ACCELERATOR_HIP)
  static rocblas_handle rocblas_h[MAX_THREADS];
#endif
  static size_t dev_workBytes[MAX_THREADS];
  static void *dev_work[MAX_THREADS];
  static DeviceMatrix<Complex> dev_tmat_store;
  static bool initialized;
public:
  int allocate(int kkrsz_max,int nspin, int numLIZ, int _nThreads);
  void free();

  static Complex* getDevM() { return dev_m[omp_get_thread_num()]; } 
  static Complex* getDevBGij() { if(!initialized) {printf("DeviceStorage not initialized\n"); exit(1);}
                                 return dev_bgij[omp_get_thread_num()]; } 
  static Complex* getDevTmatN() { return dev_tmat_n[omp_get_thread_num()]; } 
  static Complex* getDevTau() { return dev_tau[omp_get_thread_num()]; }
  static Complex* getDevTau00() { return dev_tau00[omp_get_thread_num()]; }
  static Complex* getDevT() { return dev_t[omp_get_thread_num()]; }
  static Complex* getDevT0() { return dev_t0[omp_get_thread_num()]; }
  static int* getDevIpvt() { return dev_ipvt[omp_get_thread_num()]; }
#if defined (ACCELERATOR_CUDA_C)
  static cudaStream_t getStream(int i) { return stream[omp_get_thread_num()][i]; }
  static cudaEvent_t getEvent() { return event[omp_get_thread_num()]; }
  static cublasHandle_t getCublasHandle() { return cublas_h[omp_get_thread_num()]; }
  static cusolverDnHandle_t getCusolverDnHandle() { return cusolverDnHandle[omp_get_thread_num()]; }
#endif
#if defined(ACCELERATOR_HIP)
  static rocblas_handle getRocBlasHandle() { return rocblas_h[omp_get_thread_num()]; }
#endif
  static size_t getDevWorkBytes() { return dev_workBytes[omp_get_thread_num()]; }
  static void *getDevWork() {  return dev_work[omp_get_thread_num()]; }
  static DeviceMatrix<Complex>* getDevTmatStore() { return &dev_tmat_store; }
};

DeviceMatrix<Complex>* get_dev_tmat_store();

void *allocateDStore(void);
void freeDStore(void * d_store);
int initDStore(void * d_store,int kkrsz_max, int nspin, int numLIZ, int nthreads);

class DeviceAtomCuda {
public:
  bool allocated;
  Real *LIZPos;
  int *LIZlmax;
  int *LIZStoreIdx;
  int numLIZ;
  
  int allocate(int lmax, int nspin, int _numLIZ);
  void copyFromAtom(AtomData &atom);
  void free();
};


class DeviceConstants {
  public:
//  DeviceConstants() { }
//  ~DeviceConstants() { }
  int *lofk;
  int *mofk;
  Complex *ilp1;
  // DeviceMatrix<Complex> illp;
  Complex* illp;
  // DeviceArray3d<Real> cgnt;
  Real* cgnt;

  int allocate();
  void free();
};

#endif
