// -*- mode: c++; -*-

#include <stdlib.h>
#include <iostream>

#include <hip/hip_runtime.h>

#include "Real.hpp"
#include "Complex.hpp"
#include "Matrix.hpp"

#include "DeviceMatrix.hpp"
#include "DeviceArray3d.hpp"
#include "DeviceVector.hpp"
#include "Main/SystemParameters.hpp"
#include "Misc/Coeficients.hpp"
#include "Misc/Indices.hpp"

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

#include "DeviceStorage.hpp"

// #include "cuda_error.h"
#include "deviceCheckError.hpp"

using namespace std;

//TODO move inside DeviceStorage?
//allocate a thread specific matrix on the host and pin its memory
extern "C"
Complex* get_host_m_(const int &max_nrmat_ns) {
  static Complex * m_v=0;
  static int cur_size=0;
  static deviceError_t pinned;

  /*
  if(cur_size<max_nrmat_ns) {

    //release previously allocated memory
    if(m_v!=0) {
      free(m_v);
    }

      fprintf(stderr, "Matrix not pinned\n");
      m_v = (Complex*)malloc(max_nrmat_ns*max_nrmat_ns*sizeof(Complex)*omp_get_max_threads());
    cur_size=max_nrmat_ns;
  }
  */


  if(cur_size<max_nrmat_ns) {

    //release previously allocated memory
    if(m_v!=0) {
      if(pinned) deviceFreeHost(m_v);
      else free(m_v);
    }

    //allocate new memory
    pinned = deviceMallocHost((void**)&m_v,max_nrmat_ns*max_nrmat_ns*sizeof(Complex)*omp_get_max_threads());

    if ( pinned != deviceSuccess )
    {
      fprintf(stderr, "Matrix not pinned\n");
      m_v = (Complex*)malloc(max_nrmat_ns*max_nrmat_ns*sizeof(Complex)*omp_get_max_threads());
    }
    cur_size=max_nrmat_ns;
  }

  return m_v; 
}

/*
static const int MAX_THREADS=16;
class DeviceStorage {
private:
  static int nThreads;
  static Complex *dev_m[MAX_THREADS], *dev_bgij[MAX_THREADS], *dev_tmat_n[MAX_THREADS];
  static Complex *dev_tau[MAX_THREADS], *dev_tau00[MAX_THREADS];
  static int *dev_ipvt[MAX_THREADS];
  static cublasHandle_t cublas_h[MAX_THREADS];
  static cusolverDnHandle_t cusolverDnHandle[MAX_THREADS];
  static cudaEvent_t event[MAX_THREADS];
  static cudaStream_t stream[MAX_THREADS][2];
  static size_t dev_workBytes[MAX_THREADS];
  static void *dev_work[MAX_THREADS];
  static DeviceMatrix<Complex> dev_tmat_store;
  static bool initialized;
public:
*/

int DeviceStorage::allocateAdditional(int kkrsz_max, int nspin, int numLIZ, int _nThreads) {
   // initializing a big tau and big T matrix, needed (only) for conductivity purposes
   int N = kkrsz_max * nspin * numLIZ;
   int nThreads = _nThreads;
   for (int i=0; i < nThreads; i++) {
     deviceError_t err;
     err = deviceMalloc((void **) &dev_tauFull[i],(size_t)N * (size_t)N * sizeof(Complex));
     if (err != deviceSuccess) {
        printf("failed to allocate dev_tauFull[%d], size=%d, err=%d\n",
               i, (size_t)N * (size_t)N * sizeof(Complex), err);
        exit(1);
     }
     err = deviceMalloc((void **) &dev_tFull[i], (size_t)N * (size_t)N * sizeof(Complex));
     if (err != deviceSuccess) {
        printf("failed to allocate dev_tFull[%d], size=%d, err=%d\n",
               i, (size_t)N * (size_t)N * sizeof(Complex), err);
        exit(1);
     }
   }
   return 0;
}

int DeviceStorage::allocate(int kkrsz_max,int nspin, int numLIZ, int _nThreads, int iskubo)
  {
    if(!initialized)
    {
      //printf("*************************************MEMORY IS BEING ALLOCATED\n");
      if(_nThreads>MAX_THREADS)
      {
        printf("nThreads (%d) in DeviceStorage::allocate exceeds MAX_THREADS (%d)\n",_nThreads,MAX_THREADS);
        printf("  change MAX_THREADS in src/Accelerator/DeviceStorage.cu and recompile!\n");
        exit(1);
      }
      nThreads=_nThreads;
      int N=kkrsz_max*nspin*numLIZ;
      // printf("DeviceStorage::alocate N=%d\n",N);
      for(int i=0;i<nThreads;i++)
      {
        deviceError_t err;
        err = deviceMalloc((void**)&dev_m[i],(size_t)N*(size_t)N*sizeof(Complex));
        if(err!=deviceSuccess)
        {
          printf("failed to allocate dev_m[%d], size=%zu, err=%d\n",
                i,(size_t)N*(size_t)N*sizeof(Complex),err);
          exit(1);
        }
        deviceMalloc((void**)&dev_ipvt[i],(size_t)N*sizeof(int));
	deviceMalloc((void**)&dev_info[i],(size_t)nThreads*sizeof(int));
	if (iskubo == 0){
	  err = deviceMalloc((void**)&dev_bgij[i],(size_t)N*(size_t)N*sizeof(Complex));
          if(err!=deviceSuccess)
          {
            printf("failed to allocate dev_bgij[%d], size=%zu, err=%d\n",
                  i,(size_t)N*(size_t)N*sizeof(Complex),err);
            exit(1);
          }
	}
#ifdef BUILDKKRMATRIX_GPU
        // cudaMalloc((void**)&dev_bgij[i],4*kkrsz_max*kkrsz_max*numLIZ*numLIZ*sizeof(Complex));
        deviceMalloc((void**)&dev_tmat_n[i],4*kkrsz_max*kkrsz_max*numLIZ*sizeof(Complex)); 
#endif
        deviceMalloc((void**)&dev_tau[i], 4*(size_t)N*kkrsz_max*sizeof(Complex));
        deviceMalloc((void**)&dev_tau00[i], 4*kkrsz_max*kkrsz_max*sizeof(Complex));
        deviceMalloc((void**)&dev_t[i], 4*N*kkrsz_max*sizeof(Complex));
        deviceMalloc((void**)&dev_t0[i], 4*kkrsz_max*kkrsz_max*sizeof(Complex));
        deviceStreamCreate(&stream[i][0]);
        deviceStreamCreate(&stream[i][1]);
        deviceEventCreateWithFlags(&event[i],deviceEventDisableTiming);
        hipblasCreate(&hipblas_h[i]);
        // cusolverDnCreate(&cusolverDnHandle[i]);

	int lWork = 0;
	// cusolverDnZgetrf_bufferSize(cusolverDnHandle[i], N, N,
	//			    (cuDoubleComplex *)dev_m[i], N, &lWork);
	dev_workBytes[i] = 0;
#ifndef ARCH_IBM
	// cusolverDnZZgesv_bufferSize(cusolverDnHandle[i], N, 2*kkrsz_max,
	//			    (cuDoubleComplex *)dev_m[i], N, dev_ipvt[i], (cuDoubleComplex *)dev_tau[i], N, (cuDoubleComplex *)dev_tau[i], N,
	//			    dev_work[i], &dev_workBytes[i]);
#endif

	dev_workBytes[i] = std::max(dev_workBytes[i]*sizeof(deviceDoubleComplex),
				    lWork*sizeof(deviceDoubleComplex));
	deviceMalloc((void**)&dev_work[i], dev_workBytes[i]);
        // printf("  dev_m[%d]=%zx\n",i,dev_m[i]);
      }
      deviceCheckError();
      initialized=true;
    }
    return 0;
  }
  
  void DeviceStorage::free()
  {
    if(initialized) {
   //     printf("*************************************MEMORY IS BEING FREED\n");
      // for(int i=0;i<omp_get_max_threads();i++)
      for(int i=0; i<nThreads; i++)
      {
        deviceFree(dev_m[i]);
        deviceFree(dev_ipvt[i]);
        deviceFree(dev_info[i]);
#ifdef BUILDKKRMATRIX_GPU
        deviceFree(dev_bgij[i]);
        deviceFree(dev_tmat_n[i]);
#endif
	deviceFree(dev_work[i]);
        deviceFree(dev_t0[i]);
        deviceStreamDestroy(stream[i][0]);
        deviceStreamDestroy(stream[i][1]);
        deviceEventDestroy(event[i]);
        hipblasDestroy(hipblas_h[i]);
      }
      // dev_tmat_store.clear();
      deviceFree(devTmatStore);
      deviceCheckError();
      initialized=false;
    }
  }

/*
  static Complex* getDevM() { return dev_m[omp_get_thread_num()]; } 
  static Complex* getDevBGij() { if(!initialized) {printf("DeviceStorage not initialized\n"); exit(1);}
                                 return dev_bgij[omp_get_thread_num()]; } 
  static Complex* getDevTmatN() { return dev_tmat_n[omp_get_thread_num()]; } 
  static Complex* getDevTau() { return dev_tau[omp_get_thread_num()]; }
  static Complex* getDevTau00() { return dev_tau00[omp_get_thread_num()]; }
  static int* getDevIpvt() { return dev_ipvt[omp_get_thread_num()]; } 
  static cudaStream_t getStream(int i) { return stream[omp_get_thread_num()][i]; }
  static cudaEvent_t getEvent() { return event[omp_get_thread_num()]; }
  static cublasHandle_t getCublasHandle() { return cublas_h[omp_get_thread_num()]; }
  static cusolverDnHandle_t getCusolverDnHandle() { return cusolverDnHandle[omp_get_thread_num()]; }
  static size_t getDevWorkBytes() { return dev_workBytes[omp_get_thread_num()]; }
  static void *getDevWork() {  return dev_work[omp_get_thread_num()]; }
  static DeviceMatrix<Complex>* getDevTmatStore() { return &dev_tmat_store; }
};
*/

int DeviceStorage::copyTmatStoreToDevice(Matrix<Complex> &tmatStore,
    int blkSize)
{ 
  if((tmatStoreSize > 0) && (tmatStoreSize < tmatStore.size()))
  { 
    deviceFree(devTmatStore);
    tmatStoreSize = 0;
  }
  if(tmatStoreSize == 0)
  { 
    deviceMalloc((void **)&devTmatStore, tmatStore.size()*sizeof(Complex));
    tmatStoreSize = tmatStore.size();
  }
  deviceMemcpy(devTmatStore, &tmatStore(0,0),
    tmatStore.size()*sizeof(Complex), deviceMemcpyHostToDevice);
  blkSizeTmatStore = blkSize;
  tmatStoreLDim = tmatStore.l_dim();

  return 0;
}

bool DeviceStorage::initialized = false;
Complex *DeviceStorage::dev_m[MAX_THREADS], *DeviceStorage::dev_bgij[MAX_THREADS], *DeviceStorage::dev_tmat_n[MAX_THREADS];
Complex *DeviceStorage::dev_tau[MAX_THREADS], *DeviceStorage::dev_tau00[MAX_THREADS];
Complex *DeviceStorage::dev_tauFull[MAX_THREADS], *DeviceStorage::dev_tFull[MAX_THREADS];
Complex *DeviceStorage::dev_t0[MAX_THREADS];
Complex *DeviceStorage::dev_t[MAX_THREADS];
void *DeviceStorage::dev_work[MAX_THREADS];
size_t DeviceStorage::dev_workBytes[MAX_THREADS];
int *DeviceStorage::dev_ipvt[MAX_THREADS];
int *DeviceStorage::dev_info[MAX_THREADS];
hipblasHandle_t DeviceStorage::hipblas_h[MAX_THREADS];
// cusolverDnHandle_t DeviceStorage::cusolverDnHandle[MAX_THREADS];
deviceEvent_t DeviceStorage::event[MAX_THREADS];
deviceStream_t DeviceStorage::stream[MAX_THREADS][2];
// DeviceMatrix<Complex> DeviceStorage::dev_tmat_store;
Complex *DeviceStorage::devTmatStore;
size_t DeviceStorage::tmatStoreSize = 0;
int DeviceStorage::blkSizeTmatStore = 0;
int DeviceStorage::tmatStoreLDim = 0;
int DeviceStorage::nThreads=1;
bool initialized;

std::vector<DeviceAtom> deviceAtoms;

// Device Atom
int DeviceAtom::allocate(int _lmax, int _nspin, int _numLIZ)
{
  if(allocated) free();
  allocated = true;
  numLIZ = _numLIZ;
  deviceMalloc((void**)&LIZPos,numLIZ*3*sizeof(Real));
  deviceMalloc((void**)&LIZlmax,numLIZ*sizeof(int));
  deviceMalloc((void**)&LIZStoreIdx,numLIZ*sizeof(int));
  
  return 0;
}

void DeviceAtom::free()
{
  if(allocated)
  {
    deviceFree(LIZPos);
    deviceFree(LIZlmax);
    deviceFree(LIZStoreIdx);
  }
  allocated = false;
}

void DeviceAtom::copyFromAtom(AtomData &atom)
{
  if(!allocated)
  {
    allocate(atom.lmax, atom.nspin, atom.numLIZ);
  }
  deviceMemcpy(LIZPos, &atom.LIZPos(0,0), atom.numLIZ*3*sizeof(Real), deviceMemcpyHostToDevice);
  deviceMemcpy(LIZlmax, &atom.LIZlmax[0], atom.numLIZ*sizeof(int), deviceMemcpyHostToDevice);
  deviceMemcpy(LIZStoreIdx, &atom.LIZStoreIdx[0], atom.numLIZ*sizeof(int), deviceMemcpyHostToDevice);
}

int *DeviceConstants::lofk;
int *DeviceConstants::mofk;
deviceDoubleComplex *DeviceConstants::ilp1;
// DeviceMatrix<Complex> illp(ndlj, ndlj);
deviceDoubleComplex* DeviceConstants::illp;
int DeviceConstants::ndlj_illp;
// DeviceArray3d<Real> cgnt(lmax+1,ndlj,ndlj);
Real* DeviceConstants::cgnt;
int DeviceConstants::ndlj_cgnt, DeviceConstants::lmaxp1_cgnt;

int DeviceConstants::allocate()
{
  ndlj_illp = IFactors::illp.l_dim();
  lmaxp1_cgnt = GauntCoeficients::cgnt.l_dim1();
  ndlj_cgnt = GauntCoeficients::cgnt.l_dim2();

  deviceMalloc((void**)&lofk, AngularMomentumIndices::lofk.size()*sizeof(int));
  deviceMalloc((void**)&mofk, AngularMomentumIndices::mofk.size()*sizeof(int));
  deviceMalloc((void**)&ilp1, IFactors::ilp1.size()*sizeof(deviceDoubleComplex));
  deviceMalloc((void**)&illp, IFactors::illp.size()*sizeof(deviceDoubleComplex));
  deviceMalloc((void**)&cgnt, GauntCoeficients::cgnt.size()*sizeof(double));

  deviceMemcpy(lofk, &AngularMomentumIndices::lofk[0], AngularMomentumIndices::lofk.size()*sizeof(int), deviceMemcpyHostToDevice);
  deviceMemcpy(mofk, &AngularMomentumIndices::mofk[0], AngularMomentumIndices::mofk.size()*sizeof(int), deviceMemcpyHostToDevice);
  deviceMemcpy(ilp1, &IFactors::ilp1[0], IFactors::ilp1.size()*sizeof(deviceDoubleComplex), deviceMemcpyHostToDevice);
  deviceMemcpy(illp, &IFactors::illp[0], IFactors::illp.size()*sizeof(deviceDoubleComplex), deviceMemcpyHostToDevice);
  deviceMemcpy(cgnt, &GauntCoeficients::cgnt[0], GauntCoeficients::cgnt.size()*sizeof(double), deviceMemcpyHostToDevice);

  return 0;
}

void DeviceConstants::free()
{
  deviceFree(lofk);
  deviceFree(mofk);
  deviceFree(ilp1);
  deviceFree(illp);
  deviceFree(cgnt);
}

/****************Fortran Interfaces*********************/
/*
extern "C"
Complex* get_dev_m_() {
  return DeviceStorage::getDevM();
}

extern "C"
Complex* get_dev_bgij_() {
  return DeviceStorage::getDevBGij();
}

extern "C"
Complex* get_dev_tmat_n_() {
  return DeviceStorage::getDevTmatN();
}

extern "C"
int* get_dev_ipvt_() {
  return DeviceStorage::getDevIpvt();
}

extern "C"
cudaStream_t get_stream_(const int &id) {
  return DeviceStorage::getStream(id);
}

extern "C"
cublasHandle_t get_cublas_handle_() {
  return DeviceStorage::getCublasHandle();
}

//allocate a thread specific event
extern "C"
cudaEvent_t get_cuda_event_() {
  return DeviceStorage::getEvent();
}
*/
/********************************************************/

// DeviceMatrix<Complex>* get_dev_tmat_store() {
//   return DeviceStorage::getDevTmatStore();
// }

void *allocateDStore(void)
{
  return static_cast<void *>(new DeviceStorage);
}

void freeDStore(void * d_store)
{
  static_cast<DeviceStorage*>(d_store)->free();
  delete static_cast<DeviceStorage*>(d_store);
}

int initDStore(void * d_store,int kkrsz_max, int nspin, int numLIZ, int nthreads)
{
  return (*static_cast<DeviceStorage*>(d_store)).allocate(kkrsz_max,nspin,numLIZ,nthreads);
}

// void copyTmatStoreToDevice(LocalTypeInfo &local) {
//   DeviceMatrix<Complex> &d_tmat_store=*get_dev_tmat_store();
//   d_tmat_store.copy_async(local.tmatStore,0);
// }
