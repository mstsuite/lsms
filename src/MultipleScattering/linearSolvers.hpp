/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */

#ifndef LSMS_LINEAR_SOLVERS_HPP
#define LSMS_LINEAR_SOLVERS_HPP

#include <string>
#include <utility>
#include <vector>

#include "Complex.hpp"
#include "Matrix.hpp"
#include "MultipleScattering.hpp"

#if defined(ACCELERATOR_CUDA_C) || defined(ACCELERATOR_HIP)
#include "Accelerator/DeviceStorage.hpp"
#endif

#define MST_LINEAR_SOLVER_MASK 0x00000fff
#define MST_BUILD_KKR_MATRIX_MASK 0x0000f000

#ifndef MST_BUILD_KKR_MATRIX_DEFAULT
#define MST_BUILD_KKR_MATRIX_DEFAULT 0x1000
#endif

#ifndef MST_LINEAR_SOLVER_DEFAULT
#define MST_LINEAR_SOLVER_DEFAULT 0xf00
#endif

#define MST_BUILD_KKR_MATRIX_F77 0x1000
#define MST_BUILD_KKR_MATRIX_CPP 0x2000
#define MST_BUILD_KKR_MATRIX_ACCELERATOR 0x3000

#define MST_LINEAR_SOLVER_ZGESV 1
void solveTau00zgesv(LSMSSystemParameters &lsms, LocalTypeInfo &local,
                     AtomData &atom, int iie, Matrix<Complex> &m,
                     Matrix<Complex> &tau00, int ispin);
#define MST_LINEAR_SOLVER_ZGETRF 2
void solveTau00zgetrf(LSMSSystemParameters &lsms, LocalTypeInfo &local,
                      AtomData &atom, int iie, Matrix<Complex> &m, Matrix<Complex> &tau00, int ispin);
void solveTauFullzgetrf(LSMSSystemParameters &lsms, LocalTypeInfo &local, 
                      AtomData &atom, Matrix<Complex> &m, Matrix<Complex> &tau, int ispin);

#define MST_LINEAR_SOLVER_ZCGESV 3
void solveTau00zcgesv(LSMSSystemParameters &lsms, LocalTypeInfo &local,
                      AtomData &atom, int iie, Matrix<Complex> &m,
                      Matrix<Complex> &tau00, int ispin);
#define MST_LINEAR_SOLVER_ZBLOCKLU_F77 4
void solveTau00zblocklu_f77(LSMSSystemParameters &lsms, LocalTypeInfo &local,
                            AtomData &atom, int iie, Matrix<Complex> &m,
                            Matrix<Complex> &tau00, int ispin);
#define MST_LINEAR_SOLVER_ZBLOCKLU_CPP 5
void solveTau00zblocklu_cpp(LSMSSystemParameters &lsms, LocalTypeInfo &local,
                            AtomData &atom, int iie, Matrix<Complex> &m,
                            Matrix<Complex> &tau00, int ispin);

// #ifdef ACCELERATOR_CUBLAS
#define MST_LINEAR_SOLVER_ZGETRF_CUBLAS 0x10
#define MST_LINEAR_SOLVER_ZBLOCKLU_CUBLAS 0x11
// #endif

// #ifdef ACCELERATOR_CUSOLVER
#define MST_LINEAR_SOLVER_ZZGESV_CUSOLVER 0x12
#define MST_LINEAR_SOLVER_ZGETRF_CUSOLVER 0x13
#define MST_LINEAR_SOLVER_XGETRF_CUSOLVER 0x14
#define MST_LINEAR_SOLVER_IRSXGESV_CUSOLVER 0x15
// #endif

// #ifdef ACCELERATOR_HIP
#define MST_LINEAR_SOLVER_ZGETRF_HIPBLAS 0x20
#define MST_LINEAR_SOLVER_ZGETRF_ROCSOLVER 0x20
// #endif

#ifdef ACCELERATOR_CUDA_C
void transferMatrixToGPUCuda(Complex *devM, Matrix<Complex> &m);
void transferMatrixFromGPUCuda(Matrix<Complex> &m, cuDoubleComplex *devM);
void transferT0MatrixToGPUCuda(Complex *devT0, LSMSSystemParameters &lsms, LocalTypeInfo &local,
                               AtomData &atom, int iie, int ispin);
void transferFullTMatrixToGPUCUDA(Complex *devT, LSMSSystemParameters &lsms, LocalTypeInfo &local,
                                  AtomData &atom, int ispin);

void solveTau00zgetrf_cublas(LSMSSystemParameters &lsms, LocalTypeInfo &local, DeviceStorage &d, AtomData &atom, Complex *tMatrix, Complex *devM, Matrix<Complex> &tau00);
void solveTauFullzgetrf_cublas(LSMSSystemParameters &lsms, LocalTypeInfo &local, DeviceStorage &d, AtomData &atom,
                               Complex *tMatrix, Complex *devM, Complex *devTauFull);


// CUDA Solvers
void solveTau00zzgesv_cusolver(LSMSSystemParameters &lsms, LocalTypeInfo &local, DeviceStorage &d, AtomData &atom,
                               Complex *tMatrix, Complex *devM, Matrix<Complex> &tau00, int ispin);
void solveTau00zgetrf_cusolver(LSMSSystemParameters &lsms, LocalTypeInfo &local, DeviceStorage &d, AtomData &atom,
                               Complex *tMatrix, Complex *devM, Matrix<Complex> &tau00, int ispin);
void solveTauFullzgetrf_cusolver(LSMSSystemParameters &lsms, LocalTypeInfo &local, DeviceStorage &d, AtomData &atom,
                               Complex *tMatrix, Complex *devM, Complex *devTauFull, int ispin);

#ifdef USE_XGETRF
void solveTau00Xgetrf_cusolver(LSMSSystemParameters &lsms, LocalTypeInfo &local,
                               DeviceStorage &d, AtomData &atom,
                               Complex *tMatrix, Complex *devM,
                               Matrix<Complex> &tau00, int ispin);
#endif

#ifdef USE_IRSXGESV
void solveTau00IRSXgesv_cusolver(LSMSSystemParameters &lsms,
                                 LocalTypeInfo &local, DeviceStorage &d,
                                 AtomData &atom, Complex *tMatrix,
                                 Complex *devM, Matrix<Complex> &tau00,
                                 int ispin);
#endif

#define IDX(i, j, lDim) (((j) * (lDim)) + (i))

#ifdef __CUDACC__
template <typename T>
void zeroMatrixCuda(T *devM, int lDim, int nCol) {
  //  for(int i=0; i<m.n_row(); i++)
  //    for(int j=0; j<m.n_col(); j++)
  //      m(i,j) = 0.0;
  cudaMemset(devM, 0, lDim * nCol * sizeof(T));
}

template <typename T>
__global__ void setDiagonalKernelCuda(T *devM, int lDim, int nCol, T val) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < nCol) {
    devM[IDX(i, i, lDim)] = val;
  }
}

template <typename T>
__global__ void addDiagonalKernelCuda(T *devM, int lDim, int nCol, T val) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < nCol) {
    devM[IDX(i, i, lDim)] = cuCadd(devM[IDX(i, i, lDim)], val);
  }
}

template <typename T>
void unitMatrixCuda(T *devM, int lDim, int nCol) {
  zeroMatrixCuda(devM, lDim, nCol);
  setDiagonalKernelCuda<<<nCol, 1>>>(devM, lDim, nCol, T(1.0));
}
#endif

#endif

#ifdef ACCELERATOR_HIP
void transferMatrixToGPUHip(Complex *devM, Matrix<Complex> &m);
void transferMatrixFromGPUHip(Matrix<Complex> &m, hipDoubleComplex *devM);
void transferT0MatrixToGPUHip(Complex *devT0, LSMSSystemParameters &lsms, LocalTypeInfo &local, AtomData &atom, int iie);
void transferFullTMatrixToGPUHip(Complex *devT, LSMSSystemParameters &lsms, LocalTypeInfo &local,
                                  AtomData &atom, int ispin);

void solveTau00zgetrf_rocsolver(LSMSSystemParameters &lsms,
                                LocalTypeInfo &local, DeviceStorage &d,
                                AtomData &atom, Complex *tMatrix, Complex *devM,
                                Matrix<Complex> &tau00);

void solveTauFullzgetrf_rocsolver(LSMSSystemParameters &lsms, LocalTypeInfo &local, DeviceStorage &d, AtomData &atom,
                               Complex *tMatrix, Complex *devM, Complex *devTauFull, int ispin);

#define IDX(i, j, lDim) (((j)*(lDim))+(i))

#ifdef __HIPCC__
template <typename T>
void zeroMatrixHip(T *devM, int lDim, int nCol) {
  //  for(int i=0; i<m.n_row(); i++)
  //    for(int j=0; j<m.n_col(); j++)
  //      m(i,j) = 0.0;
  hipMemset(devM, 0, lDim * nCol * sizeof(T));
}

template <typename T>
__global__ void setDiagonalKernelHip(T *devM, int lDim, int nCol, T val) {
  int i = hipBlockIdx_x * hipBlockDim_x + hipThreadIdx_x;
  if (i < nCol) {
    devM[IDX(i, i, lDim)] = val;
  }
}

template <typename T>
__global__ void addDiagonalKernelHip(T *devM, int lDim, int nCol, T val) {
  int i = hipBlockIdx_x * hipBlockDim_x + hipThreadIdx_x;
  if (i < nCol) {
    devM[IDX(i, i, lDim)] = hipCadd(devM[IDX(i, i, lDim)], val);
  }
}

template <typename T>
void unitMatrixHip(T *devM, int lDim, int nCol) {
  zeroMatrixHip(devM, lDim, nCol);
  setDiagonalKernelHip<<<nCol, 1>>>(devM, lDim, nCol, T(1.0));
}
#endif

#endif

#define MST_LINEAR_SOLVER_BLOCK_INVERSE_F77 0xf00
#define MST_LINEAR_SOLVER_BLOCK_INVERSE_CPP 0xf01
// #ifdef ACCELERATOR_CUDA_C
#define MST_LINEAR_SOLVER_BLOCK_INVERSE_CUDA 0xf10
// #endif

inline std::string linearSolverName(unsigned int solverId) {
  solverId = solverId & MST_LINEAR_SOLVER_MASK;
  std::string name("");
  char idstr[12];
  if (solverId == 0) {
    solverId = MST_LINEAR_SOLVER_DEFAULT;
    name = "default solver: ";
  }
  snprintf(idstr, 10, " (0x%03x)", solverId);
  switch (solverId) {
    case MST_LINEAR_SOLVER_ZGESV:
      name += "CPU zgesv";
      break;
    case MST_LINEAR_SOLVER_ZGETRF:
      name += "CPU zgetrf";
      break;
    case MST_LINEAR_SOLVER_ZCGESV:
      name += "CPU zcgesv";
      break;
    case MST_LINEAR_SOLVER_ZBLOCKLU_F77:
      name += "CPU zblocklu f77";
      break;
    case MST_LINEAR_SOLVER_ZBLOCKLU_CPP:
      name += "CPU zblocklu c++";
      break;

    case MST_LINEAR_SOLVER_ZGETRF_CUBLAS:
      name += "CUBLAS zgetrf";
      break;
    case MST_LINEAR_SOLVER_ZBLOCKLU_CUBLAS:
      name += "CUBLAS zblocklu";
      break;
    case MST_LINEAR_SOLVER_ZZGESV_CUSOLVER:
      name += "CUSOLVER zzgesv";
      break;
    case MST_LINEAR_SOLVER_ZGETRF_CUSOLVER:
      name += "CUSOLVER zgetrf";
      break;

    case MST_LINEAR_SOLVER_ZGETRF_ROCSOLVER:
      name += "ROCSOLVER zgetrf";
      break;

    case MST_LINEAR_SOLVER_BLOCK_INVERSE_F77:
      name += "block inverse f77";
      break;
    case MST_LINEAR_SOLVER_BLOCK_INVERSE_CPP:
      name += "block inverse c++";
      break;
    case MST_LINEAR_SOLVER_BLOCK_INVERSE_CUDA:
      name += "block inverse cuda";
      break;

    default:
      name += "unknwon solver";
  }
  name += idstr;
  return name;
}

#endif
