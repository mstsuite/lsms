/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */

#ifndef LSMS_LINEAR_SOLVERS_HPP
#define LSMS_LINEAR_SOLVERS_HPP

#include "Complex.hpp"
#include "Matrix.hpp"
#include <vector>
#include <string>
#include <utility>

#include "MultipleScattering.hpp"

#ifdef ACCELERATOR_CUDA_C
#include "Accelerator/DeviceStorage.hpp"
#endif

#define MST_LINEAR_SOLVER_MASK    0x00000fff
#define MST_BUILD_KKR_MATRIX_MASK 0x0000f000

#ifndef MST_BUILD_KKR_MATRIX_DEFAULT
#define MST_BUILD_KKR_MATRIX_DEFAULT 0x1000
#endif

#ifndef MST_LINEAR_SOLVER_DEFAULT
#define MST_LINEAR_SOLVER_DEFAULT 0x800
#endif

#define MST_BUILD_KKR_MATRIX_F77         0x1000
#define MST_BUILD_KKR_MATRIX_CPP         0x2000
#define MST_BUILD_KKR_MATRIX_ACCELERATOR 0x3000

#define MST_LINEAR_SOLVER_ZGESV 1
void solveTau00zgesv(LSMSSystemParameters &lsms, LocalTypeInfo &local, AtomData &atom, int iie, Matrix<Complex> &m, Matrix<Complex> &tau00);
#define MST_LINEAR_SOLVER_ZGETRF 2
void solveTau00zgetrf(LSMSSystemParameters &lsms, LocalTypeInfo &local, AtomData &atom, int iie, Matrix<Complex> &m, Matrix<Complex> &tau00);
#define MST_LINEAR_SOLVER_ZCGESV 3
void solveTau00zcgesv(LSMSSystemParameters &lsms, LocalTypeInfo &local, AtomData &atom, int iie, Matrix<Complex> &m, Matrix<Complex> &tau00);
#define MST_LINEAR_SOLVER_ZBLOCKLU_F77 4
void solveTau00zblocklu_f77(LSMSSystemParameters &lsms, LocalTypeInfo &local, AtomData &atom, int iie, Matrix<Complex> &m, Matrix<Complex> &tau00);
#define MST_LINEAR_SOLVER_ZBLOCKLU_CPP 5
void solveTau00zblocklu_cpp(LSMSSystemParameters &lsms, LocalTypeInfo &local, AtomData &atom, int iie, Matrix<Complex> &m, Matrix<Complex> &tau00);

// #ifdef ACCELERATOR_CUBLAS
#define MST_LINEAR_SOLVER_ZGETRF_CUBLAS 0x10
#define MST_LINEAR_SOLVER_ZBLOCKLU_CUBLAS 0x11
// #endif

// #ifdef ACCELERATOR_CUSOLVER
#define MST_LINEAR_SOLVER_ZZGESV_CUSOLVER 0x12
#define MST_LINEAR_SOLVER_ZGETRF_CUSOLVER 0x13
// #endif
#ifdef ACCELERATOR_CUDA_C
void transferMatrixToGPUCuda(Complex *devM, Matrix<Complex> &m);
void transferMatrixFromGPUCuda(Matrix<Complex> &m, cuDoubleComplex *devM);
void transferT0MatrixToGPUCuda(Complex *devT0, LSMSSystemParameters &lsms, LocalTypeInfo &local, AtomData &atom, int iie);

void solveTau00zgetrf_cublas(LSMSSystemParameters &lsms, LocalTypeInfo &local, DeviceStorage &d, AtomData &atom, Complex *tMatrix, Complex *devM, Matrix<Complex> &tau00);
void solveTau00zzgesv_cusolver(LSMSSystemParameters &lsms, LocalTypeInfo &local, DeviceStorage &d, AtomData &atom, Complex *tMatrix, Complex *devM, Matrix<Complex> &tau00);
void solveTau00zgetrf_cusolver(LSMSSystemParameters &lsms, LocalTypeInfo &local, DeviceStorage &d, AtomData &atom, Complex *tMatrix, Complex *devM, Matrix<Complex> &tau00);
#endif

#ifdef ACCELERATOR_HIP
#endif

#define MST_LINEAR_SOLVER_BLOCK_INVERSE_F77 0xf00
#define MST_LINEAR_SOLVER_BLOCK_INVERSE_CPP 0xf01
// #ifdef ACCELERATOR_CUDA_C
#define MST_LINEAR_SOLVER_BLOCK_INVERSE_CUDA 0xf10
// #endif

inline std::string linearSolverName(unsigned int solverId)
{
  solverId = solverId & 0xffff;
  std::string name("");
  char idstr[12];
  if(solverId == 0)
    {
      solverId = MST_LINEAR_SOLVER_DEFAULT;
      name = "default solver: ";
    }
  snprintf(idstr, 10, " (0x%04x)", solverId);
  switch(solverId)
    {
    case MST_LINEAR_SOLVER_ZGESV: name += "CPU zgesv"; break;
    case MST_LINEAR_SOLVER_ZGETRF: name += "CPU zgetrf"; break;
    case MST_LINEAR_SOLVER_ZCGESV: name += "CPU zcgesv"; break;
    case MST_LINEAR_SOLVER_ZBLOCKLU_F77: name += "CPU zblocklu f77"; break;
    case MST_LINEAR_SOLVER_ZBLOCKLU_CPP: name += "CPU zblocklu c++"; break;
    case MST_LINEAR_SOLVER_ZGETRF_CUBLAS: name += "CUBLAS zgetrf"; break;
    case MST_LINEAR_SOLVER_ZBLOCKLU_CUBLAS: name += "CUBLAS zblocklu"; break;
    case MST_LINEAR_SOLVER_ZZGESV_CUSOLVER: name += "CUSOLVER zzgesv"; break;
    case MST_LINEAR_SOLVER_ZGETRF_CUSOLVER: name += "CUSOLVER zgetrf"; break;
    case MST_LINEAR_SOLVER_BLOCK_INVERSE_F77: name += "block inverse f77"; break;
    case MST_LINEAR_SOLVER_BLOCK_INVERSE_CPP: name += "block inverse c++"; break;
    case MST_LINEAR_SOLVER_BLOCK_INVERSE_CUDA: name += "block inverse cuda"; break;
    default: name += "unknwon solver";
    }
  name += idstr;
  return name;
}

#endif
