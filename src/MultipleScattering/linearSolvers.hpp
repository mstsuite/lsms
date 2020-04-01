/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */

#ifndef LSMS_LINEAR_SOLVERS_HPP
#define LSMS_LINEAR_SOLVERS_HPP

#include "Complex.hpp"
#include "Matrix.hpp"
#include <vector>

#include "MultipleScattering.hpp"

#define MST_LINEAR_SOLVER_DEFAULT 10000

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

#ifdef ACCELERATOR_CUBLAS
#define MST_LINEAR_SOLVER_ZGETRF_CUBLAS 100
#define MST_LINEAR_SOLVER_ZBLOCKLU_CUBLAS 101
#endif

#ifdef ACCELERATOR_CUSOLVER
#define MST_LINEAR_SOLVER_ZZGESV_CUSOLVER 110
#define MST_LINEAR_SOLVER_ZGETRF_CUSOLVER 111
#endif

#ifdef ACCELERATOR_HIP
#endif

#define MST_LINEAR_SOLVER_BLOCK_INVERSE_F77 10000
#define MST_LINEAR_SOLVER_BLOCK_INVERSE_CPP 10001
#ifdef ACCELERATOR_CUDA_C
#define MST_LINEAR_SOLVER_BLOCK_INVERSE_CUDA 10002
#endif


#endif
