/* -*- mode: C++; c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */

#ifndef INVERSION_TEST_CUDA
#define INVERSION_TEST_CUDA

#include <stdio.h>

#include "Complex.hpp"
#include "Matrix.hpp"
#include <vector>
#include <chrono>
#include <ctime>

#include <cuda_runtime.h>
#include <cuComplex.h>
#include <cublas_v2.h>

class DeviceData
{
public:
  std::vector<cuDoubleComplex *> tMatrices;
  cuDoubleComplex *tMatrixStore;
  cuDoubleComplex *tau;
  cuDoubleComplex *tau00;
  cuDoubleComplex *m;
  cuDoubleComplex *G0;
  int *ipiv;
  int *info;
};

void allocDeviceData(DeviceData &d, int blockSize, int numBlocks);
void freeDeviceData(DeviceData &d);

void usage_cuda(const char *name);

void makeType1MatrixGPU(cublasHandle_t cublasHandle, DeviceData &d, int blockSize, int numBlocks);

void transferMatrixToGPU(cuDoubleComplex *devM, Matrix<Complex> &m);
void transferMatrixFromGPU(Matrix<Complex> &m, cuDoubleComplex *devM);

void solveTau00zgetrf_cublas(cublasHandle_t cublasHandle, DeviceData &d,
                             Matrix<Complex> &tau00, int blockSize, int numBlocks);

void initCuda(cublasHandle_t &cublasHandle);
void finalizeCuda(cublasHandle_t &cublasHandle);

#endif
