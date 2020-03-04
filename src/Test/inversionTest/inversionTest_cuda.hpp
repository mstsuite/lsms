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
#include <cusolverDn.h>

class DeviceData
{
public:
  std::vector<cuDoubleComplex *> tMatrices;
  cuDoubleComplex *tMatrixStore;
  cuDoubleComplex *tau;
  cuDoubleComplex *t;
  cuDoubleComplex *tau00;
  cuDoubleComplex *m;
  cuDoubleComplex *G0;
  int *ipiv;
  int *info;
  size_t workBytes;
  void *work;
};

class DeviceHandles
{
public:
  cublasHandle_t cublasHandle;
  cusolverDnHandle_t cusolverDnHandle;
};


void allocDeviceData(DeviceHandles &deviceHandles, DeviceData &d, int blockSize, int numBlocks);
void freeDeviceData(DeviceData &d);

void usage_cuda(const char *name);

void makeType1MatrixGPU(cublasHandle_t cublasHandle, DeviceData &d, int blockSize, int numBlocks);

void transferMatrixToGPU(cuDoubleComplex *devM, Matrix<Complex> &m);
void transferMatrixFromGPU(Matrix<Complex> &m, cuDoubleComplex *devM);

void solveTau00zgetrf_cublas(cublasHandle_t cublasHandle, DeviceData &d,
                             Matrix<Complex> &tau00, int blockSize, int numBlocks);
void solveTau00zblocklu_cublas(cublasHandle_t cublasHandle, DeviceData &d, Matrix<Complex> &tau00, Matrix<Complex> &m, std::vector<Matrix<Complex> > &tMatrices, int blockSize, int numBlocks);
void solveTau00zzgesv_cusolver(DeviceHandles &deviceHandles, DeviceData &deviceData, Matrix<Complex>
    &tau00, Matrix<Complex> &m, std::vector<Matrix<Complex> > &tMatrices, int blockSize, int
    numBlocks);

void solveTau00zgetrf_cusolver(DeviceHandles &deviceHandles, DeviceData &deviceData,
    Matrix<Complex> &tau00, int blockSize, int numBlocks);

void transferTest(DeviceHandles &deviceHandles, DeviceData &deviceData, Matrix<Complex> &tau00, int
    blockSize, int numBlocks);

void initCuda(DeviceHandles &deviceHandles);
void finalizeCuda(DeviceHandles &deviceHandles);

#endif
