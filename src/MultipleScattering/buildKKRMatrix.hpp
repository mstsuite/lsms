/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */

#ifndef LSMS_BUILD_KKR_MATRIX_HPP
#define LSMS_BUILD_KKR_MATRIX_HPP

#include "Complex.hpp"
#include "Matrix.hpp"
#include <vector>
#include <string>
#include <utility>

#include "MultipleScattering.hpp"
#include "linearSolvers.hpp"

#ifdef ACCELERATOR_CUDA_C
#include "Accelerator/DeviceStorage.hpp"
#endif

#ifndef MST_BUILD_KKR_MATRIX_DEFAULT
#define MST_BUILD_KKR_MATRIX_DEFAULT 0x1000
#endif

#define MST_BUILD_KKR_MATRIX_F77         0x1000
#define MST_BUILD_KKR_MATRIX_CPP         0x2000
void buildKKRMatrixCPU(LSMSSystemParameters &lsms, LocalTypeInfo &local, AtomData &atom, int iie, Complex energy, Complex prel,
                    Matrix<Complex> &m);
#define MST_BUILD_KKR_MATRIX_ACCELERATOR 0x3000
#ifdef ACCELERATOR_CUDA_C
#endif

#ifdef ACCELERATOR_HIP
#endif

inline std::string buildKKRMatrixName(unsigned int buildKKRMatrixId)
{
  buildKKRMatrixId = buildKKRMatrixId & MST_BUILD_KKR_MATRIX_MASK;
  std::string name("");
  char idstr[12];
  if(buildKKRMatrixId == 0)
    {
      buildKKRMatrixId = MST_BUILD_KKR_MATRIX_DEFAULT;
      name = "default buildKKRMatrix: ";
    }
  snprintf(idstr, 10, " (0x%04x)", buildKKRMatrixId);
  switch(buildKKRMatrixId)
    {
    case MST_BUILD_KKR_MATRIX_F77: name += "LSMS 1 buildKKRMatrix"; break;
    case MST_BUILD_KKR_MATRIX_CPP: name += "CPU buildKKRMatrix"; break;
    case MST_BUILD_KKR_MATRIX_ACCELERATOR: name += "Accelerator buildKKRMatrix"; break;
    default: name += "unknwon buildKKRMatrix";
    }
  name += idstr;
  return name;
}

#endif
