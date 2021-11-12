/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */
#ifndef LSMS_NEWFUNCTIONALINTERFACE_HPP
#define LSMS_NEWFUNCTIONALINTERFACE_HPP

#include "Main/SystemParameters.hpp"
#include "Real.hpp"

enum class XCFunctional
{
  VoskoWilkNusair
};

// LibxcInterface is a singleton class to provide an interface to libxc
class NewFunctionalInterface {
public:
  XCFunctional functional[numFunctionalIndices-1];
  int numFunctionals;
  bool spinPolarized;
  bool needGradients; // the functional needs gradients of the density (for GGAs)
  bool needLaplacian; // need laplacians of the density (for MetaGGAs)
  bool needKineticEnergyDensity; // for MetaGGAs
  bool needExactExchange; // for Hybrid Functionals

  int init(int nSpin, int *xcFunctional);
  void evaluate(std::vector<Real> &rMesh, Matrix<Real> &rhoIn, int jmt, int nSpin, Matrix<Real> &xcEnergyOut, Matrix<Real> &xcPotOut);
  void evaluateSingle(Real *rhoIn, int nSpin, Real *xcEnergyOut, Real *xcPotOut);
};

#endif
