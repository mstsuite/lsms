/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */
#ifndef LSMS_MIXING_H
#define LSMS_MIXING_H
#include "Real.hpp"
#include "SystemParameters.hpp"
#include "SingleSite/AtomData.hpp"
// #include "Communication/LSMSCommunication.hpp"

#include <vector>
#include <deque>
#include <cmath>

#include "LAPACK.hpp"

struct MixingParameters {

  // Different mixing quantities and algorithms
  static const int numQuantities = 5;

  enum mixQuantity {no_mixing = 0, charge = 1, potential = 2, moment_magnitude = 3,
                    moment_direction = 4};
  enum mixAlgorithm {noAlgorithm = 0, simple = 1, broyden = 2};
  
  // These parameters specify the which quantity(ies) is (are) being mixed and which algorithm(s) to used.
  // The correspondances of the indices are specified in mixQuantity.
  // bool values:
  // 0 : quantity is not used for mixing
  // 1 : quantity is used for mixing
  bool quantity[numQuantities];
  mixAlgorithm algorithm[numQuantities];
  Real mixingParameter[numQuantities];

};

#include "Communication/LSMSCommunication.hpp"

template <typename T>
void simpleMixing(T *fold, T* fnew, int n, Real alpha)
{
  if(alpha>1.0) alpha = 1.0;
  if(alpha<0.0) alpha = 0.0;
  Real beta = 1.0 - alpha;

  for(int i=0; i<n; i++)
    fold[i] = alpha * fnew[i] + beta * fold[i];
}

class MomentMixing {
public:
  virtual void update(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as) = 0;
};

class Mixing {
public:
  MomentMixing *momentMixing =  nullptr;

  virtual ~Mixing() = 0;
  // virtual void updateChargeDensity(LSMSSystemParameters &lsms, AtomData &a) = 0;
  virtual void updateChargeDensity(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as) = 0;
  // virtual void updatePotential(LSMSSystemParameters &lsms, AtomData &a) = 0;
  virtual void updatePotential(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as) = 0;
  void updateMoments(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as)
  { if(momentMixing != nullptr) momentMixing->update(comm, lsms, as);
    else {
      for(int i=0; i<as.size(); i++)
      {
        Real evecMagnitude = std::sqrt(as[i].evec[0] * as[i].evec[0] +
                                       as[i].evec[1] * as[i].evec[1] +
                                       as[i].evec[2] * as[i].evec[2]);
        as[i].evecNew[0] = as[i].evec[0] / evecMagnitude;
        as[i].evecNew[1] = as[i].evec[1] / evecMagnitude;
        as[i].evecNew[2] = as[i].evec[2] / evecMagnitude;
      }
    }
  }
  virtual void prepare(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as) = 0;
};


void setupMixing(MixingParameters &mix, Mixing* &mixing, int iprint);


#endif
