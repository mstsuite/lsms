/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */
#ifndef LSMSCOMMUNICATION_H
#define LSMSCOMMUNICATION_H

#include <mpi.h>

#include <cstdio>
#include <cstdlib>
#include <vector>

#include "Main/SystemParameters.hpp"
#include "Mixer/MixingParameter.hpp"
#include "Potential/PotentialShifter.hpp"
#include "SingleSite/AtomData.hpp"
#include "Main/mixing_params.hpp"

class TmatCommType {
 public:
  int remoteNode;
  int numTmats;
  std::vector<int> tmatStoreIdx;
  std::vector<int> globalIdx;
  std::vector<MPI_Request> communicationRequest;
};

class LSMSCommunication {
 public:
  int rank;
  int size;
  MPI_Comm comm;

  int numTmatTo, numTmatFrom;
  std::vector<TmatCommType> tmatTo, tmatFrom;
};

void initializeCommunication(LSMSCommunication &comm);

void initializeCommunication(LSMSCommunication &comm, MPI_Comm mpiCommunicator);

void finalizeCommunication();

void exitLSMS(LSMSCommunication &comm, int errorCode);
void synchronizeLSMS(LSMSCommunication &comm);

void communicateParameters(LSMSCommunication &comm, LSMSSystemParameters &lsms,
                           CrystalParameters &crystal,
                           lsms::MixingParameterPack &mix,
                           AlloyMixingDesc &alloyDesc);

void communicateParameters(LSMSCommunication &comm, LSMSSystemParameters &lsms,
                           CrystalParameters &crystal, MixingParameters &mix, AlloyMixingDesc &alloyDesc);


void communicateSingleAtomData(LSMSCommunication &comm, int from, int to,
                               int &local_id, AtomData &atom, int tag = 0);

void communicatePotentialShiftParameters(LSMSCommunication &comm, PotentialShifter &ps);

void expectTmatCommunication(LSMSCommunication &comm, LocalTypeInfo &local);

void sendTmats(LSMSCommunication &comm, LocalTypeInfo &local);

void finalizeTmatCommunication(LSMSCommunication &comm);

void printCommunicationInfo(FILE *f, LSMSCommunication &comm);

template<typename T>
void globalMax(LSMSCommunication &comm, T &a) {
  T r;
  MPI_Allreduce(&a, &r, 1, TypeTraits<T>::mpiType(), MPI_MAX, comm.comm);
  a = r;
}

void globalAnd(LSMSCommunication &comm, bool &a);

template <typename T>
void globalSum(LSMSCommunication &comm, T &a) {
  T r;
  MPI_Allreduce(&a, &r, 1, TypeTraits<T>::mpiType(), MPI_SUM, comm.comm);
  a = r;
}

template <typename T>
void globalSum(LSMSCommunication &comm, T *a, int n) {
  std::vector<T> r(n);
  MPI_Allreduce(a, r.data(), n, TypeTraits<T>::mpiType(), MPI_SUM, comm.comm);
  std::copy(r.begin(), r.end(), a);
}

double calculateFomScaleDouble(LSMSCommunication &comm, LocalTypeInfo &local);

long long calculateFomScale(LSMSCommunication &comm, LocalTypeInfo &local);

#endif
