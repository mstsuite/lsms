/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */
#ifndef LSMS_ENERGYCONTOURINTEGRATION_H
#define LSMS_ENERGYCONTOURINTEGRATION_H

#include "Complex.hpp"
#include "SystemParameters.hpp"
#include "Communication/LSMSCommunication.hpp"

// typedef enum {EnergyGridBox=1, EnergyGridGauss=2} EnergyGridType;

void energyContourIntegration(LSMSCommunication &comm,LSMSSystemParameters &lsms, LocalTypeInfo &local, bool storeLocalGreenInt = false);

extern "C"
{
  void congauss_(double *ebot,double *etop,double *eibot,Complex *egrd,Complex *dele1,int *npts,int *nume,
                 double *pi,int * ipepts, int *iprint,char *istop, int istop_len);
}

inline void copyFromGreenIntLLp(Array3d<Complex> &to, int toEnergy, Array3d<Complex> &from, int kkrsz, int n_spin_cant, int n_spin_pola)
{
  if(n_spin_cant == 2)
  {
    for(int i=0; i<kkrsz; i++)
      for(int j=0; j<kkrsz; j++)
      {
	to(i,         j,         toEnergy) = from(i, j, 0);
	to(i + kkrsz, j,         toEnergy) = from(i, j, 1);
	to(i,         j + kkrsz, toEnergy) = from(i, j, 2);
	to(i + kkrsz, j + kkrsz, toEnergy) = from(i, j, 3);
      }
  } else if(n_spin_pola == 2) {
    for(int i=0; i<kkrsz; i++)
      for(int j=0; j<kkrsz; j++)
      {
	to(i,         j,         toEnergy) = from(i, j, 0);
	to(i + kkrsz, j,         toEnergy) = 0.0; // from(i, j, 2);
	to(i,         j + kkrsz, toEnergy) = 0.0; // from(i, j, 3);
	to(i + kkrsz, j + kkrsz, toEnergy) = from(i, j, 1);
      }
  } else {
    for(int i=0; i<kkrsz; i++)
      for(int j=0; j<kkrsz; j++)
      {
	to(i,         j,         toEnergy) = from(i, j, 0);
	to(i + kkrsz, j,         toEnergy) = 0.0; // from(i, j, 2);
	to(i,         j + kkrsz, toEnergy) = 0.0; // from(i, j, 3);
	to(i + kkrsz, j + kkrsz, toEnergy) = from(i, j, 0);
      }
  }
}

inline void copyToGreenIntLLp()
{
}

#endif
