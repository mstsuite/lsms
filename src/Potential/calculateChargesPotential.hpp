/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */
#ifndef CALCULATE_CHARGES_HPP
#define CALCULATE_CHARGES_HPP

#include "Array3d.hpp"
#include "Communication/LSMSCommunication.hpp"
#include "Main/SystemParameters.hpp"
#include "Real.hpp"
#include "getvmt.hpp"

extern "C" {
void getqm_mt_(int *n_spin_pola, int *ilast, Real *rsmo, Real *sqr, Real *rho,
               int *iprpts, Real *rhot, int *i_smooth_rsmo, Real *qtotmt,
               Real *mtotmt, Real *r_sph, int *iprint);

void newexchg_(int *n_spin_pola, Real *sp, Real *rhoup, Real *rhodn, Real *vx,
               Real *enxc, Real *vxout, Real *excout, Real *ro3, Real *dz,
               Real *r_mesh, int *jmt, int *iexch);

void newpot_(int *n_spin_pola, Real *ztotss, Real *rhoup, Real *rhodn,
             Real *rhot, Real *vrold, Real *vrnew, Real *vrms, Real *vx,
             Real *vmt1, Real *vmt, Real *vxout, Real *sqr, int *jmt,
             Real *rins, Real *rmt, int *mtasa);  //, int *iexch);

void calculate_asa_ro3_(int *n_spin_pola, Real *rho1, Real *rho2, int *i_vdif,
                        Real *r_mesh, int *jmt, Real *rmt, int *iprpts,
                        Real *ro3, Real *dz);
}

/* ============================================================
   chargeSwitch:

   = 0 (default)
   * before mixing, everything calculated were put in xxxxNew
   * qtotws and mtotws calculated from xvalssNew
   * newexchg and newpot called with rhoNew
   * calculateMTZeroPotDiff calculates vdifNew

   = 1 :
   * after mixing, mixed quantities were put into xxxx (the old one)
   * qtotws and mtotws calculated from xvalss
   * newexchg and newpot called with rhotot
   * calculateMTZeroPotDiff calculates vdif
*/

void calculateChargesPotential(LSMSCommunication &comm,
                               LSMSSystemParameters &lsms, LocalTypeInfo &local,
                               CrystalParameters &crystal,
                               int chargeSwitch = 0);

void calculateCharges(LSMSCommunication &comm, LSMSSystemParameters &lsms,
                      LocalTypeInfo &local, CrystalParameters &crystal,
                      std::vector<Real> &qsub, Array3d<Real> &rhoTemp,
                      int chargeSwitch = 0);

void calculateLocalCharges(LSMSSystemParameters &lsms, LocalTypeInfo &local,
                           int chargeSwitch = 0);

void calculatePotential(LSMSCommunication &comm, LSMSSystemParameters &lsms,
                        LocalTypeInfo &local, CrystalParameters &crystal,
                        std::vector<Real> &qsub, Array3d<Real> &rhoTemp,
                        int chargeSwitch = 0);

void calculateMTZeroPotDiff(LSMSSystemParameters &lsms, LocalTypeInfo &local,
                            int chargeSwitch = 0);

#endif
