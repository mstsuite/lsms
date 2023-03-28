
#ifndef MUST_LOCALTOTALENERGY_HPP
#define MUST_LOCALTOTALENERGY_HPP

#include "Main/SystemParameters.hpp"
#include "TotalEnergy/DFTEnergy.hpp"

extern "C" {
void zeropt_(Real *ezpt, Real *tpzpt, Real *atvol, Real *ztotss);
}

void localTotalEnergy(LSMSSystemParameters &lsms, AtomData &atom, Real &energy,
                      Real &pressure);

void localTotalEnergy(LSMSSystemParameters &lsms, AtomData &atom,
                      lsms::DFTEnergy &dft_energy);

#endif  // MUST_LOCALTOTALENERGY_HPP
