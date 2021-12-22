
#ifndef MUST_LOCALTOTALENERGY_HPP
#define MUST_LOCALTOTALENERGY_HPP

#include "Main/SystemParameters.hpp"

extern "C"
{
void zeropt_(Real *ezpt, Real *tpzpt, Real *atvol, Real *ztotss);
}

void
localTotalEnergy(LSMSSystemParameters &lsms, AtomData &atom, Real &energy, Real &pressure);



#endif //MUST_LOCALTOTALENERGY_HPP
