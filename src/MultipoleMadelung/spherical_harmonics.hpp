//
// Created by F.Moitzi on 15.12.2021.
//

#ifndef MADELUNG_SPHERICAL_HARMONICS_H
#define MADELUNG_SPHERICAL_HARMONICS_H

#include <complex>

extern "C"
{

void sph_harm_1(double *vec, int *lmax, std::complex<double> *ylm);

}


#endif //MADELUNG_SPHERICAL_HARMONICS_H