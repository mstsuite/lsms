//
// Created by F.Moitzi on 15.12.2021.
//

#ifndef MADELUNG_SPHERICAL_HARMONICS_H
#define MADELUNG_SPHERICAL_HARMONICS_H

#include "Complex.hpp"

extern "C" {

void sph_harm_1_(double *vec, int *lmax, Complex *ylm);

}

#endif  // MADELUNG_SPHERICAL_HARMONICS_H