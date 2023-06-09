#ifndef LSMS_MADELUNG_H
#define LSMS_MADELUNG_H

#include "Array3d.hpp"
#include "Complex.hpp"
#include "Main/SystemParameters.hpp"
#include "Matrix.hpp"
#include "Real.hpp"

void calculateMadelungMatrices(LSMSSystemParameters &lsms,
                               CrystalParameters &crystal,
                               LocalTypeInfo &local);

extern "C" {
void cal_madelung_matrix_(int *mynod, int *num_atoms, Real *bravais_lattice_in,
                          Real *atom_posi_x_in, Real *atom_posi_y_in,
                          Real *atom_posi_z_in, Real *madmat, int *iprint,
                          char *istop, int istop_len);

void cal_madelung_matrix_j_(int *mynod, int *num_atoms,
                            Real *bravais_lattice_in, Real *rcutau, int *lmax,
                            int *ndlm, int *ndlmv, Real *atom_posi_x_in,
                            Real *atom_posi_y_in, Real *atom_posi_z_in,
                            Real *madmat, Complex *madmatj, int *iprint,
                            char *istop, int istop_len);
}

#endif
