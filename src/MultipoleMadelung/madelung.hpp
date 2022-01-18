//
// Created by F.Moitzi on 06.01.2022.
//

#ifndef SRC_MADELUNG_MADELUNG_HPP
#define SRC_MADELUNG_MADELUNG_HPP

#include <complex>
#include <vector>

#include "common.hpp"

namespace lsms {

constexpr auto EPSI = 1e-14;

/**
 * Calculate the scaling factor to get a balanced number
 * of reciprocal and real-space lattice vectors
 */
double scaling_factor(const matrix<double> &bravais, int lmax,
                      int max_iter = 100, double fstep = 0.02);

/**
 * Number of lattice vectors
 */
int num_latt_vectors(const matrix<double> &brav, double cut,
                     const std::vector<int> &nm);

/**
 * Get radius of truncation sphere
 */
double rs_trunc_radius(const matrix<double> &brav, int lmax, double eta,
                       const std::vector<int> &nm);

double kn_trunc_radius(const matrix<double> &brav, int lmax, double eta,
                       const std::vector<int> &nm);

/**
 * Get size of lattice multiplications
 */
std::vector<int> real_space_multiplication(const matrix<double> &brav, int lmax,
                                           double eta);

/**
 * Get size of reciprocal lattice multiplications
 */
std::vector<int> reciprocal_space_multiplication(const matrix<double> &brav,
                                                 int lmax, double eta);

/**
 *
 */
void calculate_madelung_matrix(
    int num_atoms, int myatom, int id,

    int lmax_mad, double eta, double a0,

    matrix<double> &r_brav, const matrix<double> &atom_position,

    matrix<double> &rslat, std::vector<double> &rslatsq,

    matrix<double> &knlat, std::vector<double> &knlatsq,

    matrix<double> &madmat, array3d<std::complex<double>> &dl_matrix);

double calculate_eta(matrix<double> &brav);

matrix<double> calculate_dl_factor(int lmax);

}  // namespace lsms

#endif  // SRC_MADELUNG_MADELUNG_HPP
