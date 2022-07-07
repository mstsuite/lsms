//
// Created by F.Moitzi on 21.12.2021.
//

#ifndef MADELUNG_LATTICE_HPP
#define MADELUNG_LATTICE_HPP

#include <complex>
#include <tuple>
#include <vector>

#include "common.hpp"

namespace lsms {

/**
 * Create all vectors in certain cutoff with a certain repeotition
 */
matrix<double> create_lattice(const matrix<double> &brav, double cutoff,
                              const std::vector<int> &nm, int size);

/**
 * Create all vectors in certain cutoff with a certain repeotition
 */
std::tuple<matrix<double>, std::vector<double>> create_lattice_and_sq(
    matrix<double> &brav, double cutoff, const std::vector<int> &nm, int size);

/**
 *  inserts a vector in a list of vectors such that they are in
 *  order of increasing length.
 */
void insert_ordered(matrix<double> &latt_vec, std::vector<double> &latt_vec_sq,
                    int len, std::vector<double> &vec, double &v_sq);

/**
 * Lattice volumes
 */
double omega(matrix<double> &bravais);

/**
 * Calculate reciprocal lattice
 */
void reciprocal_lattice(matrix<double> &bravais,
                        matrix<double> &reciprocal_bravais, double &scale);

}  // namespace lsms

#endif  // MADELUNG_LATTICE_HPP
