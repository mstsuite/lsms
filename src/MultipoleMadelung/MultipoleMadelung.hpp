//
// Created by F.Moitzi on 16.12.2021.
//

#ifndef MADELUNG_MULTIPOLEMADELUNG_HPP
#define MADELUNG_MULTIPOLEMADELUNG_HPP

#include <complex>
#include <vector>

#include "common.hpp"

namespace lsms {

class MultipoleMadelung {
 private:
  /// Global number of atoms
  int num_atoms;

  /// Local number of atoms
  int local_num_atoms;

  /// Madelung Matrix values
  matrix<double> madsum;
  array3d<std::complex<double>> dl_matrix;
  matrix<double> dl_factor;

  double scaling_factor;
  double rscut;
  double kncut;

  std::vector<int> r_nm;
  int nrslat;

  std::vector<int> k_nm;
  int nknlat;

 public:
  /// Angular-momentum index cutoff l
  int lmax;

  MultipoleMadelung(const matrix<double> &lattice,
                    const matrix<double> &atom_position, int num_atoms,
                    int lmax, std::vector<int> global_position_index);

  double getMadSum(int i, int j) const;

  std::complex<double> getDlMatrix(int g_atom, int k, int l_atom);

  double getDlFactor(int i, int j) const;

  double getScalingFactor() const;

  double getRsCut() const;

  double getKnCut() const;

  std::vector<int> getRsSize() const;

  std::vector<int> getKnSize() const;
};

}  // namespace lsms

#endif  // MADELUNG_MULTIPOLEMADELUNG_HPP
