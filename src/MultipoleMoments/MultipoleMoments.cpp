//
// Created by F.Moitzi on 24.01.2022.
//

#include "MultipoleMoments.hpp"

#include <complex>
#include <vector>

#include "Array3d.hpp"
#include "NDArray.hpp"
#include "integer_factors.hpp"

std::vector<std::complex<double>> lsms::multi_moms::multipole_mom_e(
    const NonRelativisticSingleScattererSolution &solution,
    const Matrix<std::complex<double>> &tau00_l,
    const lsms::math::GauntFactorBase &gaunt_factor, unsigned int lmax_mm,
    unsigned int lmax) {
  int ispin = 0;

  // Local parameters
  auto &rmesh = solution.atom->r_mesh;
  auto &jmt = solution.atom->jmt;

  //
  auto kmax_mm = lsms::get_kmax(lmax_mm);
  auto kmax = lsms::get_kmax(lmax);
  std::vector<std::complex<double>> multi_mom_e(kmax_mm, 0.0);

  auto lofk = lsms::get_lofk(std::max(kmax, kmax_mm));
  auto mofk = lsms::get_mofk(std::max(kmax, kmax_mm));

  auto k = 0;
  for (auto l{0}; l <= lmax_mm; l++) {
    for (auto m = -l; m <= l; m++) {

      std::complex<double> value{0.0};

      for (auto k_p{0}; k_p < kmax; k_p++) {
        auto l_p = lofk[k_p];
        auto m_p = mofk[k_p];

        for (auto k_pp{0}; k_pp < kmax; k_pp++) {

          auto l_pp = lofk[k_pp];
          auto m_pp = mofk[k_pp];

          std::vector<std::complex<double>> factor(jmt);

          auto gaunt_coeff = gaunt_factor.table(l_p, m_p, l_pp, m_pp, l, m);

          if (gaunt_coeff != 0.0) {

            std::complex<double> ss_terms{0.0};

            if (k_p == k_pp) {

              for (auto idx = 0; idx < jmt; idx++) {
                if (l != 0) {
                  factor[idx] = std::pow(rmesh[idx], l + 2) *
                                solution.zlr(idx, l_p, ispin) *
                                solution.jlr(idx, l_pp, ispin);
                } else {
                  factor[idx] = rmesh[idx] * rmesh[idx] *
                                solution.zlr(idx, l_p, ispin) *
                                solution.jlr(idx, l_pp,
                                             ispin); // Almost certain that this is already the complex conjugate
                }
              }

              ss_terms = lsms::simpson_nonuniform(rmesh, factor, jmt);

              std::printf("%2d %2d: %18.10e %18.10e\n", k_p, k_pp,
                          std::real(ss_terms), std::imag(ss_terms));

            }

            for (auto idx = 0; idx < jmt; idx++) {
              factor[idx] = std::pow(rmesh[idx], l + 2) *
                            solution.zlr(idx, l_p, ispin) *
                            solution.zlr(idx, l_pp, ispin);
            }

            std::complex<double> ms_terms = lsms::simpson_nonuniform(rmesh, factor, jmt);

            std::printf("%2d %2d: %18.10e %18.10e\n", k_p, k_pp,
                        std::real(ms_terms), std::imag(ms_terms));

            // Z * tau * Z^* - Z * J^*
            value += gaunt_coeff *
                     (-ss_terms + ms_terms * tau00_l(k_p, 0) // TODO: fix the array here
                     ) * std::sqrt(4.0 * M_PI) / (2.0 * l + 1.0);


            //std::printf("%d %d, %d %d, %d %d: %lf %lf\n", l, m, l_p, m_p, l_pp, m_pp, value, gaunt_coeff);

          }
        }

      }

      multi_mom_e[k] = value * 2.0 / M_PI;

      k++;

    }

    std::printf("DOSCK2 %lf %lf\n", std::real(multi_mom_e[0]), std::imag(multi_mom_e[0]));

  }


  return multi_mom_e;
}

void lsms::multi_moms::integrateMultipoleMoments(
    LSMSSystemParameters &lsms, int is, int ie, int nume, int lmax_mom,
    std::complex<double> energy, std::complex<double> dele1,
    std::vector<std::complex<double>> &mm_e, AtomData &atom) {
  int kmax_mom = lsms::get_kmax(lmax_mom);


  for (int isp = 0; isp < lsms.n_spin_cant * lsms.n_spin_cant; isp++) {
    // Only non-spin polarized and non-relativistic case
    for (int k = 0; k < kmax_mom; k++) {
      atom.multipole_moms[k] += std::imag(mm_e[k] * dele1);
    }
  }


  // if (ie == nume - 1) {
  //  for (int k = 0; k < kmax_mom; k++) {
  //    // Only non-spin polarized and non-relativistic case
  //   atom.multipole_moms[k] += std::imag(mm_e[k] * dele1);
  //}
  //\}
}
