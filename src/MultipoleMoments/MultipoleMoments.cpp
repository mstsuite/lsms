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
    const lsms::math::GauntFactor &gaunt_factor,
    unsigned int lmax_mm,
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

        for (auto k_pp{0}; k_pp < kmax; k_pp++) {

          auto l_pp = lofk[k_pp];

          std::vector<std::complex<double>> factor(jmt);

          for (auto idx = 0; idx < jmt; idx++) {

            factor[idx] = std::pow(rmesh[idx], l)
                          * solution.zlr(idx, l_pp, ispin)
                          * std::conj(solution.zlr(idx, l_p, ispin));

            std::complex<double> ss_terms = 4 * M_PI * lsms::simpson_nonuniform(rmesh, factor, jmt);

            factor[idx] = std::pow(rmesh[idx], l)
                          * solution.zlr(idx, l_pp, ispin)
                          * std::conj(solution.zlr(idx, l_p, ispin));

            std::complex<double> ms_terms = 4 * M_PI * lsms::simpson_nonuniform(rmesh, factor, jmt);

            /*
             *          4pi  _  m1 _     m2* _     m3 _
             *   clll = int do Y  (o) * Y   (o) * Y  (o)
             *            0      l1       l2        l3
             *
             *       L3
             *    = C
             *       L1,L2
             *
             *   l1 = 0, 1, ..., lmax_1
             *   l2 = 0, 1, ..., lmax_2
             *   l3 = 0, 1, ..., lmax_1+lmax_2
             *
             *   l2 seems to be the complex number for the conjugates
             *
             */


            // Z * tau * Z^* - Z * J^*
            value += gaunt_factor.table(k, k_p, k_pp) *
                     (-ss_terms + ms_terms * tau00_l(k_p, k_pp));


          }
        }
      }

      multi_mom_e[k] = value;

      k++;
    }
  }

  return multi_mom_e;

}

void lsms::multi_moms::integrateMultipoleMoments(LSMSSystemParameters &lsms,
                                                 int is,
                                                 int ie,
                                                 int nume,
                                                 int lmax_mom,
                                                 std::complex<double> energy,
                                                 std::complex<double> dele1,
                                                 std::vector<std::complex<double>> &mm_e,
                                                 AtomData &atom) {

  int kmax_mom = lsms::get_kmax(lmax_mom);

  std::complex<double> ede = energy * dele1;

  for (int isp = 0; isp < lsms.n_spin_cant * lsms.n_spin_cant; isp++) {

    // Only non-spin polarized and non-relativistic case
    for (int k = 0; k < kmax_mom; k++) {
      atom.multipole_moms[k] += std::imag(mm_e[k]);
    }

  }

  if (ie == nume - 1) {
    for (int k = 0; k < kmax_mom; k++) {

      // Only non-spin polarized and non-relativistic case
      atom.multipole_moms[k] += std::imag(mm_e[k]);

    }
  }


}

