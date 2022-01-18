//
// Created by F.Moitzi on 24.01.2022.
//

#include "MultipoleMoments.hpp"


std::vector<std::complex<double>> lsms::multi_moms::multipole_mom_e(
    const NonRelativisticSingleScattererSolution &solution,
    const Matrix<std::complex<double>> &tau00_l,
    const lsms::math::GauntFactor &gaunt_factor,
    unsigned int lmax) {

  int ispin = 0;


  auto kmax = lsms::get_kmax(lmax);
  std::vector<std::complex<double>> multi_mom_e(kmax, 0.0);

  auto k = 0;

  for (auto l{0}; l <= lmax; l++) {
    for (auto m = -l; m <= l; m++) {
      multi_mom_e[k] = term(ispin, k, solution, tau00_l, kmax, gaunt_factor);
      k++;
    }
  }

  return multi_mom_e;

}


std::complex<double> lsms::multi_moms::term(int ispin,
                                            int k,
                                            const NonRelativisticSingleScattererSolution &solution,
                                            const Matrix<std::complex<double>> &tau00_l,
                                            int kmax,
                                            const lsms::math::GauntFactor &gaunt_factor) {
  auto lofk = lsms::get_lofk(kmax);
  auto mofk = lsms::get_mofk(kmax);

  int lmax = lofk[kmax];
  int m = lofk[k];
  int l = mofk[k];

  //
  std::complex<double> value{0.0};

  // Local parameters
  auto &rmesh = solution.atom->r_mesh;
  auto &jmt = solution.atom->jmt;


  lsms::NDArray<std::complex<double>, 1> ss_terms(lmax); // ss_terms (l_p)
  lsms::NDArray<std::complex<double>, 2> ms_terms(lmax, lmax); // ms_terms (l_p, l_pp)

  // Single site terms
  for (auto l_p{0}; l_p < lmax; l_p++) {

    std::vector<std::complex<double>> factor(jmt);

    for (auto idx = 0; idx < jmt; idx++) {

      factor[idx] = std::pow(rmesh[idx], l)
                    * solution.zlr(idx, l_p, ispin)
                    * std::conj(solution.jlr(idx, l_p, ispin));

      ss_terms(l_p) = lsms::simpson_nonuniform(rmesh, factor, jmt);

    }
  }

  // Multiple scattering
  for (auto l_p{0}; l_p < lmax; l_p++) {
    for (auto l_pp{0}; l_pp < lmax; l_pp++) {

      std::vector<std::complex<double>> factor(jmt);

      for (auto idx = 0; idx < jmt; idx++) {

        factor[idx] = std::pow(rmesh[idx], l)
                      * solution.zlr(idx, l_p, ispin)
                      * std::conj(solution.zlr(idx, l_pp, ispin));

        ms_terms(l_p, l_pp) = lsms::simpson_nonuniform(rmesh, factor, jmt);

      }
    }
  }


  for (auto k_p{0}; k_p < kmax; k_p++) {

    auto l_p = lofk[k_p];
    auto m_p = mofk[k_p];


    for (auto k_pp{0}; k_pp < kmax; k_pp++) {

      auto l_pp = lofk[k_pp];
      auto m_pp = mofk[k_pp];


      // Single-site term
      if (k_pp == k_p) {
        value += gaunt_factor.table(k, k_pp, k_p) *
                 ss_terms(l_p);
      }

      // Multiple-scattering term
      value += gaunt_factor.table(k, k_pp, k_p) *
               ms_terms(l_p, l_pp) * tau00_l(k_p, k_pp);


    }
  }


  return 4 * M_PI * value;

}


