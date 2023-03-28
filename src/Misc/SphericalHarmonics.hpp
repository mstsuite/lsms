//
// Created by F.Moitzi on 20.09.2022.
//

#ifndef LSMS_SPHERICALHARMONICS_HPP
#define LSMS_SPHERICALHARMONICS_HPP

#include <cmath>
#include <complex>
#include <vector>

namespace lsms {

constexpr static double TOL = 0.5 * 1.0e-12;

inline int pvt(int l, int m);

inline int yvr(int l, int m);

class SphericalHarmonics {
 public:
  explicit SphericalHarmonics(int lmax);

  /* Compute an entire set of Y_{l,m}(\theta,\phi) and store in array Y */
  void computeYlm(int lmax, std::vector<double> vec,
                  std::vector<std::complex<double>> &Ylm);

 private:
  int _lmax;

  std::vector<double> A;
  std::vector<double> B;
  std::vector<double> P;

  /**
   * Compute an entire set of P_l^m(x) and store in the array P
   *
   * @param lmax
   * @param X
   */
  void computeP(int lmax, double X);
};

}  // namespace lsms

#endif  // LSMS_SPHERICALHARMONICS_HPP
