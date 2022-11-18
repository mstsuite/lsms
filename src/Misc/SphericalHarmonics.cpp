//
// Created by F.Moitzi on 20.09.2022.
//

#include "SphericalHarmonics.hpp"

int lsms::pvt(int l, int m) { return ((m) + ((l) * ((l) + 1)) / 2); }

int lsms::yvr(int l, int m) { return ((m) + (l) + ((l) * (l))); }

lsms::SphericalHarmonics::SphericalHarmonics(int lmax) : _lmax{lmax} {
  auto sizeP = (lmax + 1) * (lmax + 2) / 2;

  P.resize(sizeP);
  A.resize(sizeP);
  B.resize(sizeP);

  // Precompute coefficients a_l^m and b_l^m for all l <= L, m <= l

  for (int l = 2; l <= lmax; l++) {
    double ls = l * l;

    double lm1s = (l - 1) * (l - 1);

    for (int m = 0; m < l - 1; m++) {
      double ms = m * m;

      A[pvt(l, m)] = std::sqrt((4.0 * ls - 1.0) / (ls - ms));
      B[pvt(l, m)] = -std::sqrt((lm1s - ms) / (4.0 * lm1s - 1.0));
    }
  }
}

void lsms::SphericalHarmonics::computeP(int lmax, double cos_theta) {
  if (_lmax < lmax) {
    throw std::runtime_error("The arrays are too small!!!");
  }

  const double sintheta = sqrt(1.0 - cos_theta * cos_theta);
  double temp = 0.39894228040143267794;  // = sqrt(0.5/M_PI)
  P[pvt(0, 0)] = temp;

  if (lmax > 0) {
    constexpr double SQRT3 = 1.7320508075688772935;
    P[pvt(1, 0)] = cos_theta * SQRT3 * temp;
    constexpr double SQRT3DIV2 = -1.2247448713915890491;
    temp = SQRT3DIV2 * sintheta * temp;
    P[pvt(1, 1)] = temp;

    for (int l = 2; l <= lmax; l++) {
      for (int m = 0; m < l - 1; m++) {
        P[pvt(l, m)] = A[pvt(l, m)] * (cos_theta * P[pvt(l - 1, m)] +
                                       B[pvt(l, m)] * P[pvt(l - 2, m)]);
      }
      P[pvt(l, l - 1)] = cos_theta * sqrt(2 * (l - 1) + 3) * temp;
      temp = -sqrt(1.0 + 0.5 / l) * sintheta * temp;
      P[pvt(l, l)] = temp;
    }
  }
}

void lsms::SphericalHarmonics::computeYlm(
    int lmax, std::vector<double> vec, std::vector<std::complex<double>> &Ylm) {
  int kmax = (lmax + 1) * (lmax + 1);

  double q2 = vec[0] * vec[0] + vec[1] * vec[1];
  double r = std::sqrt(q2 + vec[2] * vec[2]);
  double q = std::sqrt(q2);

  if (r < TOL) {
    Ylm[0] = 0.28209479177387814347403972578038629292;
    std::fill(Ylm.begin() + 1, Ylm.end(), std::complex<double>(0.0, 0.0));

  } else {
    using namespace std::complex_literals;

    double cos_theta = vec[2] / r;

    std::complex<double> iphi = 1i * atan2(vec[1], vec[0]);

    std::complex<double> mphi;

    computeP(lmax, cos_theta);

    for (int l = 0; l <= lmax; l++) {
      Ylm[yvr(l, 0)] = P[pvt(l, 0)] * M_SQRT1_2;

      for (int m = 1; m <= l; m++) {
        mphi = iphi * std::complex<double>(m);

        Ylm[yvr(l, m)] = P[pvt(l, m)] * std::exp(mphi) * M_SQRT1_2;
        Ylm[yvr(l, -m)] = std::pow(-1, m) * std::conj(Ylm[yvr(l, m)]);
      }
    }
  }

  //  for (size_t l = 0; l <= L; l++) {
  //    Y[YVR(l, 0)] = P[PVT(l, 0)] * 0.5 * M_SQRT2;
  //  }
}
