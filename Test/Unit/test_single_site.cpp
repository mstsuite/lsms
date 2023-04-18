//
// Created by F.Moitzi on 15.05.2022.
//

#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/printf.h>
#include <gtest/gtest.h>
//#include <matplot/matplot.h>

#include <cmath>
#include <iostream>

#include "Array3d.hpp"
#include "Matrix.hpp"
#include "PhysicalConstants.hpp"
#include "RadialSrelSolver.hpp"
#include "accel_common.hpp"

namespace single_site_tests {

extern "C" {

void single_scatterer_nonrel_(int *nrelv, double *clight, int *lmax, int *kkrsz,
                              Complex *energy, Complex *prel, Complex *pnrel,
                              double *vr, double *r_mesh, double *h, int *jmt,
                              int *jws, Complex *tmat_l, Complex *matom,
                              Complex *zlr, Complex *jlr, double *r_sph,
                              int *iprpts, int *iprint, char *istop,
                              int istop_len);
}

TEST(SingleSiteTests, NewRadialSolver) {
  int n = 1000;
  double clight = 274.071998412000027656176826;
  int lmax = 3;
  std::complex<double> energy{(-0.299986, 0.00359092)};

  Matrix<Real> vr(n, 1);
  std::vector<Real> r_mesh(n);
  std::vector<Real> dr_mesh(n);

  Real h = 0.013;
  Real rmin = 0.0001;
  Real z = 15;

  int jmt = 780;
  int jws = 780;

  for (auto ir = 0; ir < n; ir++) {
    r_mesh[ir] = rmin * std::exp(h * ir);
    dr_mesh[ir] = rmin * h * std::exp(h * ir);
  }

  Real r_sph = r_mesh[jmt - 1];

  for (auto ir = 0; ir < n; ir++) {
    vr(ir, 0) = -z / std::pow(r_mesh[ir], 0.8);
  }

  for (auto ir = 0; ir < n; ir++) {
    vr(ir, 0) -= vr(jmt - 1, 0);
  }

  for (auto ir = jmt; ir < n; ir++) {
    vr(ir, 0) = 0.0;
  }

  double beta;

  int kappa = 1;

  std::vector<Complex> P(n);
  std::vector<Complex> Q(n);

  lsms::RadialSrelOut(energy, kappa, P.data(), Q.data(), r_mesh.data(),
                      dr_mesh.data(), &vr(0, 0), n, 0, &beta);

  std::vector<Real> y(n);
  Real result;

  // auto ax1 = matplot::subplot(2, 1, 0);

  for (auto ir = 0; ir < n; ir++) {
    y[ir] = std::real(P[ir]) * r_mesh[ir];
  }

  result = *std::max_element(y.begin(), y.end(), [](int a, int b) {
    return std::abs(a) < std::abs(b);
  });

  for (auto ir = 0; ir < n; ir++) {
    y[ir] /= result;
  }

  // matplot::plot(r_mesh, y, "--xr");
  // matplot::hold(matplot::on);

  // matplot::xlim({0, r_sph + 0.5});

  auto name = std::string(
      ::testing::UnitTest::GetInstance()->current_test_info()->name());

  // matplot::save(name + ".png");
}

TEST(SingleSiteTests, BasicSolver) {
  int n = 1000;

  // Set parameters
  int nrelv = 0;
  double clight = 274.071998412000027656176826;
  int lmax = 3;
  int kkrsz = (lmax + 1) * (lmax + 1);

  std::complex<double> energy{(-0.299986, 0.00359092)};

  Complex prel = std::sqrt(energy * (1.0 + energy * c2inv));
  Complex pnrel = std::sqrt(energy);

  Matrix<Real> vr(n, 1);
  std::vector<Real> r_mesh(n);

  Real h = 0.013;
  Real rmin = 0.0001;
  Real z = 15;

  int jmt = 780;
  int jws = 780;

  for (auto ir = 0; ir < n; ir++) {
    r_mesh[ir] = rmin * std::exp(h * ir);
  }

  Real r_sph = r_mesh[jmt - 1];

  for (auto ir = 0; ir < n; ir++) {
    vr(ir, 0) = -z / std::pow(r_mesh[ir], 0.8);
  }

  for (auto ir = 0; ir < n; ir++) {
    vr(ir, 0) -= vr(jmt - 1, 0);
  }

  for (auto ir = jmt; ir < n; ir++) {
    vr(ir, 0.0) = 0.0;
  }

  for (auto ir = 0; ir < n; ir++) {
    vr(ir, 0) = vr(ir, 0) * r_mesh[ir];
  }

  // Init arrays and matrix
  Array3d<Complex> tmat_l(kkrsz, kkrsz, 1);
  Matrix<Complex> matom(lmax + 1, 1);
  Array3d<Complex> zlr(n, lmax + 1, 1);
  Array3d<Complex> jlr(n, lmax + 1, 1);

  int iprpts = n;
  int iprint = 10;
  char istop[32];

  // Solver single site equation
  single_scatterer_nonrel_(&nrelv, &clight, &lmax, &kkrsz, &energy, &prel,
                           &pnrel, &vr(0, 0), &r_mesh[0], &h, &jmt, &jws,
                           &tmat_l(0, 0, 0), &matom(0, 0), &zlr(0, 0, 0),
                           &jlr(0, 0, 0), &r_sph, &iprpts, &iprint, istop, 32);

  for (auto ir = 0; ir < n; ir++) {
    fmt::printf("%20.16f %40.24f %40.24f\n", r_mesh[ir],
                std::real(zlr(ir, 0, 0)), std::imag(zlr(ir, 0, 0)));
  }

  std::vector<Real> y(n);
  Real result;

  // auto ax1 = matplot::subplot(2, 1, 0);

  for (auto ir = 0; ir < n; ir++) {
    y[ir] = std::real(zlr(ir, 0, 0));
  }

  result = *std::max_element(y.begin(), y.end(), [](int a, int b) {
    return std::abs(a) < std::abs(b);
  });

  for (auto ir = 0; ir < n; ir++) {
    y[ir] /= result;
  }

  // matplot::plot(r_mesh, y, "--xr");
  // matplot::hold(matplot::on);

  for (auto ir = 0; ir < n; ir++) {
    y[ir] = std::real(zlr(ir, 1, 0));
  }

  result = *std::max_element(y.begin(), y.end(), [](int a, int b) {
    return std::abs(a) < std::abs(b);
  });

  for (auto ir = 0; ir < n; ir++) {
    y[ir] /= result;
  }

  // matplot::plot(r_mesh, y, "--xg");
  // matplot::hold(matplot::on);

  for (auto ir = 0; ir < n; ir++) {
    y[ir] = std::real(zlr(ir, 2, 0));
  }

  result = *std::max_element(y.begin(), y.end(), [](int a, int b) {
    return std::abs(a) < std::abs(b);
  });

  for (auto ir = 0; ir < n; ir++) {
    y[ir] /= result;
  }

  // matplot::plot(r_mesh, y, "--xb");
  // matplot::hold(matplot::on);

  // matplot::xlim({0, r_sph + 0.5});

  // auto ax2 = matplot::subplot(2, 1, 1);

  for (auto ir = 0; ir < n; ir++) {
    y[ir] = std::real(jlr(ir, 0, 0)) * r_mesh[ir];
  }

  // matplot::plot(r_mesh, y, "--xr");
  // matplot::hold(matplot::on);

  for (auto ir = 0; ir < n; ir++) {
    y[ir] = std::real(jlr(ir, 1, 0)) * r_mesh[ir];
  }

  // matplot::plot(r_mesh, y, "--xg");

  // matplot::xlim({0, r_sph + 0.5});

  auto name = std::string(
      ::testing::UnitTest::GetInstance()->current_test_info()->name());

  // matplot::save(name + ".png");
}
}  // namespace single_site_tests