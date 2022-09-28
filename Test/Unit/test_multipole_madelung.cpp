//
// Created by F.Moitzi on 23.06.2022.
//

#include <gtest/gtest.h>

#include <iostream>
#include <vector>

#include "Madelung/Madelung.hpp"
#include "Main/SystemParameters.hpp"
#include "Misc/Coeficients.hpp"
#include "Misc/Indices.hpp"
#include "MultipoleMadelung/MultipoleMadelung.hpp"
#include "MultipoleMadelung/calculateMultipoleMadelung.hpp"
#include "accel_common.hpp"

template <typename T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &v) {
  os << "[";
  for (int i = 0; i < v.size(); ++i) {
    os << v[i];
    if (i != v.size() - 1) os << ", ";
  }
  os << "]";
  return os;
}

namespace multipole_tests {

SphericalHarmonicsCoeficients sphericalHarmonicsCoeficients;
GauntCoeficients gauntCoeficients;
IFactors iFactors;

TEST(MultiMadelungsTestSuite, Structure1) {
  int num_atoms = 2;
  int lmax = 3;

  LSMSSystemParameters lsms;
  lsms.global.iprint = -1;

  AngularMomentumIndices::init(2 * lmax);
  SphericalHarmonicsCoeficients::init(2 * lmax);
  GauntCoeficients::init(lsms);
  IFactors::init(lsms, lmax);

  LocalTypeInfo local;
  local.setNumLocal(2);
  local.global_id = {0, 1};
  //
  LocalTypeInfo local_modern = local;
  //
  CrystalParameters crystal;
  crystal.resize(num_atoms);
  crystal.num_atoms = num_atoms;
  crystal.num_types = num_atoms;

  double scale = 1;
  double a = 1.0;

  crystal.bravais(0, 0) = 1.0 * a * scale;
  crystal.bravais(0, 1) = 0.0;
  crystal.bravais(0, 2) = 0.0;

  crystal.bravais(1, 0) = 0.0;
  crystal.bravais(1, 1) = 1.0 * a * scale;
  crystal.bravais(1, 2) = 0.0;

  crystal.bravais(2, 0) = 0.0;
  crystal.bravais(2, 1) = 0.0;
  crystal.bravais(2, 2) = 1.0 * a * scale;
  //
  crystal.position(0, 0) = 0.0;
  crystal.position(1, 0) = 0.0;
  crystal.position(2, 0) = 0.0;

  crystal.position(0, 1) = 0.5 * a;
  crystal.position(1, 1) = 0.5 * a;
  crystal.position(2, 1) = 0.5 * a;

  calculateMadelungMatrices(lsms, crystal, local);

  std::cout << std::endl;
  for (int i = 0; i < num_atoms; i++) {
    std::cout << local.atom[i].madelungMatrix;
    std::cout << std::endl;
  }

  calculateMultiMadelungMatrices(lsms, crystal, local_modern, lmax);

  std::cout << std::endl;
  for (int i = 0; i < num_atoms; i++) {
    std::cout << local_modern.atom[i].madelungMatrix;
    std::cout << std::endl;
  }

  for (int i = 0; i < num_atoms; i++) {
    for (int j = 0; j < num_atoms; j++) {
      EXPECT_NEAR(local.atom[i].madelungMatrix[j],
                  local_modern.atom[i].madelungMatrix[j], 1.0e-12);
    }
  }

  EXPECT_NEAR(0.28209479177387842, lsms.dl_factor(0, 0), 1e-12);
  EXPECT_NEAR(9.4031597257959593E-002, lsms.dl_factor(0, 1), 1e-12);
  EXPECT_NEAR(-9.4031597257959454E-002, lsms.dl_factor(0, 2), 1e-12);
  EXPECT_NEAR(1.8806319451591908E-002, lsms.dl_factor(0, 3), 1e-12);
  EXPECT_NEAR(-1.8806319451591908E-002, lsms.dl_factor(0, 4), 1e-12);

  EXPECT_NEAR(-10.057957687339862,
              std::real(local_modern.atom[0].multipoleMadelung(0, 0)), 1e-12);
  EXPECT_NEAR(0.0, std::imag(local_modern.atom[0].multipoleMadelung(0, 0)),
              1e-12);
}

TEST(MultiMadelungsTestSuite, Structure2) {
  int num_atoms = 3;
  int lmax = 3;

  LSMSSystemParameters lsms;
  lsms.global.iprint = -1;

  AngularMomentumIndices::init(2 * lmax);
  SphericalHarmonicsCoeficients::init(2 * lmax);
  GauntCoeficients::init(lsms);
  IFactors::init(lsms, lmax);

  LocalTypeInfo local;
  local.setNumLocal(2);
  local.global_id = {0, 2};
  //
  LocalTypeInfo local_modern = local;
  //
  CrystalParameters crystal;
  crystal.resize(num_atoms);
  crystal.num_atoms = num_atoms;
  crystal.num_types = num_atoms;

  double scale = 1;
  double a = 1.0;

  crystal.bravais(0, 0) = 1.1;
  crystal.bravais(1, 0) = 0.2;
  crystal.bravais(2, 0) = 0.5;

  crystal.bravais(0, 1) = -0.1;
  crystal.bravais(1, 1) = 1.1;
  crystal.bravais(2, 1) = 0.1;

  crystal.bravais(0, 2) = 0.0;
  crystal.bravais(1, 2) = 0.0;
  crystal.bravais(2, 2) = 1.6;

  crystal.position(0, 0) = 0.0;
  crystal.position(1, 0) = 0.0;
  crystal.position(2, 0) = 0.1;

  crystal.position(0, 1) = 0.5;
  crystal.position(1, 1) = 0.5;
  crystal.position(2, 1) = 0.4;

  crystal.position(0, 2) = 0.75;
  crystal.position(1, 2) = 0.65;
  crystal.position(2, 2) = 0.4;

  lsms::MultipoleMadelung madelung(lsms, crystal, local_modern, lmax);

  // Test scaling factor
  EXPECT_NEAR(0.97651137065361771, madelung.getScalingFactor(), 1e-12);

  // Test of all madelung values
  EXPECT_NEAR(-0.2223187728889e+1, local_modern.atom[0].madelungMatrix[0],
              1e-12);
  EXPECT_NEAR(-0.2103456956161, local_modern.atom[0].madelungMatrix[1], 1e-12);
  EXPECT_NEAR(-0.1084495831085, local_modern.atom[0].madelungMatrix[2], 1e-12);

  EXPECT_NEAR(-0.1084495831085, local_modern.atom[1].madelungMatrix[0], 1e-12);
  EXPECT_NEAR(0.1305464656817e1, local_modern.atom[1].madelungMatrix[1], 1e-12);
  EXPECT_NEAR(-0.2223187728889e1, local_modern.atom[1].madelungMatrix[2],
              1e-12);

  // Test dl factors for lmax = 3 (kmax = 16, jmax = 10)
  {
    const std::vector<double> dl_factor{
        0.28209479177387842,     9.4031597257959454E-002,
        9.4031597257959579E-002, 9.4031597257959454E-002,
        1.8806319451591877E-002, 1.8806319451591905E-002,
        1.8806319451591908E-002, 1.8806319451591905E-002,
        1.8806319451591877E-002, 2.6866170645131245E-003,
        2.6866170645131284E-003, 2.6866170645131276E-003,
        2.6866170645131302E-003, 2.6866170645131276E-003,
        2.6866170645131284E-003, 2.6866170645131245E-003};

    for (int i = 0; i < dl_factor.size(); i++) {
      EXPECT_NEAR(dl_factor[i], lsms.dl_factor(i, 0), 1e-12);
    }
  }

  {
    const std::vector<double> dl_factor{
        9.4031597257959593E-002, 2.4278854013157353E-002,
        2.8034805800224084E-002, 2.4278854013157353E-002,
        4.1038753538304908E-003, 5.1910373406135095E-003,
        5.5059265574105955E-003, 5.1910373406135095E-003,
        4.1038753538304908E-003, 5.1703969513536556E-004,
        6.7696386864420153E-004, 7.5686861429983654E-004,
        7.8169054358676125E-004, 7.5686861429983654E-004,
        6.7696386864420153E-004, 5.1703969513536556E-004};

    for (int i = 0; i < dl_factor.size(); i++) {
      EXPECT_NEAR(dl_factor[i], lsms.dl_factor(i, 1), 1e-12);
    }
  }

  {
    const std::vector<double> dl_factor{
        -9.4031597257959454E-002, -3.4335484624283527E-002,
        -2.4278854013157353E-002, -1.4017402900111983E-002,
        -7.1081206207641032E-003, -5.8037561836757640E-003,
        -4.4955702089649052E-003, -3.1788481800593097E-003,
        -1.8353088524701872E-003, -1.0340793902707298E-003,
        -8.9553902150437618E-004, -7.5686861429983567E-004,
        -6.1798063578732195E-004, -4.7868574213659313E-004,
        -3.3848193432210006E-004, -1.9542263589668885E-004};

    for (int i = 0; i < dl_factor.size(); i++) {
      EXPECT_NEAR(dl_factor[i], lsms.dl_factor(i, 2), 1e-12);
    }
  }

  {
    const std::vector<double> dl_factor{
        -2.6866170645131245E-003, -1.0340793902707294E-003,
        -5.1703969513536556E-004, -1.9542263589668882E-004,
        -2.2046646677429237E-004, -1.3943523653932001E-004,
        -8.0502971350495422E-005, -4.0251485675247677E-005,
        -1.5213631571083300E-005, -3.2684123568626300E-005,
        -2.3111165412514750E-005, -1.5581544542482675E-005,
        -9.8546340435222894E-006, -5.6895756177928420E-006,
        -2.8447878088964045E-006, -1.0752287250125931E-006};

    for (int i = 0; i < dl_factor.size(); i++) {
      EXPECT_NEAR(dl_factor[i], lsms.dl_factor(i, 9), 1e-12);
    }
  }

  {
    using namespace std::complex_literals;

    const std::vector<std::complex<double>> dl_matrix{
        -7.8809953027095991,
        1.88966212902900115E-017,
        0.0000000000000000,
        -1.88966212902900115E-017,
        -1.8097382618160811 + 0.27561894928136865i,
        1.8021011211224087 - 0.93990993338631279i,
        -4.2970575853540689,
        -1.8021011211224087 - 0.93990993338631279i,
        -1.8097382618160811 - 0.27561894928136865i,
        -1.92304485683264783E-015 - 6.41014952277549244E-016i,
        -2.20348889845407568E-016i,
        -3.20507476138774622E-016,
        0.0000000000000000,
        3.20507476138774622E-016,
        2.20348889845407568E-016i,
        1.92304485683264783E-015 - 6.41014952277549244E-016i};

    for (int j = 0; j < dl_matrix.size(); j++) {
      auto ref_res = dl_matrix[j];
      auto res = local_modern.atom[0].multipoleMadelung(j, 0);

      EXPECT_NEAR(std::real(ref_res), std::real(res), 1e-12);
      EXPECT_NEAR(std::imag(ref_res), std::imag(res), 1e-12);
    }
  }
}

}  // namespace multipole_tests