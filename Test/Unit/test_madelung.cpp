//
// Created by F.Moitzi on 23.06.2022.
//

#include <gtest/gtest.h>

#include <iostream>
#include <vector>

#include "Madelung/Madelung.hpp"
#include "Main/SystemParameters.hpp"
#include "MultipoleMadelung/calculateMultipoleMadelung.hpp"
#include "accel_common.hpp"

namespace madlung_tests {

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

TEST(MadelungsTestSuite, SupercellLarge) {
  int num_atoms = 2;

  LSMSSystemParameters lsms;
  lsms.global.iprint = -1;

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

  double scale = 100;
  double a = 5.3821038770852745;

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

  calculateMultiMadelungMatrices(lsms, crystal, local_modern);

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
}

TEST(MadelungsTestSuite, Supercell) {
  int num_atoms = 2;

  LSMSSystemParameters lsms;
  lsms.global.iprint = -1;

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

  double scale = 16;
  double a = 5.3821038770852745;

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

  calculateMultiMadelungMatrices(lsms, crystal, local_modern);

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
}

TEST(MadelungsTestSuite, BasicStructure) {
  int num_atoms = 4;

  LSMSSystemParameters lsms;
  lsms.global.iprint = -1;

  LocalTypeInfo local;
  local.setNumLocal(4);
  local.global_id = {0, 1, 2, 3};
  //
  LocalTypeInfo local_modern = local;
  //
  CrystalParameters crystal;
  crystal.resize(num_atoms);
  crystal.num_atoms = num_atoms;
  crystal.num_types = num_atoms;

  crystal.bravais(0, 0) = 2.0;
  crystal.bravais(0, 1) = 0.0;
  crystal.bravais(0, 2) = 0.0;

  crystal.bravais(1, 0) = 0.1;
  crystal.bravais(1, 1) = 2.0;
  crystal.bravais(1, 2) = 0.0;

  crystal.bravais(2, 0) = -0.1;
  crystal.bravais(2, 1) = 0.05;
  crystal.bravais(2, 2) = 4.0;
  //
  crystal.position(0, 0) = 0.0;
  crystal.position(1, 0) = 0.01;
  crystal.position(2, 0) = 0.0;

  crystal.position(0, 1) = 0.0;
  crystal.position(1, 1) = 0.0;
  crystal.position(2, 1) = 2.0;

  crystal.position(0, 2) = 0.0;
  crystal.position(1, 2) = 1.0;
  crystal.position(2, 2) = 0.15;

  crystal.position(0, 3) = 0.0;
  crystal.position(1, 3) = 1.0;
  crystal.position(2, 3) = 2.0;

  calculateMadelungMatrices(lsms, crystal, local);

  std::cout << std::endl;
  for (int i = 0; i < num_atoms; i++) {
    std::cout << local.atom[i].madelungMatrix;
    std::cout << std::endl;
  }

  calculateMultiMadelungMatrices(lsms, crystal, local_modern);

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
}

TEST(MadelungsTestSuite, DiffStructure) {
  int global_num_atoms = 4;
  int num_atoms = 2;

  LSMSSystemParameters lsms;
  lsms.global.iprint = -1;

  LocalTypeInfo local;
  local.setNumLocal(num_atoms);
  local.global_id = {0, 2};
  //
  LocalTypeInfo local_modern = local;
  //
  CrystalParameters crystal;
  crystal.resize(global_num_atoms);
  crystal.num_atoms = global_num_atoms;
  crystal.num_types = global_num_atoms;

  crystal.bravais(0, 0) = 2.0;
  crystal.bravais(0, 1) = 0.0;
  crystal.bravais(0, 2) = 0.0;

  crystal.bravais(1, 0) = 0.1;
  crystal.bravais(1, 1) = 2.0;
  crystal.bravais(1, 2) = 0.0;

  crystal.bravais(2, 0) = -0.1;
  crystal.bravais(2, 1) = 0.05;
  crystal.bravais(2, 2) = 4.0;
  //
  crystal.position(0, 0) = 0.0;
  crystal.position(1, 0) = 0.01;
  crystal.position(2, 0) = 0.0;

  crystal.position(0, 1) = 0.0;
  crystal.position(1, 1) = 0.0;
  crystal.position(2, 1) = 2.0;

  crystal.position(0, 2) = 0.0;
  crystal.position(1, 2) = 1.0;
  crystal.position(2, 2) = 0.15;

  crystal.position(0, 3) = 0.0;
  crystal.position(1, 3) = 1.0;
  crystal.position(2, 3) = 2.0;

  calculateMadelungMatrices(lsms, crystal, local);

  std::cout << std::endl;
  for (int i = 0; i < num_atoms; i++) {
    std::cout << local.atom[i].madelungMatrix;
    std::cout << std::endl;
  }

  calculateMultiMadelungMatrices(lsms, crystal, local_modern);

  std::cout << std::endl;
  for (int i = 0; i < num_atoms; i++) {
    std::cout << local_modern.atom[i].madelungMatrix;
    std::cout << std::endl;
  }

  for (int i = 0; i < num_atoms; i++) {
    for (int j = 0; j < global_num_atoms; j++) {
      EXPECT_NEAR(local.atom[i].madelungMatrix[j],
                  local_modern.atom[i].madelungMatrix[j], 1.0e-12);
    }
  }
}

TEST(MadelungsTestSuite, AlongatedStructure) {
  int global_num_atoms = 4;
  int num_atoms = 2;

  LSMSSystemParameters lsms;
  lsms.global.iprint = -1;

  LocalTypeInfo local;
  local.setNumLocal(num_atoms);
  local.global_id = {0, 2};
  //
  LocalTypeInfo local_modern = local;
  //
  CrystalParameters crystal;
  crystal.resize(global_num_atoms);
  crystal.num_atoms = global_num_atoms;
  crystal.num_types = global_num_atoms;

  crystal.bravais(0, 0) = 50.0;
  crystal.bravais(0, 1) = 0.0;
  crystal.bravais(0, 2) = 0.0;

  crystal.bravais(1, 0) = 0.0;
  crystal.bravais(1, 1) = 50.0;
  crystal.bravais(1, 2) = 0.0;

  crystal.bravais(2, 0) = 0.0;
  crystal.bravais(2, 1) = 0.0;
  crystal.bravais(2, 2) = 2.0;
  //
  crystal.position(0, 0) = 0.0;
  crystal.position(1, 0) = 0.01;
  crystal.position(2, 0) = 0.0;

  crystal.position(0, 1) = 0.0;
  crystal.position(1, 1) = 0.0;
  crystal.position(2, 1) = 2.0;

  crystal.position(0, 2) = 0.0;
  crystal.position(1, 2) = 1.0;
  crystal.position(2, 2) = 0.15;

  crystal.position(0, 3) = 0.0;
  crystal.position(1, 3) = 1.0;
  crystal.position(2, 3) = 2.0;

  calculateMultiMadelungMatrices(lsms, crystal, local_modern);

  std::cout << std::endl;
  for (int i = 0; i < num_atoms; i++) {
    std::cout << local_modern.atom[i].madelungMatrix;
    std::cout << std::endl;
  }
}

}  // namespace madlung_tests