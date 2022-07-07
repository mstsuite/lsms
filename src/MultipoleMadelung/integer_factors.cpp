//
// Created by F.Moitzi on 24.12.2021.
//

#include "integer_factors.hpp"

std::vector<int> lsms::get_lofk(int l0) {
  std::vector<int> lofk;

  for (auto l = 0; l <= l0; l++) {
    for (auto m = -l; m <= l; m++) {
      lofk.emplace_back(l);
    }
  }

  return lofk;
}

std::vector<int> lsms::get_mofk(int l0) {
  std::vector<int> mofk;

  for (auto l = 0; l <= l0; l++) {
    for (auto m = -l; m <= l; m++) {
      mofk.emplace_back(m);
    }
  }

  return mofk;
}

std::vector<int> lsms::get_lofj(int l0) {
  std::vector<int> lofj;

  for (auto l = 0; l <= l0; l++) {
    for (auto m = 0; m <= l; m++) {
      lofj.emplace_back(l);
    }
  }

  return lofj;
}

std::vector<int> lsms::get_mofj(int l0) {
  std::vector<int> mofj;

  for (auto l = 0; l <= l0; l++) {
    for (auto m = 0; m <= l; m++) {
      mofj.emplace_back(m);
    }
  }

  return mofj;
}

std::vector<int> lsms::get_kofj(int l0) {
  std::vector<int> kofj;

  for (auto l = 0; l <= l0; l++) {
    for (auto m = 0; m <= l; m++) {
      kofj.emplace_back((l + 1) * (l + 1) - l + m - 1);
    }
  }

  return kofj;
}
