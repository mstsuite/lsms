//
// Created by F.Moitzi on 24.04.2022.
//

#include "XCBase.hpp"

lsms::XCBase::XCBase(int n_spin, std::vector<int> xc_functional)
    : _nSpin{n_spin}, _xcFunctional{xc_functional} {}

lsms::XCBase::XCBase(int n_spin, int xc_functional[3])
    : _nSpin{n_spin}, _xcFunctional(xc_functional, xc_functional + 3) {}
