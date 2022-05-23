//
// Created by F.Moitzi on 24.04.2022.
//

#include "XCBase.hpp"

lsms::XCBase::XCBase(int nSpin, std::vector<int> xcFunctional)
    : _nSpin{nSpin}, _xcFunctional{xcFunctional} {}

lsms::XCBase::XCBase(int nSpin, int xcFunctional[3])
    : _nSpin{nSpin}, _xcFunctional(xcFunctional, xcFunctional + 3) {}
