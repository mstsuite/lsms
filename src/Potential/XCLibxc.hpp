//
// Created by F.Moitzi on 24.04.2022.
//

#ifndef LSMS_XCLIBXC_HPP
#define LSMS_XCLIBXC_HPP

#include "XCBase.hpp"

#ifdef USE_LIBXC

#include <xc.h>

namespace lsms {

class XCFuncType {
 private:
  bool isInitialised{false};
  xc_func_type _func_type;

 public:
  XCFuncType() = default;

  XCFuncType(int nspin, int type) : _func_type{} {
    isInitialised = true;

    if (xc_func_init(&_func_type, type, nspin) != 0) {
      isInitialised = false;
      throw std::runtime_error(
          "Specified `libxc` functional type doesn't exists");
    }
  }

  ~XCFuncType() {
    if (isInitialised) {
      //    xc_func_end(&_func_type);
    }
  }

  const xc_func_type &get_functional() const;
};

class XCLibxc : public XCBase {
 private:
  std::vector<XCFuncType> functionals;  // functionals
  std::size_t numFunctionals;
  bool needGradients;  // the functional needs gradients of the density (for
  // GGAs)
  bool needLaplacian;  // need laplacians of the density (for MetaGGAs)
  bool needKineticEnergyDensity;  // for MetaGGAs
  bool needExactExchange;         // for Hybrid Functionals

 public:
  XCLibxc(int nSpin, std::vector<int> xcFunctional);

  XCLibxc(int nSpin, int xcFunctional[3]);

  void evaluate(std::vector<Real> &rMesh,
                std::vector<Real> &drMesh,
                const Matrix<Real> &rhoIn,
                int jmt, Matrix<Real> &xcEnergyOut,
                Matrix<Real> &xcPotOut) override;

  void evaluate(const Real rhoIn[2], Real &xcEnergyOut,
                Real xcPotOut[2]) override;

  const std::vector<XCFuncType> &get_functionals() const;


  template <typename Cnt1, typename Cnt2>
  void addData(const Cnt1& yData, Cnt2 xData) // is pass-by-value intended?
  {
    using std::begin;
    using std::end;
    typedef decltype(*begin(yData)) T;
  }

};

}  // namespace lsms

#endif

#endif  // LSMS_XCLIBXC_HPP
