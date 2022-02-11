//
// Created by F.Moitzi on 25.01.2022.
//

#ifndef LSMS_GAUNTFACTOR_HPP
#define LSMS_GAUNTFACTOR_HPP

#include "Array3d.hpp"
#include "NDArray.hpp"

namespace lsms {

  namespace math {

    class GauntFactorBase {
    public:
      virtual double table(int l1, int m1, int l2, int m2, int l3,
                           int m3) const = 0;
    };

    class GauntFactor : public GauntFactorBase {
    public:
      lsms::NDArray<double, 3> cgnt;
      lsms::NDArray<int, 2> nj3;
      lsms::NDArray<int, 3> kj3;

    public:
      explicit GauntFactor(int lmax);

      double table(int l1, int m1, int l2, int m2, int l3, int m3) const;
    };

    class GauntFactorWrapper : public GauntFactorBase {
    private:
      const Array3d<double> &cgnt;

    public:
      explicit GauntFactorWrapper(const Array3d<double> & cgnt) : cgnt( cgnt ) {};

      double table(int l1, int m1, int l2, int m2, int l3, int m3) const;
    };

  }  // namespace math

}  // namespace lsms

#endif  // LSMS_GAUNTFACTOR_HPP
