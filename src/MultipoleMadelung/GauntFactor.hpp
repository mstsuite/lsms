//
// Created by F.Moitzi on 25.01.2022.
//

#ifndef LSMS_GAUNTFACTOR_HPP
#define LSMS_GAUNTFACTOR_HPP

#include "NDArray.hpp"

namespace lsms {

  namespace math {

    class GauntFactor {

    public:

      lsms::NDArray<double, 3> cgnt;
      lsms::NDArray<int, 2> nj3;
      lsms::NDArray<int, 3> kj3;


    public:

      explicit GauntFactor(int lmax);

      double table(int k1, int k2, int k3) const;

    };

  }

}

#endif //LSMS_GAUNTFACTOR_HPP
