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

      /*
       *
       *                   L       4pi  ->   ->   * ->     ->
       *  cgnt(j,L',L") = C      = int do*Y (o )*Y (o )*Y (o )
       *                   L',L"    0      L      L'     L"
       *
       *  L', L" <= lmax
       *
       *  L = kj3(j,L',L"), j=1,2,...,nj3(L',L") <= lmax+1
       *
       *  Important is that the second index is the complex one
       *
       */
      double table(int k1, int k2, int k3) const;

    };

  }

}

#endif //LSMS_GAUNTFACTOR_HPP
