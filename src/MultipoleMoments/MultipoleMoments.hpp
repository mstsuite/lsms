//
// Created by F.Moitzi on 24.01.2022.
//

#ifndef LSMS_MULTIPOLEMOMENTS_HPP
#define LSMS_MULTIPOLEMOMENTS_HPP

#include <complex>

#include "SingleSite/SingleSiteScattering.hpp"

#include "utils.hpp"
#include "integrator.hpp"
#include "NDArray.hpp"
#include "integer_factors.hpp"
#include "common.hpp"
#include "GauntFactor.hpp"

namespace lsms {

  namespace multi_moms {

    using size_t = std::size_t;

    std::vector<std::complex<double>> multipole_mom_e(
        const NonRelativisticSingleScattererSolution &solution,
        const Matrix<std::complex<double>> &tau00_l,
        const lsms::math::GauntFactor &gaunt_factor,
        unsigned int lmax);


    std::complex<double> term(
        int ispin,
        int k,
        const NonRelativisticSingleScattererSolution &solution,
        const Matrix<std::complex<double>> &tau00_l,
        int kmax,
        const lsms::math::GauntFactor &gaunt_factor);


  }


}


#endif //LSMS_MULTIPOLEMOMENTS_HPP
