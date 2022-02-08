//
// Created by F.Moitzi on 24.01.2022.
//

#ifndef LSMS_MULTIPOLEMOMENTS_HPP
#define LSMS_MULTIPOLEMOMENTS_HPP

#include <complex>

#include "SingleSite/SingleSiteScattering.hpp"

#include "Array3d.hpp"
#include "utils.hpp"
#include "integrator.hpp"
#include "NDArray.hpp"
#include "integer_factors.hpp"
#include "common.hpp"
#include "GauntFactor.hpp"

namespace lsms {

  namespace multi_moms {

    using size_t = std::size_t;

    /**
     * Energy-resolved multipole moments
     */
    std::vector<std::complex<double>> multipole_mom_e(
        const NonRelativisticSingleScattererSolution &solution,
        const Matrix<std::complex<double>> &tau00_l,
        const lsms::math::GauntFactor &gaunt_factor,
        unsigned int lmax_mm,
        unsigned int lmax);


    /**
     * Integrated multipole moments
     */
    void integrateMultipoleMoments(LSMSSystemParameters &lsms,
                                   int is,
                                   int ie,
                                   int nume,
                                   int lmax_mom,
                                   std::complex<double> energy,
                                   std::complex<double> dele1,
                                   std::vector<std::complex<double>> &mm_e,
                                   AtomData &atom);





  }


}


#endif //LSMS_MULTIPOLEMOMENTS_HPP
