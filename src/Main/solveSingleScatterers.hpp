
#include <vector>
#include <cmath>

#include "Real.hpp"
#include "Complex.hpp"
#include "PhysicalConstants.hpp"

#include "SingleSite/SingleSiteScattering.hpp"
#include "Misc/Indices.hpp"
#include "Misc/Coeficients.hpp"

#ifndef LSMS_SOLVE_SINGLE_SCATTERERS
#define LSMS_SOLVE_SINGLE_SCATTERERS

void solveSingleScatterers(LSMSSystemParameters &lsms,
                           LocalTypeInfo &local,
                           std::vector<Matrix<Real>> &vr,
                           Complex energy,
                           std::vector<NonRelativisticSingleScattererSolution> &solution,
                           int iie);

void solveSingleScatterers(LSMSSystemParameters &lsms,
                           LocalTypeInfo &local,
                           std::vector<Matrix<Real>> &vr,
                           Complex energy,
                           std::vector<RelativisticSingleScattererSolution> &solution,
                           int iie);

#endif // LSMS_SOLVE_SINGLE_SCATTERERS