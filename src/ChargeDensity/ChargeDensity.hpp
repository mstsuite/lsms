//
// Created by F.Moitzi on 03.01.2023.
//

#ifndef LSMS_SRC_CHARGEDENSITY_CHARGEDENSITY_HPP_
#define LSMS_SRC_CHARGEDENSITY_CHARGEDENSITY_HPP_

#include "AtomData.hpp"
#include "SystemParameters.hpp"

namespace lsms {

/**
 * Calculated the local radial charge density
 *
 * ONLY ASA version for nspin == 1 and nspin == 2
 *
 */
void calculateRadialChargeDensity(LSMSSystemParameters &lsms,
                                  LocalTypeInfo &local);

/**
 * Calculate integrate local charges
 *
 * ONLY ASA version for nspin = 1 and nspin = 2
 */
void calculateCharge(LSMSSystemParameters &lsms, LocalTypeInfo &local,
                     std::vector<Real> &qsub);

/**
 * Check radial charge density
 */
void checkRadialChargeDensity(LSMSSystemParameters &lsms,
                              LocalTypeInfo &local);

/**
 *
 */
void copyChargesAndPotential(LSMSSystemParameters &lsms,
                             LocalTypeInfo &local);

/**
 * Calculate RMS
 */
void calculateQRMS(LSMSSystemParameters &lsms, LocalTypeInfo &local);

}  // namespace lsms

#endif  // LSMS_SRC_CHARGEDENSITY_CHARGEDENSITY_HPP_
