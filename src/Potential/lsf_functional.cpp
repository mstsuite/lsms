//
// Created by F.Moitzi on 28.08.2021.
//

#include "lsf_functional.hpp"

#include "PhysicalConstants.hpp"
#include "Real.hpp"

lsms::LSFFunctional::LSFFunctional(Real temperature, LSFFunctionalType type)
    : temperature{temperature * convertKtoRydberg}, type{type} {}

Real lsms::LSFFunctional::exchange_field(Real mag_mom) {
  Real field = lsms::signum(mag_mom) * temperature;

  if (type == LSFFunctionalType::Heisenberg) {
    field *= 1.0 / (1.0 + std::abs(mag_mom));
  } else if (type == LSFFunctionalType::Order1) {
    field *= 1.0 / (0.01 + std::abs(mag_mom));
  } else if (type == LSFFunctionalType::Order2) {
    field *= 2.0 / (0.01 + std::abs(mag_mom));
  } else if (type == LSFFunctionalType::Order3) {
    field *= 3.0 / (0.01 + std::abs(mag_mom));
  } else {
    field = 0.0;
  }

  return field;
}

Real lsms::LSFFunctional::entropy(Real mag_mom) {
  if (type == LSFFunctionalType::Heisenberg) {
    return 1.0 * std::log(1.0 + std::abs(mag_mom));
  } else if (type == LSFFunctionalType::Order1) {
    return 1.0 * std::log(0.01 + std::abs(mag_mom));
  } else if (type == LSFFunctionalType::Order2) {
    return 2.0 * std::log(0.01 + std::abs(mag_mom));
  } else if (type == LSFFunctionalType::Order3) {
    return 3.0 * std::log(0.01 + std::abs(mag_mom));
  }

  return 0.0;
}
