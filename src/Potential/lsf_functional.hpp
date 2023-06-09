//
// Created by F.Moitzi on 28.08.2021.
//

#ifndef LSMS_LSF_FUNCTIONAL_HPP
#define LSMS_LSF_FUNCTIONAL_HPP

#include <map>

#include "Real.hpp"

namespace lsms {

enum class LSFFunctionalType {
  Heisenberg,
  Order1,
  Order2,
  Order3,
  Non

};

static std::map<int, LSFFunctionalType> LSFTypeMap{
    {-1, LSFFunctionalType::Heisenberg}, {1, LSFFunctionalType::Order1},
    {2, LSFFunctionalType::Order2},      {3, LSFFunctionalType::Order3},
    {0, LSFFunctionalType::Non},
};

class LSFFunctional {
 private:
  Real temperature;
  LSFFunctionalType type;

 public:
  explicit LSFFunctional(Real temperature = 0.0,
                         LSFFunctionalType type = LSFFunctionalType::Non);

  Real exchange_field(Real mag_mom);

  Real entropy(Real mag_mom);
};

template <typename T>
int signum(T val) {
  static_assert(std::is_signed<T>::value, "Signed values only!!!");
  return (T(0) < val) - (val < T(0));
}

}  // namespace lsms

#endif  // LSMS_LSF_FUNCTIONAL_HPP
