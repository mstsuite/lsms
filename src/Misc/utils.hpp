
#ifndef LSMS_UTILS_HPP
#define LSMS_UTILS_HPP

#include <cmath>
#include <type_traits>
#include <vector>

namespace lsms {

/**
 * Calculate euclidian norm
 */
template <typename T>
T static inline norm(std::vector<T> distance) {
  return std::sqrt(distance[0] * distance[0] + distance[1] * distance[1] +
                   distance[2] * distance[2]);
}

/**
 * Calculate euclidian norm
 */
template <typename T>
T static inline norm(T *distance) {
  return std::sqrt(distance[0] * distance[0] + distance[1] * distance[1] +
                   distance[2] * distance[2]);
}

/**
 * Exponential factor
 */
template <typename T, typename = typename std::enable_if<
                          std::is_arithmetic<T>::value, T>::type>
inline constexpr T pow(const T base, unsigned const exponent) {
  return (exponent == 0) ? 1 : (base * pow(base, exponent - 1));
}
}  // namespace lsms

#endif  // LSMS_UTILS_HPP