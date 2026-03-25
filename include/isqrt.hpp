// provide isqrt until it becomes part of the C++ standard
// requires unsigned type

#ifndef LSMS_ISQRT_HPP
#define LSMS_ISQRT_HPP

#include <bit>

template<class T>
constexpr T isqrt(const T n) noexcept {
if (n <= T{1})
  return n;

T i_current{0}, i_next{T(T{1} << ((std::bit_width(T(n - 1)) + 1) >> 1))};
do {
  i_current = i_next;
  i_next = T((i_current + n / i_current) >> 1);
} while (i_next < i_current);

return i_current;
}

#endif
