//
// Created by F.Moitzi on 24.12.2021.
//

#ifndef MADELUNG_INTEGER_FACTORS_HPP
#define MADELUNG_INTEGER_FACTORS_HPP

#include <type_traits>
#include <vector>

namespace lsms {

std::vector<int> get_lofk(int l0);

std::vector<int> get_mofk(int l0);

std::vector<int> get_lofj(int l0);

std::vector<int> get_mofj(int l0);

std::vector<int> get_kofj(int l0);

template <typename Integer,
          typename = std::enable_if_t<std::is_integral<Integer>::value>>
Integer get_jmax(Integer lmax) {
  return (lmax + 1) * (lmax + 2) / 2;
}

template <typename Integer,
          typename = std::enable_if_t<std::is_integral<Integer>::value>>
Integer get_kmax(Integer lmax) {
  return (lmax + 1) * (lmax + 1);
}

}  // namespace lsms

#endif  // MADELUNG_INTEGER_FACTORS_HPP
