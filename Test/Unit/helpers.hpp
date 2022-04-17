//
// Created by F.Moitzi on 17.04.2022.
//

#ifndef LSMS_HELPERS_HPP
#define LSMS_HELPERS_HPP

#include <vector>

namespace lsms_helper {

template <typename T>
std::vector<T> linspace(double start, double end, std::size_t num) {
  std::vector<T> linspaced;

  if (0 != num) {
    if (1 == num) {
      linspaced.push_back(static_cast<T>(start));
    } else {
      double delta = (end - start) / (num - 1);

      for (auto i = 0; i < (num - 1); ++i) {
        linspaced.push_back(static_cast<T>(start + delta * i));
      }
      // ensure that start and end are exactly the same as the input
      linspaced.push_back(static_cast<T>(end));
    }
  }
  return linspaced;
}

}  // namespace lsms_helper

#endif  // LSMS_HELPERS_HPP
