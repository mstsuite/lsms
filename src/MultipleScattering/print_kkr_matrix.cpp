//
// Created by F.Moitzi on 15.03.2023.
//

#include <complex>
#include <fstream>

#include "print_kkr_matrix.hpp"

#include "fmt/core.h"

namespace lsms {

void print_kkr_matrix(FILE *file,
                      AtomData &atom,
                      const Matrix<Complex> &m) {

  for (int ir1 = 0; ir1 < atom.numLIZ; ir1++) {
    int kkr1 = (atom.LIZlmax[ir1] + 1) * (atom.LIZlmax[ir1] + 1);

    for (int ir2 = 0; ir2 < atom.numLIZ; ir2++) {

      int kkr2 = (atom.LIZlmax[ir2] + 1) * (atom.LIZlmax[ir2] + 1);

      auto x = atom.LIZPos(0, ir1) - atom.LIZPos(0, ir2);
      auto y = atom.LIZPos(1, ir1) - atom.LIZPos(1, ir2);
      auto z = atom.LIZPos(2, ir1) - atom.LIZPos(2, ir2);

//      for (auto i = 0; i < kkr1; i++) {
//        for (auto j = 0; j < kkr2; j++) {
      fmt::print(file, "{} {} {:6.3f} {:6.3f} {:6.3f}\n", ir1, ir2, x, y, z);
//        }
//      }

    }

  }

}

} // lsms