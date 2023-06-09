//
// Created by F.Moitzi on 15.03.2023.
//

#ifndef LSMS_SRC_MULTIPLESCATTERING_PRINT_KKR_MATRIX_HPP_
#define LSMS_SRC_MULTIPLESCATTERING_PRINT_KKR_MATRIX_HPP_

#include "Complex.hpp"
#include "Matrix.hpp"
#include "AtomData.hpp"

namespace lsms {

void print_kkr_matrix(FILE *file,
                      AtomData &atom,
                      const Matrix<Complex> &m);

} // lsms

#endif //LSMS_SRC_MULTIPLESCATTERING_PRINT_KKR_MATRIX_HPP_
