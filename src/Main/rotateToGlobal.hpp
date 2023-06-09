// rotate results from the relativistic Green's function calculation
// green_function_rel.f to the global frame of reference
// this is from LSMS_1.9 gettau_c.f after the call to green_function_rel.

// this is rotr from LSMS_1.9
// ! rotation matrix for 3D vectors
// ! Simon L. Altmann: Rotations,Quaternions,..., p.75, Eq.(3.3.11)
// ! input:
// !        tvec   normal vector of axis
// !        phi    angle of rotation
// ! output:
// !        drot   matrix of rotation
// !

#include "Array3d.hpp"
#include "Complex.hpp"
#include "Matrix.hpp"
#include "SingleSite/AtomData.hpp"

#ifndef LSMS_ROTATE_TO_GLOBAL_HPP
#define LSMS_ROTATE_TO_GLOBAL_HPP

void rotateToGlobal(AtomData &atom, Matrix<Complex> &dos,
                    Matrix<Complex> &dosck, Matrix<Complex> &dos_orb,
                    Matrix<Complex> &dosck_rob, Array3d<Complex> &green,
                    Array3d<Complex> &dens_orb, int i);

#endif  // LSMS_ROTATE_TO_GLOBAL_HPP
