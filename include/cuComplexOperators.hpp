/* -*- mode: C++; c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */

#ifndef LSMS_CUCOMPLEX_OPERATORS_HPP
#define LSMS_CUCOMPLEX_OPERATORS_HPP

#include <cuComplex.h>

inline cuDoubleComplex operator+(cuDoubleComplex a, cuDoubleComplex b)
{ return cuCadd(a, b); }

inline cuDoubleComplex operator+(double a, cuDoubleComplex b)
{ return cuCadd(make_cuDoubleComplex(a, 0), b); }

inline cuDoubleComplex operator+(cuDoubleComplex a, double b)
{ return cuCadd(a, make_cuDoubleComplex(b, 0)); }

inline cuDoubleComplex operator-(cuDoubleComplex a, cuDoubleComplex b)
{ return cuCsub(a, b); }

inline cuDoubleComplex operator-(double a, cuDoubleComplex b)
{ return cuCsub(make_cuDoubleComplex(a, 0), b); }

inline cuDoubleComplex operator-(cuDoubleComplex a, double b)
{ return cuCsub(a, make_cuDoubleComplex(b, 0)); }

inline cuDoubleComplex operator*(cuDoubleComplex a, cuDoubleComplex b)
{ return cuCmul(a, b); }

inline cuDoubleComplex operator*(double a, cuDoubleComplex b)
{ return cuCmul(make_cuDoubleComplex(a, 0), b); }

inline cuDoubleComplex operator*(cuDoubleComplex a, double b)
{ return cuCmul(a, make_cuDoubleComplex(b, 0)); }

inline cuDoubleComplex operator/(cuDoubleComplex a, cuDoubleComplex b)
{ return cuCdiv(a, b); }

inline cuDoubleComplex operator/(double a, cuDoubleComplex b)
{ return cuCdiv(make_cuDoubleComplex(a, 0), b); }

 inline cuDoubleComplex operator/(cuDoubleComplex a, double b)
{ return cuCdiv(a, make_cuDoubleComplex(b, 0)); }

inline cuDoubleComplex operator-(cuDoubleComplex a)
{ return make_cuDoubleComplex(-a.x, -a.y); }

inline cuDoubleComplex exp(cuDoubleComplex a)
{ return make_cuDoubleComplex(a.x * cos(a.y), a.x * sin(a.y)); }

#endif
