/* -*- mode: C++; c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */

#ifndef LSMS_CUCOMPLEX_OPERATORS_HPP
#define LSMS_CUCOMPLEX_OPERATORS_HPP

#include <cuComplex.h>

__device__
inline cuDoubleComplex operator+(cuDoubleComplex a, cuDoubleComplex b)
{ return cuCadd(a, b); }

__device__
inline cuDoubleComplex operator+(double a, cuDoubleComplex b)
{ return cuCadd(make_cuDoubleComplex(a, 0), b); }

__device__
inline cuDoubleComplex operator+(cuDoubleComplex a, double b)
{ return cuCadd(a, make_cuDoubleComplex(b, 0)); }

__device__
inline cuDoubleComplex operator-(cuDoubleComplex a, cuDoubleComplex b)
{ return cuCsub(a, b); }

__device__
inline cuDoubleComplex operator-(double a, cuDoubleComplex b)
{ return cuCsub(make_cuDoubleComplex(a, 0), b); }

__device__
inline cuDoubleComplex operator-(cuDoubleComplex a, double b)
{ return cuCsub(a, make_cuDoubleComplex(b, 0)); }

__device__
inline cuDoubleComplex operator*(cuDoubleComplex a, cuDoubleComplex b)
{ return cuCmul(a, b); }

__device__
inline cuDoubleComplex operator*(double a, cuDoubleComplex b)
{ return cuCmul(make_cuDoubleComplex(a, 0), b); }

__device__
inline cuDoubleComplex operator*(cuDoubleComplex a, double b)
{ return cuCmul(a, make_cuDoubleComplex(b, 0)); }

__device__
inline cuDoubleComplex operator/(cuDoubleComplex a, cuDoubleComplex b)
{ return cuCdiv(a, b); }

__device__
inline cuDoubleComplex operator/(double a, cuDoubleComplex b)
{ return cuCdiv(make_cuDoubleComplex(a, 0), b); }

__device__
 inline cuDoubleComplex operator/(cuDoubleComplex a, double b)
{ return cuCdiv(a, make_cuDoubleComplex(b, 0)); }

__device__
inline cuDoubleComplex operator-(cuDoubleComplex a)
{ return make_cuDoubleComplex(-a.x, -a.y); }

__device__
inline cuDoubleComplex exp(cuDoubleComplex a)
{ return make_cuDoubleComplex(exp(a.x) * cos(a.y), exp(a.x) * sin(a.y)); }

#endif
