#ifndef LSMS_CALCULATEEVEC_HPP
#define LSMS_CALCULATEEVEC_HPP

#include <cmath>

#include "Main/SystemParameters.hpp"

void calculateEvec(LSMSSystemParameters &lsms, LocalTypeInfo &local);

void mixEvec(LSMSSystemParameters &lsms, LocalTypeInfo &local, Real alpev);

#endif
