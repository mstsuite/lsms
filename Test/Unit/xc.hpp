//
// Created by F.Moitzi on 24.04.2022.
//

#ifndef LSMS_XC_HPP
#define LSMS_XC_HPP

extern "C" {

void getvxc_scalar(double *n, bool *relat, double *c_light, double *exc,
                   double *vxc);
}

#endif  // LSMS_XC_HPP
