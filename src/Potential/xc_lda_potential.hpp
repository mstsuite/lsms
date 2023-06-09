
#ifndef RADIALDFTLIB_XCPOTENTIALS_H
#define RADIALDFTLIB_XCPOTENTIALS_H

#include <cmath>
#include <vector>

namespace lsms {

template<typename T>
inline T get_Y(T y, T b, T c) {
  return y * y + b * y + c;
}

template<typename T>
inline void VxcLDA(const T &Rgrid,
                   const T &density,
                   T &Vxc, T &exc) {

  constexpr T y0 = -0.10498;
  constexpr T b = 3.72744;
  constexpr T c = 12.9352;
  constexpr T A = 0.0621814;

  T Q, rs, y, ec, ex, Vc, Vx, beta, mu, R, S;

  if (density < 1.0e-14) {
    Vxc = 0;
    exc = 0;
    return;
  }

  /*
   * Correlation:
   */

  Q = sqrt(4 * c - b * b);
  rs = cbrt(3.0 / (4.0 * M_PI * density));
  y = sqrt(rs);

  ec = A / 2.0 * (log(y * y / get_Y(y, b, c)) + 2 * b / Q * atan(Q / (2 * y + b))
      - b * y0 / get_Y(y0, b, c) *
          (log((y - y0) * (y - y0) / get_Y(y, b, c))
              + 2 * (b + 2 * y0) / Q * atan(Q / (2 * y + b))
          ));

  Vc = ec - A / 6.0 * (c * (y - y0) - b * y0 * y) / ((y - y0) * get_Y(y, b, c));

  /*
   * Exchange:
   */

  // exchange energy density
  ex = -3.0 / (4.0 * M_PI) * cbrt(3 * M_PI * M_PI * density);
  // exchange potential
  Vx = 4.0 * ex / 3.0;

  /*
   *
   */

  exc = ex + ec;
  Vxc = Vx + Vc;

  Vxc = Vxc * 2.0;
  exc = exc * 2.0;

}

template<typename T>
inline void VxcRLDA(const T &Rgrid,
                    const T &density,
                    T &Vxc, T &exc) {

  constexpr T c_light = 137.03599908421;
  constexpr T y0 = -0.10498;
  constexpr T b = 3.72744;
  constexpr T c = 12.9352;
  constexpr T A = 0.0621814;

  T Q, rs, y, ec, ex, Vc, Vx, beta, mu, R, S;

  if (density < 1.0e-14) {
    Vxc = 0;
    exc = 0;
    return;
  }

  /*
   * Correlation:
   */

  Q = sqrt(4 * c - b * b);
  rs = cbrt(3.0 / (4.0 * M_PI * density));
  y = sqrt(rs);

  ec = A / 2.0 * (log(y * y / get_Y(y, b, c)) + 2 * b / Q * atan(Q / (2 * y + b))
      - b * y0 / get_Y(y0, b, c) *
          (log((y - y0) * (y - y0) / get_Y(y, b, c))
              + 2 * (b + 2 * y0) / Q * atan(Q / (2 * y + b))
          ));

  Vc = ec - A / 6.0 * (c * (y - y0) - b * y0 * y) / ((y - y0) * get_Y(y, b, c));

  /*
   * Exchange:
   */

  // exchange energy density
  ex = -3.0 / (4.0 * M_PI) * cbrt(3 * M_PI * M_PI * density);
  // exchange potential
  Vx = 4.0 * ex / 3.0;

  /*
   *
   */

  beta = -4 * M_PI * ex / (3 * c_light);
  mu = sqrt(1 + beta * beta);
  R = 1 - 3 * pow(((beta * mu - log(beta + mu)) / (beta * beta)), 2) / 2.0;
  S = 3 * log(beta + mu) / (2 * beta * mu) - 1.0 / 2.0;

  ex = ex * R;
  Vx = Vx * S;
  exc = ex + ec;
  Vxc = Vx + Vc;

  Vxc = Vxc * 2.0;
  exc = exc * 2.0;

};

}

#endif //RADIALDFTLIB_XCPOTENTIALS_H
