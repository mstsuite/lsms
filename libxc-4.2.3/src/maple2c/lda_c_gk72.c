/* 
  This file was generated automatically with ../scripts/maple2c.pl.
  Do not edit this file directly as it can be overwritten!!

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.

  Maple version     : Maple 2016 (X86 64 LINUX)
  Maple source      : ..//maple/lda_c_gk72.mpl
  Type of functional: work_lda
*/

static void
func0(const xc_func_type *p, xc_lda_work_t *r)
{
  double t2, t3, t5, t6, t7, t11, t12, t17;
  double t18, t20, t22, t25, t26, t27, t30, t31;
  double t37, t39, t40, t46, t54, t59, t60, t61;
  double t63, t66, t67, t69, t74, t75, t82, t91;
  double t93, t95, t96, t97, t100, t101, t103, t105;
  double t111, t115, t117, t119, t121, t123, t124, t144;
  double t145, t165, t182;


  t2 = log(r->rs);
  t3 = r->rs * t2;
  t5 = -0.7e0 + r->rs;
  t6 = Heaviside(t5);
  t7 = t2 * t6;
  t11 = -0.10e2 + r->rs;
  t12 = Heaviside(t11);
  t17 = 0.1e1 / r->rs;
  t18 = t17 * t12;
  t20 = sqrt(r->rs);
  t22 = 0.1e1 / t20 / r->rs;
  t25 = r->rs * r->rs;
  t26 = 0.1e1 / t25;
  t27 = t26 * t12;
  t30 = 0.1e1 / t20 / t25;
  t31 = t30 * t12;
  r->f = -0.48e-1 - 0.17e-1 * r->rs + 0.9e-2 * t3 - 0.1212e-1 * t7 + 0.17e-1 * r->rs * t6 - 0.1898e-1 * t2 * t12 - 0.9e-2 * t3 * t6 + 0.438e0 * t18 + 0.1325e1 * t22 * t12 - 0.147e1 * t27 - 0.4e0 * t31 + 0.311e-1 * t2 - 0.1356e-1 * t6 + 0.6156e-1 * t12;

  if(r->order < 1) return;

  t37 = t17 * t6;
  t39 = 0.0;
  t40 = t2 * t39;
  t46 = 0.0;
  t54 = t17 * t46;
  t59 = t25 * r->rs;
  t60 = 0.1e1 / t59;
  t61 = t60 * t12;
  t63 = t26 * t46;
  t66 = 0.1e1 / t20 / t59;
  t67 = t66 * t12;
  t69 = t30 * t46;
  t74 = -0.438e0 * t27 + 0.438e0 * t54 - 0.19875000000000000000e1 * t31 + 0.1325e1 * t22 * t46 + 0.294e1 * t61 - 0.147e1 * t63 + 0.10000000000000000000e1 * t67 - 0.4e0 * t69 + 0.311e-1 * t17 - 0.1356e-1 * t39 + 0.6156e-1 * t46;
  r->dfdrs = t74 - 0.8e-2 + 0.9e-2 * t2 - 0.1212e-1 * t37 - 0.1212e-1 * t40 + 0.8e-2 * t6 + 0.17e-1 * r->rs * t39 - 0.1898e-1 * t18 - 0.1898e-1 * t2 * t46 - 0.9e-2 * t7 - 0.9e-2 * t3 * t39;

  if(r->order < 2) return;

  t75 = 0.0;
  t82 = 0.0;
  t91 = -0.9e-2 * t3 * t75 + 0.9e-2 * t17 - 0.311e-1 * t26 + 0.16e-1 * t39 - 0.1356e-1 * t75 + 0.6156e-1 * t82 + 0.1898e-1 * t27 - 0.9e-2 * t37 - 0.18e-1 * t40 - 0.3796e-1 * t54 + 0.876e0 * t61 - 0.876e0 * t63 + 0.49687500000000000000e1 * t67;
  t93 = t26 * t6;
  t95 = t25 * t25;
  t96 = 0.1e1 / t95;
  t97 = t96 * t12;
  t100 = 0.1e1 / t20 / t95;
  t101 = t100 * t12;
  t103 = t17 * t39;
  t105 = t2 * t75;
  t111 = t17 * t82;
  t115 = t60 * t46;
  t117 = t26 * t82;
  t119 = t66 * t46;
  t121 = t30 * t82;
  t123 = -0.39750000000000000000e1 * t69 + 0.1212e-1 * t93 - 0.882e1 * t97 - 0.35000000000000000000e1 * t101 - 0.2424e-1 * t103 - 0.1212e-1 * t105 + 0.17e-1 * r->rs * t75 - 0.1898e-1 * t2 * t82 + 0.438e0 * t111 + 0.1325e1 * t22 * t82 + 0.588e1 * t115 - 0.147e1 * t117 + 0.20000000000000000000e1 * t119 - 0.4e0 * t121;
  r->d2fdrs2 = t123 + t91;

  if(r->order < 3) return;

  t124 = 0.0;
  t144 = -0.9e-2 * t3 * t124 - 0.9e-2 * t26 + 0.622e-1 * t60 + 0.24e-1 * t75 - 0.3796e-1 * t61 + 0.5694e-1 * t63 + 0.9e-2 * t93 - 0.2628e1 * t97 - 0.17390625000000000000e2 * t101 - 0.27e-1 * t103 - 0.27e-1 * t105 - 0.5694e-1 * t111 + 0.2628e1 * t115 - 0.1314e1 * t117 + 0.14906250000000000000e2 * t119 - 0.59625000000000000000e1 * t121 - 0.2424e-1 * t60 * t6;
  t145 = t95 * r->rs;
  t165 = 0.0;
  t182 = 0.3528e2 / t145 * t12 + 0.15750000000000000000e2 / t20 / t145 * t12 + 0.3636e-1 * t26 * t39 - 0.2646e2 * t96 * t46 - 0.10500000000000000000e2 * t100 * t46 - 0.3636e-1 * t17 * t75 - 0.1212e-1 * t2 * t124 + 0.17e-1 * r->rs * t124 - 0.1898e-1 * t2 * t165 + 0.438e0 * t17 * t165 + 0.1325e1 * t22 * t165 + 0.882e1 * t60 * t82 - 0.147e1 * t26 * t165 + 0.30000000000000000000e1 * t66 * t82 - 0.4e0 * t30 * t165 - 0.1356e-1 * t124 + 0.6156e-1 * t165;
  r->d3fdrs3 = t182 + t144;

  if(r->order < 4) return;


}

static void
func1(const xc_func_type *p, xc_lda_work_t *r)
{
  double t2, t3, t5, t6, t7, t11, t12, t17;
  double t18, t20, t22, t25, t26, t27, t30, t31;
  double t37, t39, t40, t46, t54, t59, t60, t61;
  double t63, t66, t67, t69, t74, t75, t77, t91;
  double t93, t95, t96, t97, t100, t101, t103, t105;
  double t111, t115, t117, t119, t121, t123, t124, t126;
  double t144, t148, t182;


  t2 = log(r->rs);
  t3 = r->rs * t2;
  t5 = -0.7e0 + r->rs;
  t6 = Heaviside(t5);
  t7 = t2 * t6;
  t11 = -0.10e2 + r->rs;
  t12 = Heaviside(t11);
  t17 = 0.1e1 / r->rs;
  t18 = t17 * t12;
  t20 = sqrt(r->rs);
  t22 = 0.1e1 / t20 / r->rs;
  t25 = r->rs * r->rs;
  t26 = 0.1e1 / t25;
  t27 = t26 * t12;
  t30 = 0.1e1 / t20 / t25;
  t31 = t30 * t12;
  r->f = -0.48e-1 - 0.17e-1 * r->rs + 0.9e-2 * t3 - 0.1212e-1 * t7 + 0.17e-1 * r->rs * t6 - 0.1898e-1 * t2 * t12 - 0.9e-2 * t3 * t6 + 0.438e0 * t18 + 0.1325e1 * t22 * t12 - 0.147e1 * t27 - 0.4e0 * t31 + 0.311e-1 * t2 - 0.1356e-1 * t6 + 0.6156e-1 * t12;

  if(r->order < 1) return;

  t37 = t17 * t6;
  t39 = 0.0;
  t40 = t2 * t39;
  t46 = 0.0;
  t54 = t17 * t46;
  t59 = t25 * r->rs;
  t60 = 0.1e1 / t59;
  t61 = t60 * t12;
  t63 = t26 * t46;
  t66 = 0.1e1 / t20 / t59;
  t67 = t66 * t12;
  t69 = t30 * t46;
  t74 = -0.438e0 * t27 + 0.438e0 * t54 - 0.19875000000000000000e1 * t31 + 0.1325e1 * t22 * t46 + 0.294e1 * t61 - 0.147e1 * t63 + 0.10000000000000000000e1 * t67 - 0.4e0 * t69 + 0.311e-1 * t17 - 0.1356e-1 * t39 + 0.6156e-1 * t46;
  r->dfdrs = t74 - 0.8e-2 + 0.9e-2 * t2 - 0.1212e-1 * t37 - 0.1212e-1 * t40 + 0.8e-2 * t6 + 0.17e-1 * r->rs * t39 - 0.1898e-1 * t18 - 0.1898e-1 * t2 * t46 - 0.9e-2 * t7 - 0.9e-2 * t3 * t39;
  r->dfdz = 0.0e0;

  if(r->order < 2) return;

  t75 = 0.0;
  t77 = 0.0;
  t91 = 0.6156e-1 * t75 - 0.1356e-1 * t77 + 0.16e-1 * t39 + 0.9e-2 * t17 - 0.311e-1 * t26 - 0.9e-2 * t3 * t77 + 0.1898e-1 * t27 - 0.9e-2 * t37 - 0.18e-1 * t40 - 0.3796e-1 * t54 + 0.876e0 * t61 - 0.876e0 * t63 + 0.49687500000000000000e1 * t67;
  t93 = t26 * t6;
  t95 = t25 * t25;
  t96 = 0.1e1 / t95;
  t97 = t96 * t12;
  t100 = 0.1e1 / t20 / t95;
  t101 = t100 * t12;
  t103 = t17 * t39;
  t105 = t2 * t77;
  t111 = t17 * t75;
  t115 = t60 * t46;
  t117 = t26 * t75;
  t119 = t66 * t46;
  t121 = t30 * t75;
  t123 = -0.39750000000000000000e1 * t69 + 0.1212e-1 * t93 - 0.882e1 * t97 - 0.35000000000000000000e1 * t101 - 0.2424e-1 * t103 - 0.1212e-1 * t105 + 0.17e-1 * r->rs * t77 - 0.1898e-1 * t2 * t75 + 0.438e0 * t111 + 0.1325e1 * t22 * t75 + 0.588e1 * t115 - 0.147e1 * t117 + 0.20000000000000000000e1 * t119 - 0.4e0 * t121;
  r->d2fdrs2 = t123 + t91;
  r->d2fdrsz = 0.0e0;
  r->d2fdz2 = 0.0e0;

  if(r->order < 3) return;

  t124 = 0.0;
  t126 = 0.0;
  t144 = 0.6156e-1 * t124 - 0.1356e-1 * t126 + 0.24e-1 * t77 + 0.622e-1 * t60 - 0.9e-2 * t26 - 0.9e-2 * t3 * t126 - 0.3796e-1 * t61 + 0.5694e-1 * t63 + 0.9e-2 * t93 - 0.2628e1 * t97 - 0.17390625000000000000e2 * t101 - 0.27e-1 * t103 - 0.27e-1 * t105 - 0.5694e-1 * t111 + 0.2628e1 * t115 - 0.1314e1 * t117 + 0.14906250000000000000e2 * t119;
  t148 = t95 * r->rs;
  t182 = -0.59625000000000000000e1 * t121 - 0.2424e-1 * t60 * t6 + 0.3528e2 / t148 * t12 + 0.15750000000000000000e2 / t20 / t148 * t12 + 0.3636e-1 * t26 * t39 - 0.2646e2 * t96 * t46 - 0.10500000000000000000e2 * t100 * t46 - 0.3636e-1 * t17 * t77 - 0.1212e-1 * t2 * t126 + 0.17e-1 * r->rs * t126 - 0.1898e-1 * t2 * t124 + 0.438e0 * t17 * t124 + 0.1325e1 * t22 * t124 + 0.882e1 * t60 * t75 - 0.147e1 * t26 * t124 + 0.30000000000000000000e1 * t66 * t75 - 0.4e0 * t30 * t124;
  r->d3fdrs3 = t182 + t144;
  r->d3fdrs2z = 0.0e0;
  r->d3fdrsz2 = 0.0e0;
  r->d3fdz3 = 0.0e0;

  if(r->order < 4) return;


}

void 
xc_lda_c_gk72_func(const xc_func_type *p, xc_lda_work_t *r)
{
  if(p->nspin == XC_UNPOLARIZED)
    func0(p, r);
  else
    func1(p, r);
}

#define maple2c_order 3
#define maple2c_func  xc_lda_c_gk72_func
