/*
  This file was generated automatically with ./scripts/maple2c.pl.
  Do not edit this file directly as it can be overwritten!!

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.

  Maple version     : Maple 2020 (X86 64 LINUX)
  Maple source      : ./maple/lda_exc/lda_x_1d_soft.mpl
  Type of functional: lda_exc
*/

#define maple2c_order 4
#define MAPLE2C_FLAGS (XC_FLAGS_I_HAVE_EXC | XC_FLAGS_I_HAVE_VXC | XC_FLAGS_I_HAVE_FXC | XC_FLAGS_I_HAVE_KXC | XC_FLAGS_I_HAVE_LXC)


static inline void
func_unpol(const xc_func_type *p, int order, const double *rho , double *zk LDA_OUT_PARAMS_NO_EXC(XC_COMMA double *, ))
{

#ifndef XC_DONT_COMPILE_EXC
  double t3, t4, t5, t7, t8, t11, t12, t14;
  double t15, t16, t17, t18, t19, t24;

#ifndef XC_DONT_COMPILE_VXC
  double t25, t26, t27, t28, t32;

#ifndef XC_DONT_COMPILE_FXC
  double t36, t37, t38, t42, t47;

#ifndef XC_DONT_COMPILE_KXC
  double t52, t53, t54, t60, t66;

#ifndef XC_DONT_COMPILE_LXC
  double t70, t77, t93;
#endif

#endif

#endif

#endif

#endif


  lda_x_1d_exponential_params *params;

  assert(p->params != NULL);
  params = (lda_x_1d_exponential_params * )(p->params);

  t3 = 0.1e1 <= p->zeta_threshold;
  t4 = rho[0] / 0.2e1 <= p->dens_threshold || t3;
  t5 = p->zeta_threshold - 0.1e1;
  t7 = my_piecewise5(t3, t5, t3, -t5, 0);
  t8 = 0.1e1 + t7;
  t11 = t8 * M_PI * params->beta * rho[0];
  t12 = xc_integrate(func1, NULL, 0.0, t11);
  t14 = xc_integrate(func2, NULL, 0.0, t11);
  t15 = 0.1e1 / M_PI;
  t16 = t14 * t15;
  t17 = 0.1e1 / params->beta;
  t18 = 0.1e1 / rho[0];
  t19 = t17 * t18;
  t24 = my_piecewise3(t4, 0, -0.79577471545947667883e-1 * (t8 * t12 - t16 * t19) * t17);
  if(zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    zk[0] = 0.2e1 * t24;

#ifndef XC_DONT_COMPILE_VXC

  if(order < 1) return;


  t25 = params->beta * params->beta;
  t26 = 0.1e1 / t25;
  t27 = rho[0] * rho[0];
  t28 = 0.1e1 / t27;
  t32 = my_piecewise3(t4, 0, -0.79577471545947667883e-1 * t16 * t26 * t28);
  if(vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    vrho[0] = 0.2e1 * rho[0] * t32 + 0.2e1 * t24;

#ifndef XC_DONT_COMPILE_FXC

  if(order < 2) return;


  t36 = t8 * t8;
  t37 = xc_bessel_K0( t11);
  t38 = t36 * t37;
  t42 = 0.1e1 / t27 / rho[0];
  t47 = my_piecewise3(t4, 0, -0.50000000000000000000e0 * t38 * t18 + 0.15915494309189533577e0 * t16 * t26 * t42);
  if(v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    v2rho2[0] = 0.2e1 * rho[0] * t47 + 0.4e1 * t32;

#ifndef XC_DONT_COMPILE_KXC

  if(order < 3) return;


  t52 = xc_bessel_K1( t11);
  t53 = t36 * t8 * t52;
  t54 = M_PI * params->beta;
  t60 = t27 * t27;
  t66 = my_piecewise3(t4, 0, 0.50000000000000000000e0 * t53 * t54 * t18 + 0.15000000000000000000e1 * t38 * t28 - 0.47746482927568600731e0 * t16 * t26 / t60);
  if(v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    v3rho3[0] = 0.2e1 * rho[0] * t66 + 0.6e1 * t47;

#ifndef XC_DONT_COMPILE_LXC

  if(order < 4) return;


  t70 = t36 * t36;
  t77 = M_PI * M_PI;
  t93 = my_piecewise3(t4, 0, 0.50000000000000000000e0 * t70 * (-t37 - 0.1e1 / t8 * t15 * t19 * t52) * t77 * t25 * t18 - 0.20000000000000000000e1 * t53 * t54 * t28 - 0.60000000000000000001e1 * t38 * t42 + 0.19098593171027440292e1 * t16 * t26 / t60 / rho[0]);
  if(v4rho4 != NULL && (p->info->flags & XC_FLAGS_HAVE_LXC))
    v4rho4[0] = 0.2e1 * rho[0] * t93 + 0.8e1 * t66;

#ifndef XC_DONT_COMPILE_MXC

  if(order < 5) return;


#endif

#endif

#endif

#endif

#endif


}


static inline void
func_pol(const xc_func_type *p, int order, const double *rho , double *zk LDA_OUT_PARAMS_NO_EXC(XC_COMMA double *, ))
{

#ifndef XC_DONT_COMPILE_EXC
  double t2, t3, t4, t5, t7, t8, t9, t11;
  double t12, t13, t14, t15, t16, t17, t18, t20;
  double t21, t22, t23, t24, t29, t31, t32, t33;
  double t34, t35, t36, t38, t39, t44;

#ifndef XC_DONT_COMPILE_VXC
  double t45, t46, t47, t48, t49, t51, t52, t56;
  double t58, t60, t64, t67, t68, t73, t75, t80;

#ifndef XC_DONT_COMPILE_FXC
  double t86, t87, t89, t90, t92, t94, t95, t96;
  double t97, t100, t101, t103, t104, t106, t110, t112;
  double t114, t116, t117, t118, t119, t122, t123, t125;
  double t127, t131, t134, t135, t137, t143, t144, t146;
  double t152, t158, t159, t161, t163, t164, t167, t173;
  double t175, t177, t179, t180, t183, t189;

#ifndef XC_DONT_COMPILE_KXC
  double t194, t195, t196, t198, t199, t204, t206, t208;
  double t209, t212, t214, t217, t219, t220, t222, t223;
  double t225, t226, t227, t228, t229, t231, t235, t237;
  double t242, t244, t246, t247, t250, t252, t255, t257;
  double t258, t260, t261, t263, t264, t265, t266, t268;
  double t272, t275, t276, t277, t278, t279, t280, t285;
  double t294, t296, t301, t310, t313, t314, t319, t322;
  double t324, t325, t326, t329, t332, t335, t336, t337;
  double t342, t345, t348, t350, t355, t358, t360, t361;
  double t362, t365, t368, t371, t372, t373, t378, t381;
  double t384, t390, t391, t393, t396, t399, t400, t403;
  double t404, t407, t410, t413, t420, t422, t424, t427;
  double t430, t431, t434, t435, t438, t441, t444, t451;

#ifndef XC_DONT_COMPILE_LXC
  double t457, t458, t460, t461, t463, t464, t465, t467;
  double t468, t469, t470, t475, t478, t479, t480, t481;
  double t482, t483, t484, t485, t487, t499, t501, t512;
  double t514, t516, t517, t519, t520, t523, t525, t527;
  double t528, t529, t531, t532, t533, t534, t539, t542;
  double t543, t544, t545, t546, t547, t548, t550, t562;
  double t564, t575, t577, t579, t580, t582, t583, t586;
  double t608, t609, t610, t611, t613, t616, t635, t637;
  double t640, t665, t666, t674, t676, t678, t682, t694;
  double t705, t715, t717, t721, t740, t741, t752, t754;
  double t755, t757, t772, t783, t790, t792, t796, t815;
  double t816, t827, t829, t832, t834, t857, t869, t870;
  double t872, t876, t878, t894, t895, t906, t908, t910;
  double t912, t936, t948, t952, t965, t970, t971, t981;
  double t1004, t1007, t1016, t1021, t1031, t1054, t1057;
#endif

#endif

#endif

#endif

#endif


  lda_x_1d_exponential_params *params;

  assert(p->params != NULL);
  params = (lda_x_1d_exponential_params * )(p->params);

  t2 = rho[0] - rho[1];
  t3 = rho[0] + rho[1];
  t4 = 0.1e1 / t3;
  t5 = t2 * t4;
  t7 = 0.1e1 + t5 <= p->zeta_threshold;
  t8 = rho[0] <= p->dens_threshold || t7;
  t9 = p->zeta_threshold - 0.1e1;
  t11 = 0.1e1 - t5 <= p->zeta_threshold;
  t12 = -t9;
  t13 = my_piecewise5(t7, t9, t11, t12, t5);
  t14 = 0.1e1 + t13;
  t15 = t14 * M_PI;
  t16 = params->beta * t3;
  t17 = t15 * t16;
  t18 = xc_integrate(func1, NULL, 0.0, t17);
  t20 = xc_integrate(func2, NULL, 0.0, t17);
  t21 = 0.1e1 / M_PI;
  t22 = t20 * t21;
  t23 = 0.1e1 / params->beta;
  t24 = t23 * t4;
  t29 = my_piecewise3(t8, 0, -0.79577471545947667883e-1 * (t14 * t18 - t22 * t24) * t23);
  t31 = rho[1] <= p->dens_threshold || t11;
  t32 = my_piecewise5(t11, t9, t7, t12, -t5);
  t33 = 0.1e1 + t32;
  t34 = t33 * M_PI;
  t35 = t34 * t16;
  t36 = xc_integrate(func1, NULL, 0.0, t35);
  t38 = xc_integrate(func2, NULL, 0.0, t35);
  t39 = t38 * t21;
  t44 = my_piecewise3(t31, 0, -0.79577471545947667883e-1 * (-t39 * t24 + t33 * t36) * t23);
  if(zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    zk[0] = t29 + t44;

#ifndef XC_DONT_COMPILE_VXC

  if(order < 1) return;


  t45 = t3 * t3;
  t46 = 0.1e1 / t45;
  t47 = t2 * t46;
  t48 = t4 - t47;
  t49 = my_piecewise5(t7, 0, t11, 0, t48);
  t51 = t23 * t46;
  t52 = t22 * t51;
  t56 = my_piecewise3(t8, 0, -0.79577471545947667883e-1 * (t49 * t18 + t52) * t23);
  t58 = my_piecewise5(t11, 0, t7, 0, -t48);
  t60 = t39 * t51;
  t64 = my_piecewise3(t31, 0, -0.79577471545947667883e-1 * (t58 * t36 + t60) * t23);
  if(vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    vrho[0] = t29 + t44 + t3 * (t56 + t64);

  t67 = -t4 - t47;
  t68 = my_piecewise5(t7, 0, t11, 0, t67);
  t73 = my_piecewise3(t8, 0, -0.79577471545947667883e-1 * (t68 * t18 + t52) * t23);
  t75 = my_piecewise5(t11, 0, t7, 0, -t67);
  t80 = my_piecewise3(t31, 0, -0.79577471545947667883e-1 * (t75 * t36 + t60) * t23);
  if(vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    vrho[1] = t29 + t44 + t3 * (t73 + t80);

#ifndef XC_DONT_COMPILE_FXC

  if(order < 2) return;


  t86 = 0.1e1 / t45 / t3;
  t87 = t2 * t86;
  t89 = -0.2e1 * t46 + 0.2e1 * t87;
  t90 = my_piecewise5(t7, 0, t11, 0, t89);
  t92 = t49 * M_PI;
  t94 = t15 * params->beta;
  t95 = t92 * t16 + t94;
  t96 = t49 * t95;
  t97 = xc_bessel_K0( t17);
  t100 = t95 * t97;
  t101 = t14 * t4;
  t103 = 0.20000000000000000000e1 * t100 * t101;
  t104 = t23 * t86;
  t106 = 0.2e1 * t22 * t104;
  t110 = my_piecewise3(t8, 0, -0.79577471545947667883e-1 * (t90 * t18 + 0.20e1 * t96 * t97 + t103 - t106) * t23);
  t112 = my_piecewise5(t11, 0, t7, 0, -t89);
  t114 = t58 * M_PI;
  t116 = t34 * params->beta;
  t117 = t114 * t16 + t116;
  t118 = t58 * t117;
  t119 = xc_bessel_K0( t35);
  t122 = t117 * t119;
  t123 = t33 * t4;
  t125 = 0.20000000000000000000e1 * t122 * t123;
  t127 = 0.2e1 * t39 * t104;
  t131 = my_piecewise3(t31, 0, -0.79577471545947667883e-1 * (t112 * t36 + 0.20e1 * t118 * t119 + t125 - t127) * t23);
  if(v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    v2rho2[0] = 0.2e1 * t56 + 0.2e1 * t64 + t3 * (t110 + t131);

  t134 = 0.2e1 * t87;
  t135 = my_piecewise5(t7, 0, t11, 0, t134);
  t137 = t68 * t95;
  t143 = my_piecewise3(t8, 0, -0.79577471545947667883e-1 * (t135 * t18 + 0.20e1 * t137 * t97 + t103 - t106) * t23);
  t144 = my_piecewise5(t11, 0, t7, 0, -t134);
  t146 = t75 * t117;
  t152 = my_piecewise3(t31, 0, -0.79577471545947667883e-1 * (t144 * t36 + 0.20e1 * t146 * t119 + t125 - t127) * t23);
  if(v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    v2rho2[1] = t56 + t64 + t73 + t80 + t3 * (t143 + t152);

  t158 = 0.2e1 * t46 + 0.2e1 * t87;
  t159 = my_piecewise5(t7, 0, t11, 0, t158);
  t161 = t68 * M_PI;
  t163 = t161 * t16 + t94;
  t164 = t68 * t163;
  t167 = t163 * t97;
  t173 = my_piecewise3(t8, 0, -0.79577471545947667883e-1 * (t159 * t18 + 0.20e1 * t164 * t97 + 0.20000000000000000000e1 * t167 * t101 - t106) * t23);
  t175 = my_piecewise5(t11, 0, t7, 0, -t158);
  t177 = t75 * M_PI;
  t179 = t177 * t16 + t116;
  t180 = t75 * t179;
  t183 = t179 * t119;
  t189 = my_piecewise3(t31, 0, -0.79577471545947667883e-1 * (t175 * t36 + 0.20e1 * t180 * t119 + 0.20000000000000000000e1 * t183 * t123 - t127) * t23);
  if(v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    v2rho2[2] = 0.2e1 * t73 + 0.2e1 * t80 + t3 * (t173 + t189);

#ifndef XC_DONT_COMPILE_KXC

  if(order < 3) return;


  t194 = t45 * t45;
  t195 = 0.1e1 / t194;
  t196 = t2 * t195;
  t198 = 0.6e1 * t86 - 0.6e1 * t196;
  t199 = my_piecewise5(t7, 0, t11, 0, t198);
  t204 = t90 * M_PI;
  t206 = t92 * params->beta;
  t208 = t204 * t16 + 0.2e1 * t206;
  t209 = t49 * t208;
  t212 = t95 * t95;
  t214 = xc_bessel_K1( t17);
  t217 = t208 * t97;
  t219 = 0.20000000000000000000e1 * t217 * t101;
  t220 = t212 * t214;
  t222 = 0.20000000000000000000e1 * t220 * t101;
  t223 = t49 * t4;
  t225 = 0.20000000000000000000e1 * t100 * t223;
  t226 = t14 * t46;
  t227 = t100 * t226;
  t228 = 0.60000000000000000000e1 * t227;
  t229 = t23 * t195;
  t231 = 0.6e1 * t22 * t229;
  t235 = my_piecewise3(t8, 0, -0.79577471545947667883e-1 * (t199 * t18 + 0.40e1 * t90 * t95 * t97 + 0.20e1 * t209 * t97 - 0.20e1 * t49 * t212 * t214 + t219 - t222 + t225 - t228 + t231) * t23);
  t237 = my_piecewise5(t11, 0, t7, 0, -t198);
  t242 = t112 * M_PI;
  t244 = t114 * params->beta;
  t246 = t242 * t16 + 0.2e1 * t244;
  t247 = t58 * t246;
  t250 = t117 * t117;
  t252 = xc_bessel_K1( t35);
  t255 = t246 * t119;
  t257 = 0.20000000000000000000e1 * t255 * t123;
  t258 = t250 * t252;
  t260 = 0.20000000000000000000e1 * t258 * t123;
  t261 = t58 * t4;
  t263 = 0.20000000000000000000e1 * t122 * t261;
  t264 = t33 * t46;
  t265 = t122 * t264;
  t266 = 0.60000000000000000000e1 * t265;
  t268 = 0.6e1 * t39 * t229;
  t272 = my_piecewise3(t31, 0, -0.79577471545947667883e-1 * (t237 * t36 + 0.40e1 * t112 * t117 * t119 + 0.20e1 * t247 * t119 - 0.20e1 * t58 * t250 * t252 + t257 - t260 + t263 - t266 + t268) * t23);
  if(v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    v3rho3[0] = 0.3e1 * t110 + 0.3e1 * t131 + t3 * (t235 + t272);

  t275 = 0.2e1 * t143;
  t276 = 0.2e1 * t152;
  t277 = 0.2e1 * t86;
  t278 = 0.6e1 * t196;
  t279 = t277 - t278;
  t280 = my_piecewise5(t7, 0, t11, 0, t279);
  t285 = t68 * t208;
  t294 = my_piecewise3(t8, 0, -0.79577471545947667883e-1 * (t280 * t18 + 0.40e1 * t135 * t95 * t97 + 0.20e1 * t285 * t97 - 0.20e1 * t68 * t212 * t214 + t219 - t222 + t225 - t228 + t231) * t23);
  t296 = my_piecewise5(t11, 0, t7, 0, -t279);
  t301 = t75 * t246;
  t310 = my_piecewise3(t31, 0, -0.79577471545947667883e-1 * (t296 * t36 + 0.40e1 * t144 * t117 * t119 + 0.20e1 * t301 * t119 - 0.20e1 * t75 * t250 * t252 + t257 - t260 + t263 - t266 + t268) * t23);
  if(v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    v3rho3[1] = t110 + t131 + t275 + t276 + t3 * (t294 + t310);

  t313 = -t277 - t278;
  t314 = my_piecewise5(t7, 0, t11, 0, t313);
  t319 = t135 * t163;
  t322 = t135 * M_PI;
  t324 = t161 * params->beta;
  t325 = t322 * t16 + t206 + t324;
  t326 = t68 * t325;
  t329 = t214 * t95;
  t332 = t325 * t97;
  t335 = t163 * t214;
  t336 = t95 * t14;
  t337 = t336 * t4;
  t342 = t167 * t226;
  t345 = t314 * t18 + 0.20e1 * t159 * t95 * t97 + 0.20e1 * t319 * t97 + 0.20e1 * t326 * t97 - 0.20e1 * t164 * t329 + 0.20000000000000000000e1 * t332 * t101 - 0.20000000000000000000e1 * t335 * t337 + 0.20000000000000000000e1 * t167 * t223 - 0.20000000000000000000e1 * t342 - 0.40000000000000000000e1 * t227 + t231;
  t348 = my_piecewise3(t8, 0, -0.79577471545947667883e-1 * t345 * t23);
  t350 = my_piecewise5(t11, 0, t7, 0, -t313);
  t355 = t144 * t179;
  t358 = t144 * M_PI;
  t360 = t177 * params->beta;
  t361 = t358 * t16 + t244 + t360;
  t362 = t75 * t361;
  t365 = t252 * t117;
  t368 = t361 * t119;
  t371 = t179 * t252;
  t372 = t117 * t33;
  t373 = t372 * t4;
  t378 = t183 * t264;
  t381 = t350 * t36 + 0.20e1 * t175 * t117 * t119 + 0.20e1 * t355 * t119 + 0.20e1 * t362 * t119 - 0.20e1 * t180 * t365 + 0.20000000000000000000e1 * t368 * t123 - 0.20000000000000000000e1 * t371 * t373 + 0.20000000000000000000e1 * t183 * t261 - 0.20000000000000000000e1 * t378 - 0.40000000000000000000e1 * t265 + t268;
  t384 = my_piecewise3(t31, 0, -0.79577471545947667883e-1 * t381 * t23);
  if(v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    v3rho3[2] = t275 + t276 + t173 + t189 + t3 * (t348 + t384);

  t390 = -0.6e1 * t86 - 0.6e1 * t196;
  t391 = my_piecewise5(t7, 0, t11, 0, t390);
  t393 = t159 * t163;
  t396 = t159 * M_PI;
  t399 = t396 * t16 + 0.2e1 * t324;
  t400 = t68 * t399;
  t403 = t163 * t163;
  t404 = t68 * t403;
  t407 = t399 * t97;
  t410 = t403 * t214;
  t413 = t68 * t4;
  t420 = my_piecewise3(t8, 0, -0.79577471545947667883e-1 * (t391 * t18 + 0.40e1 * t393 * t97 + 0.20e1 * t400 * t97 - 0.20e1 * t404 * t214 + 0.20000000000000000000e1 * t407 * t101 - 0.20000000000000000000e1 * t410 * t101 + 0.20000000000000000000e1 * t167 * t413 - 0.60000000000000000000e1 * t342 + t231) * t23);
  t422 = my_piecewise5(t11, 0, t7, 0, -t390);
  t424 = t175 * t179;
  t427 = t175 * M_PI;
  t430 = t427 * t16 + 0.2e1 * t360;
  t431 = t75 * t430;
  t434 = t179 * t179;
  t435 = t75 * t434;
  t438 = t430 * t119;
  t441 = t434 * t252;
  t444 = t75 * t4;
  t451 = my_piecewise3(t31, 0, -0.79577471545947667883e-1 * (t422 * t36 + 0.40e1 * t424 * t119 + 0.20e1 * t431 * t119 - 0.20e1 * t435 * t252 + 0.20000000000000000000e1 * t438 * t123 - 0.20000000000000000000e1 * t441 * t123 + 0.20000000000000000000e1 * t183 * t444 - 0.60000000000000000000e1 * t378 + t268) * t23);
  if(v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    v3rho3[3] = 0.3e1 * t173 + 0.3e1 * t189 + t3 * (t420 + t451);

#ifndef XC_DONT_COMPILE_LXC

  if(order < 4) return;


  t457 = 0.1e1 / t194 / t3;
  t458 = t2 * t457;
  t460 = -0.24e2 * t195 + 0.24e2 * t458;
  t461 = my_piecewise5(t7, 0, t11, 0, t460);
  t463 = t217 * t226;
  t464 = 0.80000000000000000000e1 * t463;
  t465 = t208 * t214;
  t467 = 0.60000000000000000000e1 * t465 * t337;
  t468 = t220 * t226;
  t469 = 0.80000000000000000000e1 * t468;
  t470 = t212 * t95;
  t475 = -t97 - 0.1e1 / t14 * t21 * t24 * t214;
  t478 = 0.20000000000000000000e1 * t470 * t475 * t101;
  t479 = t49 * t46;
  t480 = t100 * t479;
  t481 = 0.80000000000000000000e1 * t480;
  t482 = t14 * t86;
  t483 = t100 * t482;
  t484 = 0.24000000000000000000e2 * t483;
  t485 = t23 * t457;
  t487 = 0.24e2 * t22 * t485;
  t499 = t204 * params->beta;
  t501 = t199 * M_PI * t16 + 0.3e1 * t499;
  t512 = 0.20000000000000000000e1 * t501 * t97 * t101;
  t514 = 0.40000000000000000000e1 * t217 * t223;
  t516 = 0.40000000000000000000e1 * t220 * t223;
  t517 = t90 * t4;
  t519 = 0.20000000000000000000e1 * t100 * t517;
  t520 = t461 * t18 - t464 - t467 + t469 - t478 - t481 + t484 - t487 + 0.60e1 * t199 * t95 * t97 + 0.60e1 * t90 * t208 * t97 - 0.60e1 * t90 * t212 * t214 + 0.20e1 * t49 * t501 * t97 - 0.60e1 * t209 * t329 - 0.20e1 * t49 * t470 * t475 + t512 + t514 - t516 + t519;
  t523 = my_piecewise3(t8, 0, -0.79577471545947667883e-1 * t520 * t23);
  t525 = my_piecewise5(t11, 0, t7, 0, -t460);
  t527 = t255 * t264;
  t528 = 0.80000000000000000000e1 * t527;
  t529 = t246 * t252;
  t531 = 0.60000000000000000000e1 * t529 * t373;
  t532 = t258 * t264;
  t533 = 0.80000000000000000000e1 * t532;
  t534 = t250 * t117;
  t539 = -t119 - 0.1e1 / t33 * t21 * t24 * t252;
  t542 = 0.20000000000000000000e1 * t534 * t539 * t123;
  t543 = t58 * t46;
  t544 = t122 * t543;
  t545 = 0.80000000000000000000e1 * t544;
  t546 = t33 * t86;
  t547 = t122 * t546;
  t548 = 0.24000000000000000000e2 * t547;
  t550 = 0.24e2 * t39 * t485;
  t562 = t242 * params->beta;
  t564 = t237 * M_PI * t16 + 0.3e1 * t562;
  t575 = 0.20000000000000000000e1 * t564 * t119 * t123;
  t577 = 0.40000000000000000000e1 * t255 * t261;
  t579 = 0.40000000000000000000e1 * t258 * t261;
  t580 = t112 * t4;
  t582 = 0.20000000000000000000e1 * t122 * t580;
  t583 = t525 * t36 - t528 - t531 + t533 - t542 - t545 + t548 - t550 + 0.60e1 * t237 * t117 * t119 + 0.60e1 * t112 * t246 * t119 - 0.60e1 * t112 * t250 * t252 + 0.20e1 * t58 * t564 * t119 - 0.60e1 * t247 * t365 - 0.20e1 * t58 * t534 * t539 + t575 + t577 - t579 + t582;
  t586 = my_piecewise3(t31, 0, -0.79577471545947667883e-1 * t583 * t23);
  if(v4rho4 != NULL && (p->info->flags & XC_FLAGS_HAVE_LXC))
    v4rho4[0] = 0.4e1 * t235 + 0.4e1 * t272 + t3 * (t523 + t586);

  t608 = 0.12e2 * t195;
  t609 = 0.24e2 * t458;
  t610 = -t608 + t609;
  t611 = my_piecewise5(t7, 0, t11, 0, t610);
  t613 = -t464 - t467 + t469 - t478 - t481 + t484 - t487 + 0.60e1 * t280 * t95 * t97 + 0.60e1 * t135 * t208 * t97 - 0.60e1 * t135 * t212 * t214 + 0.20e1 * t68 * t501 * t97 + t512 + t514 - t516 + t519 - 0.60e1 * t285 * t329 - 0.20e1 * t68 * t470 * t475 + t611 * t18;
  t616 = my_piecewise3(t8, 0, -0.79577471545947667883e-1 * t613 * t23);
  t635 = my_piecewise5(t11, 0, t7, 0, -t610);
  t637 = -t528 - t531 + t533 - t542 - t545 + t548 - t550 + 0.60e1 * t296 * t117 * t119 + 0.60e1 * t144 * t246 * t119 - 0.60e1 * t144 * t250 * t252 + 0.20e1 * t75 * t564 * t119 + t575 + t577 - t579 + t582 - 0.60e1 * t301 * t365 - 0.20e1 * t75 * t534 * t539 + t635 * t36;
  t640 = my_piecewise3(t31, 0, -0.79577471545947667883e-1 * t637 * t23);
  if(v4rho4 != NULL && (p->info->flags & XC_FLAGS_HAVE_LXC))
    v4rho4[1] = t235 + t272 + 0.3e1 * t294 + 0.3e1 * t310 + t3 * (t616 + t640);

  t665 = 0.2e1 * t322 * params->beta;
  t666 = t280 * M_PI * t16 + t499 + t665;
  t674 = t332 * t226;
  t676 = t167 * t479;
  t678 = t167 * t482;
  t682 = 0.40e1 * t314 * t95 * t97 + 0.20e1 * t159 * t208 * t97 - 0.20e1 * t159 * t212 * t214 + 0.20e1 * t280 * t163 * t97 + 0.40e1 * t135 * t325 * t97 + 0.20e1 * t68 * t666 * t97 + 0.20000000000000000000e2 * t483 - t487 - 0.20e1 * t164 * t475 * t212 - 0.40000000000000000000e1 * t674 - 0.40000000000000000000e1 * t676 + 0.40000000000000000000e1 * t678 - 0.40e1 * t319 * t329;
  t694 = t325 * t214;
  t705 = t335 * t336 * t46;
  t715 = my_piecewise5(t7, 0, t11, 0, t609);
  t717 = -0.40e1 * t326 * t329 - 0.20e1 * t164 * t465 + 0.20000000000000000000e1 * t666 * t97 * t101 + 0.40000000000000000000e1 * t332 * t223 + 0.20000000000000000000e1 * t167 * t517 - 0.40000000000000000000e1 * t694 * t337 - 0.20000000000000000000e1 * t335 * t208 * t14 * t4 - 0.40000000000000000000e1 * t335 * t96 * t4 + 0.40000000000000000000e1 * t705 - 0.20000000000000000000e1 * t163 * t475 * t212 * t14 * t4 - 0.40000000000000000000e1 * t463 + 0.40000000000000000000e1 * t468 - 0.40000000000000000000e1 * t480 + t715 * t18;
  t721 = my_piecewise3(t8, 0, -0.79577471545947667883e-1 * (t682 + t717) * t23);
  t740 = 0.2e1 * t358 * params->beta;
  t741 = t296 * M_PI * t16 + t562 + t740;
  t752 = t368 * t264;
  t754 = 0.40e1 * t350 * t117 * t119 + 0.20e1 * t175 * t246 * t119 - 0.20e1 * t175 * t250 * t252 + 0.20e1 * t296 * t179 * t119 + 0.40e1 * t144 * t361 * t119 + 0.20e1 * t75 * t741 * t119 - 0.40000000000000000000e1 * t527 + 0.40000000000000000000e1 * t532 - 0.40000000000000000000e1 * t544 + 0.20000000000000000000e2 * t547 - t550 - 0.20e1 * t180 * t539 * t250 - 0.40000000000000000000e1 * t752;
  t755 = t183 * t543;
  t757 = t183 * t546;
  t772 = t361 * t252;
  t783 = t371 * t372 * t46;
  t790 = my_piecewise5(t11, 0, t7, 0, -t609);
  t792 = -0.40000000000000000000e1 * t755 + 0.40000000000000000000e1 * t757 - 0.40e1 * t355 * t365 - 0.40e1 * t362 * t365 - 0.20e1 * t180 * t529 + 0.20000000000000000000e1 * t741 * t119 * t123 + 0.40000000000000000000e1 * t368 * t261 + 0.20000000000000000000e1 * t183 * t580 - 0.40000000000000000000e1 * t772 * t373 - 0.20000000000000000000e1 * t371 * t246 * t33 * t4 - 0.40000000000000000000e1 * t371 * t118 * t4 + 0.40000000000000000000e1 * t783 - 0.20000000000000000000e1 * t179 * t539 * t250 * t33 * t4 + t790 * t36;
  t796 = my_piecewise3(t31, 0, -0.79577471545947667883e-1 * (t754 + t792) * t23);
  if(v4rho4 != NULL && (p->info->flags & XC_FLAGS_HAVE_LXC))
    v4rho4[2] = 0.2e1 * t294 + 0.2e1 * t310 + 0.2e1 * t348 + 0.2e1 * t384 + t3 * (t721 + t796);

  t815 = t396 * params->beta;
  t816 = t314 * M_PI * t16 + t665 + t815;
  t827 = t407 * t226;
  t829 = t410 * t226;
  t832 = t167 * t68 * t46;
  t834 = 0.20e1 * t391 * t95 * t97 + 0.40e1 * t314 * t163 * t97 + 0.40e1 * t159 * t325 * t97 + 0.20e1 * t135 * t399 * t97 + 0.20e1 * t68 * t816 * t97 - 0.20e1 * t135 * t403 * t214 + 0.12000000000000000000e2 * t483 - t487 - 0.60000000000000000000e1 * t674 - 0.60000000000000000000e1 * t676 + 0.12000000000000000000e2 * t678 - 0.20000000000000000000e1 * t827 + 0.20000000000000000000e1 * t829 - 0.20000000000000000000e1 * t832;
  t857 = t399 * t214;
  t869 = t608 + t609;
  t870 = my_piecewise5(t7, 0, t11, 0, t869);
  t872 = -0.40e1 * t393 * t329 - 0.20e1 * t400 * t329 - 0.40e1 * t164 * t694 - 0.20e1 * t404 * t475 * t95 + 0.20000000000000000000e1 * t816 * t97 * t101 + 0.20000000000000000000e1 * t407 * t223 - 0.20000000000000000000e1 * t410 * t223 + 0.20000000000000000000e1 * t332 * t413 + 0.20000000000000000000e1 * t167 * t135 * t4 + 0.60000000000000000000e1 * t705 - 0.20000000000000000000e1 * t857 * t337 - 0.40000000000000000000e1 * t335 * t101 * t325 - 0.20000000000000000000e1 * t403 * t475 * t337 - 0.20000000000000000000e1 * t335 * t137 * t4 + t870 * t18;
  t876 = my_piecewise3(t8, 0, -0.79577471545947667883e-1 * (t834 + t872) * t23);
  t878 = my_piecewise5(t11, 0, t7, 0, -t869);
  t894 = t427 * params->beta;
  t895 = t350 * M_PI * t16 + t740 + t894;
  t906 = t438 * t264;
  t908 = t441 * t264;
  t910 = t878 * t36 + 0.20e1 * t422 * t117 * t119 + 0.40e1 * t350 * t179 * t119 + 0.40e1 * t175 * t361 * t119 + 0.20e1 * t144 * t430 * t119 + 0.20e1 * t75 * t895 * t119 - 0.20e1 * t144 * t434 * t252 + 0.12000000000000000000e2 * t547 - t550 - 0.60000000000000000000e1 * t752 - 0.60000000000000000000e1 * t755 + 0.12000000000000000000e2 * t757 - 0.20000000000000000000e1 * t906 + 0.20000000000000000000e1 * t908;
  t912 = t183 * t75 * t46;
  t936 = t430 * t252;
  t948 = -0.20000000000000000000e1 * t912 - 0.40e1 * t424 * t365 - 0.20e1 * t431 * t365 - 0.40e1 * t180 * t772 - 0.20e1 * t435 * t539 * t117 + 0.20000000000000000000e1 * t895 * t119 * t123 + 0.20000000000000000000e1 * t438 * t261 - 0.20000000000000000000e1 * t441 * t261 + 0.20000000000000000000e1 * t368 * t444 + 0.20000000000000000000e1 * t183 * t144 * t4 + 0.60000000000000000000e1 * t783 - 0.20000000000000000000e1 * t936 * t373 - 0.40000000000000000000e1 * t371 * t123 * t361 - 0.20000000000000000000e1 * t434 * t539 * t373 - 0.20000000000000000000e1 * t371 * t146 * t4;
  t952 = my_piecewise3(t31, 0, -0.79577471545947667883e-1 * (t910 + t948) * t23);
  if(v4rho4 != NULL && (p->info->flags & XC_FLAGS_HAVE_LXC))
    v4rho4[3] = 0.3e1 * t348 + 0.3e1 * t384 + t420 + t451 + t3 * (t876 + t952);

  t965 = t403 * t163;
  t970 = 0.24e2 * t195 + 0.24e2 * t458;
  t971 = my_piecewise5(t7, 0, t11, 0, t970);
  t981 = t391 * M_PI * t16 + 0.3e1 * t815;
  t1004 = -t487 + 0.24000000000000000000e2 * t678 - 0.80000000000000000000e1 * t827 + 0.80000000000000000000e1 * t829 - 0.80000000000000000000e1 * t832 - 0.60000000000000000000e1 * t857 * t163 * t14 * t4 - 0.20000000000000000000e1 * t965 * t475 * t101 + t971 * t18 - 0.60e1 * t400 * t335 - 0.20e1 * t68 * t965 * t475 + 0.20000000000000000000e1 * t981 * t97 * t101 + 0.40000000000000000000e1 * t407 * t413 - 0.40000000000000000000e1 * t410 * t413 + 0.20000000000000000000e1 * t167 * t159 * t4 + 0.60e1 * t391 * t163 * t97 + 0.60e1 * t159 * t399 * t97 - 0.60e1 * t159 * t403 * t214 + 0.20e1 * t68 * t981 * t97;
  t1007 = my_piecewise3(t8, 0, -0.79577471545947667883e-1 * t1004 * t23);
  t1016 = t434 * t179;
  t1021 = my_piecewise5(t11, 0, t7, 0, -t970);
  t1031 = t422 * M_PI * t16 + 0.3e1 * t894;
  t1054 = -t550 + 0.24000000000000000000e2 * t757 - 0.80000000000000000000e1 * t906 + 0.80000000000000000000e1 * t908 - 0.80000000000000000000e1 * t912 - 0.60000000000000000000e1 * t936 * t179 * t33 * t4 - 0.20000000000000000000e1 * t1016 * t539 * t123 + t1021 * t36 - 0.60e1 * t431 * t371 - 0.20e1 * t75 * t1016 * t539 + 0.20000000000000000000e1 * t1031 * t119 * t123 + 0.40000000000000000000e1 * t438 * t444 - 0.40000000000000000000e1 * t441 * t444 + 0.20000000000000000000e1 * t183 * t175 * t4 + 0.60e1 * t422 * t179 * t119 + 0.60e1 * t175 * t430 * t119 - 0.60e1 * t175 * t434 * t252 + 0.20e1 * t75 * t1031 * t119;
  t1057 = my_piecewise3(t31, 0, -0.79577471545947667883e-1 * t1054 * t23);
  if(v4rho4 != NULL && (p->info->flags & XC_FLAGS_HAVE_LXC))
    v4rho4[4] = 0.4e1 * t420 + 0.4e1 * t451 + t3 * (t1007 + t1057);

#ifndef XC_DONT_COMPILE_MXC

  if(order < 5) return;


#endif

#endif

#endif

#endif

#endif


}
