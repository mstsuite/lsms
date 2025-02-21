#ifndef LSMS_PHYSICALCONSTANTS_H
#define LSMS_PHYSICALCONSTANTS_H

// fine-structure constant alpha
constexpr double alphaInverse = 137.035999206;
constexpr double alpha = 1.0/alphaInverse;
constexpr double cphot = 2.0 * alphaInverse; // 274.072;
constexpr double c2inv = 1.0/(cphot*cphot);

// Joule in Rydberg
constexpr double convertJouleToRydberg = 1.380649e-23;
// Rydberg in eV
constexpr double convertRydbergToeV = 13.605698066;

// Boltzman constant in Ry/K
// defined by SI to be 1.380649e-23 J/K

// const double kBoltzmann = 6.3336e-6;
constexpr double kBoltzmann = 1.380649e-23 * 4.5874208973812E+17;

constexpr double convertKtoRydberg = 6.33361706838587e-06;

constexpr long double PI = 3.141592653589793238462643383279502884L;

#endif
