#ifndef LSMS_PHYSICALCONSTANTS_H
#define LSMS_PHYSICALCONSTANTS_H

// fine-structure constant alpha
const double alphaInverse = 137.035999206;
const double alpha = 1.0/alphaInverse;
const double cphot = 2.0 * alphaInverse; // 274.072;
const double c2inv = 1.0/(cphot*cphot);

// Joule in Rydberg
const double convertJouleToRydberg = 1.380649e-23;
// Rydberg in eV
const double convertRydbergToeV = 13.605698066;

// Boltzman constant in Ry/K
// defined by SI to be 1.380649e-23 J/K

// const double kBoltzmann = 6.3336e-6;
const double kBoltzmann = 1.380649e-23 * 4.5874208973812E+17;

#endif
