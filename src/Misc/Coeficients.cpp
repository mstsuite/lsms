#include "Coeficients.hpp"

int SphericalHarmonicsCoeficients::lmax;
std::vector<Real> SphericalHarmonicsCoeficients::clm;

int GauntCoeficients::lmax;
Array3d<Real> GauntCoeficients::cgnt;

int IFactors::lmax;
Matrix<Complex> IFactors::illp;
std::vector<Complex> IFactors::ilp1;
