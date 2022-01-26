#include "Coeficients.hpp"
#include "Indices.hpp"

int SphericalHarmonicsCoeficients::lmax;
std::vector<Real> SphericalHarmonicsCoeficients::clm;

int GauntCoeficients::lmax;
Array3d<Real> GauntCoeficients::cgnt;

int IFactors::lmax;
Matrix<Complex> IFactors::illp;
std::vector<Complex> IFactors::ilp1;

int AngularMomentumIndices::lmax;
int AngularMomentumIndices::ndlj;
int AngularMomentumIndices::ndlm;
std::vector<int> AngularMomentumIndices::lofk;
std::vector<int> AngularMomentumIndices::mofk;
std::vector<int> AngularMomentumIndices::lofj;
std::vector<int> AngularMomentumIndices::mofj;

