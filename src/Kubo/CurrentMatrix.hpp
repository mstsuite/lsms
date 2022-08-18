#include "Matrix.hpp"
#include "Complex.hpp"
#include "Array3d.hpp"
#include "Misc/ClebschGordan.hpp"
#include "Misc/Indices.hpp"
#include "Misc/integrator.hpp"
#include "Misc/integrateOneDim.hpp"
#include "SingleSite/AtomData.hpp"
#include "Main/SystemParameters.hpp"
#include "PhysicalConstants.hpp"
#include "SingleSite/SingleSiteScattering.hpp"

class CurrentMatrix {
public:
    Complex sqrtm1 = (0.0,1.0);
    Complex energy,prel,pnrel;
    int local_index, is,lmax_cg,kkrsz;
    AtomData *atom;
    NonRelativisticSingleScattererSolution solutionNonRel;
    Matrix<Complex> Jx, Jy, Jz;
    Array3d<Complex> zlrd;
    Array3d<Complex> tau1,tau2;
    CurrentMatrix(){};
    void init(LSMSSystemParameters &lsms, LocalTypeInfo &local, AtomData &a, int lindex, Complex en, int ispin);
    void calRadialSolutionDerivative(AtomData &a);
    Complex calPrefactor(int L, int Lp, int dir, int choice);
    Complex calRadialIntegral(AtomData &a, int L, int Lp, int dir, int choice);
    void assembleJxFromRadialIntegral(AtomData &a);
    void calJyzFromJx();
    void calculateNegativeEnergyCounterparts();
};
