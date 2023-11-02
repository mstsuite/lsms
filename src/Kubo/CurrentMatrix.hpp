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
#include "MultipleScattering/linearSolvers.hpp"
#include "MultipleScattering/MultipleScattering.hpp"
#include "MultipleScattering/buildKKRMatrix.hpp"
#if defined(ACCELERATOR_LIBSCI) || defined(ACCELERATOR_CUDA_C) || defined(ACCELERATOR_HIP)
#include "Accelerator/DeviceStorage.hpp"
extern DeviceStorage *deviceStorage;
#endif

class CurrentMatrix {
public:
    const Complex sqrtm1 = Complex(0.0,1.0);
    Complex energy,prel,pnrel;
    int local_index, is,lmax_cg,kkrsz,nrmat;
    AtomData *atom;
    NonRelativisticSingleScattererSolution solutionNonRel;
    Matrix<Complex> Jx, Jy, Jz;
    Array3d<Complex> zlrd;
    Matrix<Complex> tau0,tau1,m,bigT;
    Complex *devM, *devT, *devTauFull;
    CurrentMatrix(){};
    void init(LSMSSystemParameters &lsms, LocalTypeInfo &local, AtomData &a, int lindex, Complex en, int ispin);
    void calRadialSolutionDerivative(AtomData &a);
    Complex calPrefactor(int L, int Lp, int dir, int choice);
    Complex calRadialIntegral(AtomData &a, int L, int Lp, int dir, int choice);
    void assembleJxFromRadialIntegral(AtomData &a);
    void calJyzFromJx();
    void calTauFull(LSMSSystemParameters &lsms, LocalTypeInfo &local, AtomData &a);
};
