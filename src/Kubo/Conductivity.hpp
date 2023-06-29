#include "Matrix.hpp"
#include "Complex.hpp"
#include "Array3d.hpp"
#include "Misc/Indices.hpp"
#include "Main/SystemParameters.hpp"
#include "Communication/LSMSCommunication.hpp"
#include "PhysicalConstants.hpp"
#include "SingleSite/SingleSiteScattering.hpp"
#include "Kubo/CurrentMatrix.hpp"

class Conductivity {
public:
     int pola;
     Real omega;
     Complex efermi;
     Array3d <Complex> sigmatilde1, sigmatilde2, sigmatilde3, sigmatilde4;
     Array3d <Real> sigma, rho;
     Matrix <Real> spin_summed_rho;
     Matrix <CurrentMatrix> cm;
     Conductivity(LSMSSystemParameters &lsms, LSMSCommunication &comm, LocalTypeInfo &local, Real volume);
     Complex calSigmaTilde(LocalTypeInfo &local, int dir1, int dir2, int is, int etype);
     void calSigma(LSMSCommunication &comm, LocalTypeInfo &local);
     void invertConductivityMatrix(int is);
     void writeSigmaTildeMat(Array3d <Complex> &stm, std::string matname);
     void writeRhoMat();
     void processTauMatrix(Matrix <Complex> &tau1, Matrix <Complex> &tau2, int m, int n, int is, int kkrsz, int etype);
     void processCurrentMatrix(LocalTypeInfo &local, int index, int dir, Matrix <Complex> &Jout1, int kkrsz, int etype);
};
