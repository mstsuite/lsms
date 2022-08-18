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
     Complex efermi;
     Array3d <Complex> sigmatilde1, sigmatilde2, sigmatilde3, sigmatilde4;
     Array3d <Real> sigma;
     Matrix <CurrentMatrix> cm;
     Conductivity(LSMSSystemParameters &lsms, LSMSCommunication &comm, LocalTypeInfo &local);
};
