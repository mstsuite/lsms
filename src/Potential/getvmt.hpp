#include "Communication/LSMSCommunication.hpp"
#include "Main/SystemParameters.hpp"
#include "Real.hpp"
#include "SingleSite/AtomData.hpp"

// getvmt.f in LSMS 1
void getvmt(LSMSSystemParameters &lsms, AtomData &atom,
            CrystalParameters &crystal, std::vector<Real> &qsub, int &mytype,
            Real &vmt, Real &vmt1, Real &u0, Real &u0MT);
