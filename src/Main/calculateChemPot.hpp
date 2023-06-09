#ifndef LSMS_CALCULATECEMPOT_HPP
#define LSMS_CALCULATECEMPOT_HPP

#include "SystemParameters.hpp"
#include "Communication/LSMSCommunication.hpp"

namespace lsms {

void calculateChemPot(LSMSCommunication &comm, LSMSSystemParameters &lsms, LocalTypeInfo &local,
                      Real &eigensum);

}

#endif
