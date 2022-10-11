#ifndef LSMS_POTENTIALIO_H
#define LSMS_POTENTIALIO_H

#include "Communication/LSMSCommunication.hpp"
#include "Main/SystemParameters.hpp"

int loadPotentials(LSMSCommunication &comm, LSMSSystemParameters &lsms,
                   CrystalParameters &crystal, LocalTypeInfo &local);
int writePotentials(LSMSCommunication &comm, LSMSSystemParameters &lsms,
                    CrystalParameters &crystal, LocalTypeInfo &local);
void initialAtomSetup(LSMSCommunication &comm, LSMSSystemParameters &lsms,
                      CrystalParameters &crystal, LocalTypeInfo &local);

int loadAlloyBank(LSMSCommunication &comm, LSMSSystemParameters &lsms,
                  AlloyMixingDesc &alloyDesc, AlloyAtomBank &alloyBank);

#endif
