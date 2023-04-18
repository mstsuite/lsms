#ifndef LSMS_INITIALIZEATOMS_H
#define LSMS_INITIALIZEATOMS_H
#include <stdio.h>

#include "Communication/LSMSCommunication.hpp"
#include "Main/SystemParameters.hpp"
#include "SingleSite/AtomData.hpp"

void initializeAtom(AtomData &a);
void initializeNewPotentials(LSMSCommunication &comm, LSMSSystemParameters &lsms,
                            CrystalParameters &crystal, LocalTypeInfo &local);
int initializeNewAlloyBank(LSMSCommunication &comm, LSMSSystemParameters &lsms,
                           AlloyMixingDesc &alloyDesc,
                           AlloyAtomBank &alloyBank);

#endif
