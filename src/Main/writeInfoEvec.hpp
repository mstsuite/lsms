//
// Created by F.Moitzi on 13.12.2021.
//

#ifndef LSMS_WRITEINFOEVEC_HPP
#define LSMS_WRITEINFOEVEC_HPP

#include "Communication/LSMSCommunication.hpp"
#include "Real.hpp"

int writeInfoEvec(LSMSCommunication &comm, LSMSSystemParameters &lsms,
                  CrystalParameters &crystal, LocalTypeInfo &local, Real eband,
                  const char *name);

int writeLocalAtomData(LSMSCommunication &comm, LSMSSystemParameters &lsms,
                       CrystalParameters &crystal, LocalTypeInfo &local,
                       Real eband, const char *name);

int readInfoEvec(LSMSCommunication &comm, LSMSSystemParameters &lsms,
                 CrystalParameters &crystal, LocalTypeInfo &local,
                 const char *name);

#endif  // LSMS_WRITEINFOEVEC_HPP
