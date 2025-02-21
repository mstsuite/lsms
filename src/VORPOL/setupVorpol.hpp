/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */

#ifndef LSMS_SETUP_VORPOL_H
#define LSMS_SETUP_VORPOL_H

#include "Communication/LSMSCommunication.hpp"
#include "Main/SystemParameters.hpp"
#include "Misc/Coeficients.hpp"
#include "VORPOL.hpp"

void setupVorpol(LSMSSystemParameters &lsms, CrystalParameters &crystal,
                 LocalTypeInfo &local);

void calculateVolumes(LSMSCommunication &comm, LSMSSystemParameters &lsms,
                      CrystalParameters &crystal, LocalTypeInfo &local);

#endif  // LSMS_SETUP_VORPOL_H