//
// Created by F.Moitzi on 11.10.2022.
//

#ifndef LSMS_SETUPCELL_HPP
#define LSMS_SETUPCELL_HPP

#include "Communication/LSMSCommunication.hpp"
#include "Main/SystemParameters.hpp"
#include "SingleSite/AtomData.hpp"

namespace lsms {

void setupCell(LSMSSystemParameters &lsms, CrystalParameters &crystal,
               LocalTypeInfo &local);

void printCell(LSMSSystemParameters &lsms, CrystalParameters &crystal,
               LocalTypeInfo &local, LSMSCommunication &comm);

}  // namespace lsms

#endif  // LSMS_SETUPCELL_HPP
