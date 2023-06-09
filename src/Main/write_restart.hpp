/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */
#ifndef WRITE_RESTART_HPP
#define WRITE_RESTART_HPP

#include "MixingParameter.hpp"
#include "Potential/PotentialShifter.hpp"
#include "SystemParameters.hpp"
#include "mixing_params.hpp"

int writeRestart(const char *restartName,
                 LSMSSystemParameters &lsms,
                 CrystalParameters &crystal,
                 lsms::MixingParameterPack &mix,
                 PotentialShifter &potentialShifter,
                 AlloyMixingDesc &alloyDesc);

int writeRestart(const char *restartName,
                 LSMSSystemParameters &lsms, CrystalParameters &crystal, MixingParameters &mix,
                 PotentialShifter &potentialShifter, AlloyMixingDesc &alloyDesc);

#endif
