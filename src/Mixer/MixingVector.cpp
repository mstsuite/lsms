//
// Created by F.Moitzi on 23.12.2022.
//

#include "MixingVector.hpp"

namespace lsms {

/**
 * Mix potential and difference
 */
PotentialMixingVector::PotentialMixingVector(LSMSSystemParameters &lsms,
                                             LocalTypeInfo &local) {
  rms = 0.0;
  size = 0;
  starts.resize(local.atom.size());

  for (auto i = 0; i < local.atom.size(); i++) {
    starts[i] = size;
    size += local.atom[i].jmt;
  }

  total_size = lsms.n_spin_pola * size + 1;

  vec_old.resize(total_size);
  vec_new.resize(total_size);
}

/**
 * Mix charge density
 */
ChargeMixingVector::ChargeMixingVector(LSMSSystemParameters &lsms,
                                       LocalTypeInfo &local,
                                       bool no_spin_split) : no_spin_split{no_spin_split} {
  rms = 0.0;
  size = 0;
  starts.resize(local.atom.size());

  for (auto i = 0; i < local.atom.size(); i++) {
    starts[i] = size;
    size += local.atom[i].jmt;
  }

  if (no_spin_split) {
    total_size = lsms.n_spin_pola * size;
  } else {
    total_size = size;
  }

  vec_old.resize(total_size);
  vec_new.resize(total_size);
}

/**
 * Mix spin density
 */
SpinMixingVector::SpinMixingVector(LSMSSystemParameters &lsms,
                                   LocalTypeInfo &local) {
  rms = 0.0;
  size = 0;
  starts.resize(local.atom.size());

  for (auto i = 0; i < local.atom.size(); i++) {
    starts[i] = size;
    size += local.atom[i].jmt;
  }

  total_size = size;

  vec_old.resize(total_size);
  vec_new.resize(total_size);
}

void PotentialMixingVector::copyToVector(LSMSCommunication &comm,
                                         LSMSSystemParameters &lsms,
                                         LocalTypeInfo &local) {
  Real local_rms = 0.0;

  for (auto i = 0; i < local.atom.size(); i++) {
    for (auto j = 0; j < local.atom[i].jmt; j++) {
      vec_new[starts[i] + j] = local.atom[i].vrNew(j, 0);

      if (lsms.n_spin_pola == 2) {
        vec_new[starts[i] + j + size] = local.atom[i].vrNew(j, 1);
      }

      vec_old[starts[i] + j] = local.atom[i].vr(j, 0);

      if (lsms.n_spin_pola == 2) {
        vec_old[starts[i] + j + size] = local.atom[i].vr(j, 1);
      }
    }

    vec_new[lsms.n_spin_pola * size] = local.atom[i].vdifNew;
    vec_old[lsms.n_spin_pola * size] = local.atom[i].vdif;

    for (auto is = 0; is < lsms.n_spin_pola; is++) {
      local_rms += local.atom[i].vrms[is];
    }
  }

  globalSum(comm, local_rms);
  rms = local_rms / (1.0 * lsms.num_atoms * lsms.n_spin_pola);
}

void PotentialMixingVector::copyFromVector(LSMSCommunication &comm,
                                           LSMSSystemParameters &lsms,
                                           LocalTypeInfo &local) {
  for (auto i = 0; i < local.atom.size(); i++) {
    for (auto j = 0; j < local.atom[i].jmt; j++) {
      local.atom[i].vrNew(j, 0) = vec_new[starts[i] + j];

      if (lsms.n_spin_pola == 2) {
        local.atom[i].vrNew(j, 1) = vec_new[starts[i] + j + size];
      }
    }

    local.atom[i].vdifNew = vec_new[lsms.n_spin_pola * size];
  }
}

void ChargeMixingVector::copyToVector(LSMSCommunication &comm,
                                      LSMSSystemParameters &lsms,
                                      LocalTypeInfo &local) {
  Real local_rms = 0.0;

  for (auto i = 0; i < local.atom.size(); i++) {
    for (auto j = 0; j < local.atom[i].jmt; j++) {
      vec_new[starts[i] + j] = local.atom[i].rhoNew(j, 0);

      if (lsms.n_spin_pola == 2 && no_spin_split) {
        vec_new[starts[i] + j + size] = local.atom[i].rhoNew(j, 1);
      }

      vec_old[starts[i] + j] = local.atom[i].rhotot(j, 0);

      if (lsms.n_spin_pola == 2 && no_spin_split) {
        vec_old[starts[i] + j + size] = local.atom[i].rhotot(j, 1);
      }
    }

    for (auto is = 0; is < lsms.n_spin_pola; is++) {
      local_rms += local.atom[i].qrms[is];
    }
  }

  globalSum(comm, local_rms);
  rms = local_rms / (1.0 * lsms.num_atoms * lsms.n_spin_pola);
}

void ChargeMixingVector::copyFromVector(LSMSCommunication &comm,
                                        LSMSSystemParameters &lsms,
                                        LocalTypeInfo &local) {
  for (auto i = 0; i < local.atom.size(); i++) {
    for (auto j = 0; j < local.atom[i].jmt; j++) {

      local.atom[i].rhoNew(j, 0) = vec_new[starts[i] + j];

      if (lsms.n_spin_pola == 2 && no_spin_split) {
        local.atom[i].rhoNew(j, 1) = vec_new[starts[i] + j + size];
      }

    }
  }
}

void SpinMixingVector::copyToVector(LSMSCommunication &comm,
                                    LSMSSystemParameters &lsms,
                                    LocalTypeInfo &local) {
  Real local_rms = 0.0;

  for (auto i = 0; i < local.atom.size(); i++) {
    for (auto j = 0; j < local.atom[i].jmt; j++) {

      if (lsms.n_spin_pola == 2) {
        vec_new[starts[i] + j] = local.atom[i].rhoNew(j, 1);
        vec_old[starts[i] + j] = local.atom[i].rhotot(j, 1);
      }
    }

    for (auto is = 0; is < lsms.n_spin_pola; is++) {
      local_rms += local.atom[i].qrms[is];
    }
  }

  globalSum(comm, local_rms);
  rms = local_rms / (1.0 * lsms.num_atoms * lsms.n_spin_pola);
}

void SpinMixingVector::copyFromVector(LSMSCommunication &comm,
                                      LSMSSystemParameters &lsms,
                                      LocalTypeInfo &local) {
  for (auto i = 0; i < local.atom.size(); i++) {
    for (auto j = 0; j < local.atom[i].jmt; j++) {

      if (lsms.n_spin_pola == 2) {
        local.atom[i].rhoNew(j, 1) = vec_new[starts[i] + j];
      }

    }
  }
}

DefaultMixerVector::DefaultMixerVector(std::size_t t) {
  rms = 0.0;
  size = t;
  total_size = t;
  vec_new.resize(t);
  vec_old.resize(t);
}
void DefaultMixerVector::copyToVector(LSMSCommunication &comm,
                                      LSMSSystemParameters &lsms,
                                      LocalTypeInfo &local) {}

void DefaultMixerVector::copyFromVector(LSMSCommunication &comm,
                                        LSMSSystemParameters &lsms,
                                        LocalTypeInfo &local) {}
}  // namespace lsms