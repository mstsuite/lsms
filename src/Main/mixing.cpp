#include "mixing.hpp"
#include "Communication/LSMSCommunication.hpp"
#include "mixing_params.hpp"


void
FrozenPotential::updateChargeDensity(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as) {
  for (int i = 0; i < as.size(); i++) {
    as[i].rhotot = as[i].rhoNew;
    as[i].xvalws[0] = as[i].xvalwsNew[0];
    as[i].xvalws[1] = as[i].xvalwsNew[1];
  }
}

void SimpleChargeDensityMixing::updateChargeDensity(LSMSCommunication &comm, LSMSSystemParameters &lsms,
                                                    std::vector<AtomData> &as) {
  for (int i = 0; i < as.size(); i++) {
    simpleMixing(&as[i].rhotot(0, 0),
                 &as[i].rhoNew(0, 0), as[i].rhotot.size(), alpha);
    simpleMixing(&as[i].xvalws[0],
                 &as[i].xvalwsNew[0], 2, alpha);
  }
}

void SimpleChargeDensityMixing::updatePotential(LSMSCommunication &comm, LSMSSystemParameters &lsms,
                                                std::vector<AtomData> &as) {
  for (int i = 0; i < as.size(); i++) {
    as[i].vr = as[i].vrNew;
    as[i].vdif = as[i].vdifNew;

  }

}

void SimplePotentialMixing::updateChargeDensity(LSMSCommunication &comm, LSMSSystemParameters &lsms,
                                                std::vector<AtomData> &as) {
  for (int i = 0; i < as.size(); i++) {
    as[i].rhotot = as[i].rhoNew;
    as[i].xvalws[0] = as[i].xvalwsNew[0];
    as[i].xvalws[1] = as[i].xvalwsNew[1];
  }
}

void
SimplePotentialMixing::updatePotential(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as) {
  for (int i = 0; i < as.size(); i++) {
    simpleMixing(&as[i].vr(0, 0), &as[i].vrNew(0, 0), as[i].vr.size(), alpha);
    simpleMixing(&as[i].vdif, &as[i].vdifNew, 1, alpha);

  }
}


void EfMixing::updateChargeDensity(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as) {
  lsms.chempot = alpha * lsms.chempot + (1.0 - alpha) * efOld;

}


void EfMixing::prepare(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as) {
  efOld = lsms.chempot;
}


void setupMixing(MixingParameters &mix, Mixing *&mixing, int iprint) {

  if (iprint >= 0)
    printf("\n");
  mixing = NULL;

  // frozen potential by default
  if (!mix.quantity[MixingParameters::no_mixing] &&
      !mix.quantity[MixingParameters::charge] &&
      !mix.quantity[MixingParameters::potential] &&
      !mix.quantity[MixingParameters::moment_magnitude] &&
      !mix.quantity[MixingParameters::moment_direction]) {
    mixing = new FrozenPotential;
    if (iprint >= 1)
      printf("Mixing method     : frozen potential (default)\n");
  }
    // no mixing
  else if (mix.quantity[MixingParameters::no_mixing]) {
    mixing = new NoMixing;
    if (iprint >= 1)
      printf("Mixing method     : no mixing\n");
  }
    // charge mixing
  else if (!mix.quantity[MixingParameters::no_mixing] &&
           mix.quantity[MixingParameters::charge] &&
           !mix.quantity[MixingParameters::potential] &&
           !mix.quantity[MixingParameters::moment_magnitude] &&
           !mix.quantity[MixingParameters::moment_direction]) {
    switch (mix.algorithm[MixingParameters::charge]) {
      case 1 :
        mixing = new SimpleChargeDensityMixing(mix.mixingParameter[MixingParameters::charge]);
        if (iprint >= 1)
          printf("Mixing method     : simple\n");
        break;
      case 2 :
        mixing = new BroydenChargeDensityMixing(mix.mixingParameter[MixingParameters::charge]);
        if (iprint >= 1)
          printf("Mixing method     : broyden\n");
        break;
      default :
        mixing = new NoMixing;
        if (iprint >= 1)
          printf("Mixing method     : no mixing\n");
    }
    if (iprint >= 1) {
      printf("Mixing quantity   : charge\n");
      printf("Mixing parameters : %4.2f\n", mix.mixingParameter[MixingParameters::charge]);
    }
  }
    // potential mixing
  else if (!mix.quantity[MixingParameters::no_mixing] &&
           !mix.quantity[MixingParameters::charge] &&
           mix.quantity[MixingParameters::potential] &&
           !mix.quantity[MixingParameters::moment_magnitude] &&
           !mix.quantity[MixingParameters::moment_direction]) {
    switch (mix.algorithm[MixingParameters::potential]) {
      case 1 :
        if (mix.mixingParameter[MixingParameters::potential] == 0.0) {
          mixing = new FrozenPotential;
          if (iprint >= 1)
            printf("Mixing method     : frozen potential\n");
        } else {
          mixing = new SimplePotentialMixing(mix.mixingParameter[MixingParameters::potential]);
          if (iprint >= 1)
            printf("Mixing method     : simple\n");
        }
        break;
      case 2 :
        mixing = new BroydenPotentialMixing(mix.mixingParameter[MixingParameters::potential]);
        if (iprint >= 1)
          printf("Mixing method     : broyden\n");
        break;
      default :
        mixing = new FrozenPotential;
        if (iprint >= 1)
          printf("Mixing method     : frozen potential\n");
    }
    if (iprint >= 1) {
      printf("Mixing quantity   : potential\n");
      printf("Mixing parameters : %4.2f\n", mix.mixingParameter[MixingParameters::potential]);
    }
  } else {
    if (iprint >= 1) {
      printf("Type of mixing is not supported.\n");
      for (int i = 0; i < mix.numQuantities; i++) {
        printf("quantity = %5d, algorithm = %5d, mixing parameter = %6.3f\n",
               mix.quantity[i], mix.algorithm[i], mix.mixingParameter[i]);
      }
      exit(1);
    }
  }

}

void BroydenPotentialMixing::updateChargeDensity(LSMSCommunication &comm, LSMSSystemParameters &lsms,
                                                 std::vector<AtomData> &as) {
  for (int i = 0; i < as.size(); i++) {
    as[i].rhotot = as[i].rhoNew;
    as[i].xvalws[0] = as[i].xvalwsNew[0];
    as[i].xvalws[1] = as[i].xvalwsNew[1];
  }
}

void BroydenPotentialMixing::prepare(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as) {
  vSize = 0;
  vStarts.resize(as.size());

  for (int i = 0; i < as.size(); i++) {
    vStarts[i] = vSize;
    vSize += as[i].vr.n_row();
  }
  mixer.init(alpha, 2 * vSize + 1);
  fNew.resize(2 * vSize + 1);
  fOld.resize(2 * vSize + 1);
}

void
BroydenPotentialMixing::updatePotential(LSMSCommunication &comm, LSMSSystemParameters &lsms,
                                        std::vector<AtomData> &as) {
  Real rms = 0.0;

  // first: copy potentials into fNew vector before mixing
  for (int i = 0; i < as.size(); i++) {
    for (int j = 0; j < as[i].vr.n_row(); j++) {
      fNew[vStarts[i] + j] = as[i].vrNew(j, 0);
      fNew[vStarts[i] + j + vSize] = as[i].vrNew(j, 1);
      fOld[vStarts[i] + j] = as[i].vr(j, 0);
      fOld[vStarts[i] + j + vSize] = as[i].vr(j, 1);
    }
    fNew[2 * vSize] = as[i].vdifNew;
    fOld[2 * vSize] = as[i].vdif;

    rms += as[i].vrms[0] + as[i].vrms[1];
  }
  rms = rms / (2.0 * as.size());
  globalSum(comm, rms);
  rms /= comm.size;

  // Broyden mixing
  mixer.mix(comm, fOld, fNew, rms);

  // copy mixed results back
  for (int i = 0; i < as.size(); i++) {
    for (int j = 0; j < as[i].vr.n_row(); j++) {
      as[i].vr(j, 0) = fNew[vStarts[i] + j];
      as[i].vr(j, 1) = fNew[vStarts[i] + j + vSize];
    }
    as[i].vdif = fNew[2 * vSize];
  }
}

void BroydenChargeDensityMixing::updateChargeDensity(LSMSCommunication &comm, LSMSSystemParameters &lsms,
                                                     std::vector<AtomData> &as) {
  Real rms = 0.0;

  // first: copy potentials into fNew vector before mixing
  for (int i = 0; i < as.size(); i++) {
    for (int j = 0; j < as[i].rhotot.n_row(); j++) {
      fNew[vStarts[i] + j] = as[i].rhoNew(j, 0);
      fNew[vStarts[i] + j + vSize] = as[i].rhoNew(j, 1);
      fOld[vStarts[i] + j] = as[i].rhotot(j, 0);
      fOld[vStarts[i] + j + vSize] = as[i].rhotot(j, 1);

    }
    fNew[2 * vSize] = as[i].xvalwsNew[0];
    fNew[2 * vSize + 1] = as[i].xvalwsNew[1];
    fOld[2 * vSize] = as[i].xvalws[0];
    fOld[2 * vSize + 1] = as[i].xvalws[1];

    rms += as[i].qrms[0] + as[i].qrms[1];
  }
  rms = rms / (2.0 * as.size());
  globalSum(comm, rms);
  rms /= comm.size;

  // Broyden mixing
  mixer.mix(comm, fOld, fNew, rms);

  // copy mixed results back
  for (int i = 0; i < as.size(); i++) {
    for (int j = 0; j < as[i].rhotot.n_row(); j++) {
      as[i].rhotot(j, 0) = fNew[vStarts[i] + j];
      as[i].rhotot(j, 1) = fNew[vStarts[i] + j + vSize];
    }
    as[i].xvalws[0] = fNew[2 * vSize];
    as[i].xvalws[1] = fNew[2 * vSize + 1];
  }
}

void BroydenChargeDensityMixing::updatePotential(LSMSCommunication &comm, LSMSSystemParameters &lsms,
                                                 std::vector<AtomData> &as) {
  for (int i = 0; i < as.size(); i++) {
    as[i].vr = as[i].vrNew;
    as[i].vdif = as[i].vdifNew;
  }
}

void
BroydenChargeDensityMixing::prepare(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as) {
  vSize = 0;
  vStarts.resize(as.size());

  for (int i = 0; i < as.size(); i++) {
    vStarts[i] = vSize;
    vSize += as[i].rhotot.n_row();
  }
  mixer.init(alpha, 2 * vSize + 2);
  fNew.resize(2 * vSize + 2);
  fOld.resize(2 * vSize + 2);
}
