//
// Created by F.Moitzi on 12.12.2022.
//

#include "Mixer.hpp"

#include "fmt/core.h"
#include "fmt/format.h"

namespace lsms {

std::unique_ptr<AbstractMixer> AbstractMixer::generateMixer(
    unsigned int mixerType, const RealMixingVector &mixVector,
    const MixingParameter &params) {
  if (mixerType == MixerType::SIMPLE_MIXER) {
    return std::make_unique<SimpleMixer>(params.alpha);
  } else if (mixerType == MixerType::BROYDEN_MIXER) {
    return std::make_unique<BroydenMixer>(mixVector.total_size, params.alpha,
                                          params.max_broyden, params.iter_reset,
                                          params.w0);
  } else if (mixerType == MixerType::NO_MIXER) {
    return std::make_unique<NoMixer>();
  } else {
    throw std::runtime_error("Mixing type doesn't exists");
  }
}

void NoMixer::mix(LSMSCommunication &comm, RealMixingVector &mix_vector) {}

void NoMixer::mix(RealMixingVector &mix_vector) {}

SimpleMixer::SimpleMixer(Real alpha) : alpha{alpha} {}

void SimpleMixer::mix(LSMSCommunication &comm, RealMixingVector &mix) {
  this->mix(mix);
}

void SimpleMixer::mix(RealMixingVector &mix) {
  Real beta = 1.0 - alpha;

  for (auto i = 0; i < mix.total_size; i++) {
    mix.vec_new[i] = beta * mix.vec_old[i] + alpha * mix.vec_new[i];
  }
}

BroydenMixer::BroydenMixer(std::size_t size, Real alpha, std::size_t maxBroyden,
                           std::size_t iterReset, Real w0)
    : vectorSize{size},
      iterationReset{iterReset},
      maxBroydenLength{maxBroyden},
      alpha{alpha},
      w0{w0} {
  currentIteration = 0;

  vOld.resize(vectorSize);
  F.resize(vectorSize);
  dF.resize(vectorSize);
  w.resize(maxBroydenLength);
  cm.resize(maxBroydenLength);
  a.resize(maxBroyden, maxBroyden);
  b.resize(maxBroyden, maxBroyden);
  ipiv.resize(maxBroyden + 1);
  vt.resize(maxBroydenLength);
  u.resize(maxBroydenLength);

  for (auto i = 0; i < maxBroydenLength; i++) {
    vt[i].resize(vectorSize);
    u[i].resize(vectorSize);
  }

  auto LWORK = maxBroydenLength * maxBroydenLength;
  WORK.resize(LWORK);
}

void BroydenMixer::save(std::vector<Real> &fOld, std::vector<Real> &fNew,
                        Real wtmp) {
  if (currentIteration < maxBroydenLength + 1) {
    for (int i = 0; i < vectorSize; i++) {
      u[currentIteration - 1][i] = fNew[i];
      vt[currentIteration - 1][i] = fOld[i];
    }

    w[currentIteration - 1] = wtmp;

  } else {
    for (int j = 0; j < maxBroydenLength - 1; j++) {
      for (int i = 0; i < vectorSize; i++) {
        u[j][i] = u[j + 1][i];
        vt[j][i] = vt[j + 1][i];
      }

      w[j] = w[j + 1];
    }

    for (int i = 0; i < vectorSize; i++) {
      u[maxBroydenLength - 1][i] = fNew[i];
      vt[maxBroydenLength - 1][i] = fOld[i];
    }
    w[maxBroydenLength - 1] = wtmp;
  }
}

void BroydenMixer::invert(Matrix<Real> &A, int nn) {
  int LDIM = A.l_dim();
  int LWORK = LDIM * LDIM;
  int INFO;

  // computes an LU factorization
  LAPACK::dgetrf_(&nn, &nn, &A(0, 0), &LDIM, ipiv.data(), &INFO);

  // computes the inverse of a matrix using the computed LU factorization
  LAPACK::dgetri_(&nn, &A(0, 0), &LDIM, ipiv.data(), WORK.data(), &LWORK,
                  &INFO);
}

void BroydenMixer::mix(LSMSCommunication &comm, RealMixingVector &mix) {
  auto &fOld = mix.vec_old;
  auto &fNew = mix.vec_new;
  auto &rms = mix.rms;

  if (currentIteration == 0) {
    // first iteration: perform linear mixing, set up internal storage
    for (auto i = 0; i < vectorSize; i++) {
      F[i] = fNew[i] - fOld[i];
      vOld[i] = fOld[i];
      fNew[i] = fOld[i] + alpha * F[i];
    }

  } else {
    auto nn = std::min(currentIteration, maxBroydenLength);
    Real dFnorm = 0.0;
    Real Fnorm = 0.0;
    for (auto i = 0; i < vectorSize; i++) {
      dF[i] = fNew[i] - fOld[i] - F[i];
      F[i] = fNew[i] - fOld[i];
      dFnorm += dF[i] * dF[i];
      Fnorm += F[i] * F[i];
    }



    // Realry without communication!
    globalSum(comm, dFnorm);
    globalSum(comm, Fnorm);

    dFnorm = std::sqrt(dFnorm);
    Fnorm = std::sqrt(Fnorm);

    Real fac2 = 1.0 / dFnorm;
    Real fac1 = alpha * fac2;

    for (int i = 0; i < vectorSize; i++) {
      fNew[i] = fac1 * dF[i] + fac2 * (fOld[i] - vOld[i]);
      vOld[i] = fOld[i];
      fOld[i] = fac2 * dF[i];
    }

//    fmt::print("{:13.8f}\n", fmt::join(fNew, " "));
//    fmt::print("{:13.8f}\n", fmt::join(vOld, " "));
//    fmt::print("{:13.8f}\n", fmt::join(fOld, " "));

    Real wtmp = 0.0;

    if (rms > 1.0e-9) wtmp = 2.0 * std::sqrt(0.01 / rms);
    if (wtmp < 1.0) wtmp = 1.0;

    save(fOld, fNew, wtmp);

    // off diagonal part of a(i,j)
    for (auto j = 0; j < nn - 1; j++) {
      for (auto i = j + 1; i < nn; i++) {
        Real aij = 0.0;
        for (auto k = 0; k < vectorSize; k++) aij += (vt[j][k]) * (vt[i][k]);
        a(i, j) = a(j, i) = aij;
      }
    }

    // diagonal elements a(i,i) and cm(i)
    for (auto i = 0; i < nn; i++) {
      Real aij = 0.0;
      Real cmj = 0.0;
      for (auto k = 0; k < vectorSize; k++) {
        cmj += vt[i][k] * F[k];
        aij += vt[i][k] * vt[i][k];
      }
      a(i, i) = aij;
      cm[i] = cmj;
    }
    // sum over all sites
    globalSum(comm, &a(0, 0), maxBroydenLength * maxBroydenLength);
    globalSum(comm, cm.data(), maxBroydenLength);

//    std::cout << "cm" << std::endl;
//    std::cout << cm[0] << std::endl;
//    fmt::print("{:13.8f}\n", fmt::join(vt[0], " "));
//    fmt::print("{:13.8f}\n", fmt::join(F, " "));



    // now calculate the b-matrix
    //  b = [w(0)*w(0)*delta(i,j) + w(i)*w(j)*a(i,j)]^-1
    //
    for (int i = 0; i < nn; i++) {
      for (int j = 0; j < nn; j++) {
        b(i, j) = a(i, j) * w[j] * w[i];
      }
      b(i, i) = w0 * w0 + a(i, i) * w[i] * w[i];
    }

    invert(b, nn);

//    std::cout << "vOld" << std::endl;
//    fmt::print("{:13.8f}\n", fmt::join(vOld, " "));
//    std::cout << "F" << std::endl;
//    fmt::print("{:13.8f}\n", fmt::join(F, " "));

    // mix vectors
    for (int k = 0; k < vectorSize; k++) fNew[k] = vOld[k] + alpha * F[k];

//    std::cout << "FNew" << std::endl;
//    fmt::print("{:13.8f}\n", fmt::join(fNew, " "));

    for (int i = 0; i < nn; i++) {
      Real gmi = 0.0;
      for (int j = 0; j < nn; j++) gmi += cm[j] * b(j, i) * w[j];
      for (int k = 0; k < vectorSize; k++)
        fNew[k] = fNew[k] - gmi * u[i][k] * w[i];
    }

//    std::cout << "FNew" << std::endl;
//    fmt::print("{:13.8f}\n", fmt::join(fNew, " "));
  }
  currentIteration++;
  if (iterationReset > 0 && currentIteration > iterationReset)
    currentIteration = 0;
}

void BroydenMixer::mix(RealMixingVector &mix) {
  auto &fOld = mix.vec_old;
  auto &fNew = mix.vec_new;
  auto &rms = mix.rms;

  if (currentIteration == 0) {
    // first iteration: perform linear mixing, set up internal storage
    for (auto i = 0; i < vectorSize; i++) {
      F[i] = fNew[i] - fOld[i];
      vOld[i] = fOld[i];
      fNew[i] = fOld[i] + alpha * F[i];
    }

  } else {
    auto nn = std::min(currentIteration, maxBroydenLength);
    Real dFnorm = 0.0;
    Real Fnorm = 0.0;
    for (auto i = 0; i < vectorSize; i++) {
      dF[i] = fNew[i] - fOld[i] - F[i];
      F[i] = fNew[i] - fOld[i];
      dFnorm += dF[i] * dF[i];
      Fnorm += F[i] * F[i];
    }

    dFnorm = std::sqrt(dFnorm);
    Fnorm = std::sqrt(Fnorm);

    Real fac2 = 1.0 / dFnorm;
    Real fac1 = alpha * fac2;

    for (int i = 0; i < vectorSize; i++) {
      fNew[i] = fac1 * dF[i] + fac2 * (fOld[i] - vOld[i]);
      vOld[i] = fOld[i];
      fOld[i] = fac2 * dF[i];
    }

    Real wtmp = 0.0;

    if (rms > 1.0e-9) wtmp = 2.0 * std::sqrt(0.01 / rms);
    if (wtmp < 1.0) wtmp = 1.0;

    save(fOld, fNew, wtmp);

    // off diagonal part of a(i,j)
    for (auto j = 0; j < nn - 1; j++)
      for (auto i = j + 1; i < nn; i++) {
        Real aij = 0.0;
        for (auto k = 0; k < vectorSize; k++) aij += (vt[j][k]) * (vt[i][k]);
        a(i, j) = a(j, i) = aij;
      }
    // diagonal elements a(i,i) and cm(i)
    for (auto i = 0; i < nn; i++) {
      Real aij = 0.0;
      Real cmj = 0.0;
      for (auto k = 0; k < vectorSize; k++) {
        cmj += vt[i][k] * F[k];
        aij += vt[i][k] * vt[i][k];
      }
      a(i, i) = aij;
      cm[i] = cmj;
    }

    // now calculate the b-matrix
    //  b = [w(0)*w(0)*delta(i,j) + w(i)*w(j)*a(i,j)]^-1
    //
    for (int i = 0; i < nn; i++) {
      for (int j = 0; j < nn; j++) {
        b(i, j) = a(i, j) * w[j] * w[i];
      }
      b(i, i) = w0 * w0 + a(i, i) * w[i] * w[i];
    }

    invert(b, nn);

    // mix vectors
    for (int k = 0; k < vectorSize; k++) fNew[k] = vOld[k] + alpha * F[k];

    for (int i = 0; i < nn; i++) {
      Real gmi = 0.0;
      for (int j = 0; j < nn; j++) gmi += cm[j] * b(j, i) * w[j];
      for (int k = 0; k < vectorSize; k++)
        fNew[k] = fNew[k] - gmi * u[i][k] * w[i];
    }
  }
  currentIteration++;
  if (iterationReset > 0 && currentIteration > iterationReset)
    currentIteration = 0;
}

void printMixingParameters(const MixingParameterPack &mix, LSMSCommunication &comm, LSMSSystemParameters &lsms) {

  if (comm.rank == 0) {
    fmt::print("============ {} ============\n", "Mixing parameters");

    fmt::print("Charge mixing method:    {}\n",
               lsms::getMixerName(mix.chd_mixer_type));

    if (lsms.n_spin_pola == 2) {
      fmt::print("Spin mixing method:      {}\n",
                 lsms::getMixerName(mix.spd_mixer_type));
    }

    fmt::print("Potential mixing method: {}\n",
               lsms::getMixerName(mix.pot_mixer_type));

    fmt::print("Initial mixing method:   {}\n",
               lsms::getMixerName(mix.init_mixer_type));

    fmt::print("===========================================\n");
  }

}

}  // namespace lsms