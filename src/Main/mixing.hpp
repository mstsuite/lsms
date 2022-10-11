#ifndef LSMS_MIXING_H
#define LSMS_MIXING_H

#include "Real.hpp"
#include "SystemParameters.hpp"
#include "mixing_params.hpp"
#include "SingleSite/AtomData.hpp"
#include "Communication/LSMSCommunication.hpp"

#include <vector>
#include <deque>
#include <cmath>

#include "LAPACK.hpp"

template<typename T>
void simpleMixing(T *fold, T *fnew, int n, Real alpha) {
  if (alpha > 1.0) alpha = 1.0;
  if (alpha < 0.0) alpha = 0.0;
  Real beta = 1.0 - alpha;

  for (int i = 0; i < n; i++)
    fold[i] = alpha * fnew[i] + beta * fold[i];
}


/**
 * Modified Broyden method follows D. D. Johnson, PRB 38, 12807
**/
template<typename T>
class BroydenMixing {

  int vectorSize{};       // size of vectors to be mixed
  int currentIteration{}; // currently used number of iterations
  int iterationReset{};   // number of iterations after which the Broyden mixing is reset
  int maxBroydenLength{}; // maximum number of iterations that are used
  std::vector<T> vOld, F, dF;
  std::vector<Real> w, cm;
  std::vector<std::vector<T> > u, vt;
  std::vector<int> ipiv;
  Real alpha{}, w0{};
  Matrix<Real> a, b;

public:
  void save(std::vector<T> &fOld, std::vector<T> &fNew, Real wtmp) {
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


  // a <- a^-1
  void invert(Matrix<Real> &A, int nn) {
    int LDIM = A.l_dim();
    int LWORK = LDIM * LDIM;
    double *WORK = new double[LWORK];
    int INFO;

    LAPACK::dgetrf_(&nn, &nn, &A(0, 0), &LDIM, &ipiv[0], &INFO);
    LAPACK::dgetri_(&nn, &A(0, 0), &LDIM, &ipiv[0], WORK, &LWORK, &INFO);

    delete[] WORK;
  }

  void init(Real a_, int vs, int mbl = 10, int itres = 25) {
    alpha = a_;
    vectorSize = vs;
    currentIteration = 0;
    maxBroydenLength = mbl;
    iterationReset = itres;
    vOld.resize(vectorSize);
    F.resize(vectorSize);
    dF.resize(vectorSize);
    w.resize(maxBroydenLength);
    cm.resize(maxBroydenLength);
    a.resize(mbl, mbl);
    b.resize(mbl, mbl);
    ipiv.resize(mbl + 1);
    vt.resize(maxBroydenLength);
    u.resize(maxBroydenLength);
    for (int i = 0; i < maxBroydenLength; i++) {
      vt[i].resize(vectorSize);
      u[i].resize(vectorSize);
    }
    w0 = 0.01;
  }

  void mix(LSMSCommunication &comm, std::vector<T> &fOld, std::vector<T> &fNew, Real rms) {
    if (currentIteration == 0) {
      // first iteration: perform linear mixing, set up internal storage
      for (int i = 0; i < vectorSize; i++) {
        F[i] = fNew[i] - fOld[i];
        vOld[i] = fOld[i];
        fNew[i] = fOld[i] + alpha * F[i];
      }
    } else {
      int nn = std::min(currentIteration, maxBroydenLength);
      Real dFnorm = 0.0;
      Real Fnorm = 0.0;
      for (int i = 0; i < vectorSize; i++) {
        dF[i] = fNew[i] - fOld[i] - F[i];
        F[i] = fNew[i] - fOld[i];
        dFnorm += dF[i] * dF[i];
        Fnorm += F[i] * F[i];
      }
      ///* // Try without communication!
      globalSum(comm, dFnorm);
      globalSum(comm, Fnorm);
      //*/
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
      if (rms < 1.0e-9) wtmp = 2.0 * std::sqrt(0.01 / rms);
      if (wtmp < 1.0) wtmp = 1.0;

      save(fOld, fNew, wtmp);

      // off diagonal part of a(i,j)
      for (int j = 0; j < nn - 1; j++)
        for (int i = j + 1; i < nn; i++) {
          Real aij = 0.0;
          for (int k = 0; k < vectorSize; k++)
            aij += (vt[j][k]) * (vt[i][k]);
          a(i, j) = a(j, i) = aij;
        }
      // diagonal elements a(i,i) and cm(i)
      for (int i = 0; i < nn; i++) {
        Real aij = 0.0;
        Real cmj = 0.0;
        for (int k = 0; k < vectorSize; k++) {
          cmj += vt[i][k] * F[k];
          aij += vt[i][k] * vt[i][k];
        }
        a(i, i) = aij;
        cm[i] = cmj;
      }
      // sum over all sites
      globalSum(comm, &a(0, 0), maxBroydenLength * maxBroydenLength);
      globalSum(comm, &cm[0], maxBroydenLength);

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
      for (int k = 0; k < vectorSize; k++)
        fNew[k] = vOld[k] + alpha * F[k];

      for (int i = 0; i < nn; i++) {
        Real gmi = 0.0;
        for (int j = 0; j < nn; j++)
          gmi += cm[j] * b(j, i) * w[j];
        for (int k = 0; k < vectorSize; k++)
          fNew[k] = fNew[k] - gmi * u[i][k] * w[i];
      }

    }
    currentIteration++;
    if (iterationReset > 0 && currentIteration > iterationReset)
      currentIteration = 0;
  }
};

class MomentMixing {
public:
  virtual void update(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as) = 0;
};

class Mixing {
public:
  MomentMixing *momentMixing = nullptr;

  virtual ~Mixing() = default;

  // virtual void updateChargeDensity(LSMSSystemParameters &lsms, AtomData &a) = 0;
  virtual void
  updateChargeDensity(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as) = 0;

  // virtual void updatePotential(LSMSSystemParameters &lsms, AtomData &a) = 0;
  virtual void updatePotential(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as) = 0;

  virtual void prepare(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as) = 0;

  void updateMoments(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as);

};

class NoMixing : public Mixing {
public:

  void updateChargeDensity(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as) override {}

  void updatePotential(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as) override {}

  void prepare(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as) override {}
};


class FrozenPotential : public Mixing {

public:

  // void updateChargeDensity(LSMSSystemParameters &lsms, AtomData &a)
  // {
  //   a.rhotot = a.rhoNew;
  //   a.xvalws[0] = a.xvalwsNew[0];
  //   a.xvalws[1] = a.xvalwsNew[1];
  // }

  void updateChargeDensity(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as) override;

  void updatePotential(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as) override {}

  void prepare(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as) override {}

};


class SimpleChargeDensityMixing : public Mixing {

public:
  Real alpha;

  explicit SimpleChargeDensityMixing(Real _alpha) : alpha(_alpha) {}

  // void updateChargeDensity(LSMSSystemParameters &lsms, AtomData &a)
  // {
  //   simpleMixing( &a.rhotot(0,0), &a.rhoNew(0,0), a.rhotot.size(), alpha);
  //   simpleMixing( &a.xvalws[0], &a.xvalwsNew[0], 2, alpha);
  // }

  void updateChargeDensity(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as) override;

  void updatePotential(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as) override;

  void prepare(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as) override {}

};

class SimplePotentialMixing : public Mixing {

public:
  Real alpha;

  explicit SimplePotentialMixing(Real _alpha) : alpha(_alpha) {}

  void updateChargeDensity(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as) override;


  void updatePotential(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as) override;


  void prepare(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as) override {}

};


class EfMixing : public Mixing {
  Real efOld{};
  Real alpha{};

public:

  explicit EfMixing(Real _alpha) : alpha(_alpha) {}

  void updateChargeDensity(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as) override;

  void updatePotential(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as) override {};

  void prepare(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as) override;

};

class BroydenChargeDensityMixing : public Mixing {
  BroydenMixing<Real> mixer;
  std::vector<Real> fNew, fOld;
  int vSize{};
  std::vector<size_t> vStarts;
  Real alpha{};

public:
  explicit BroydenChargeDensityMixing(Real _alpha) : alpha(_alpha) {}

  void updateChargeDensity(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as) override;

  void updatePotential(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as) override;

  void prepare(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as) override;

};


class BroydenPotentialMixing : public Mixing {
  BroydenMixing<Real> mixer;
  std::vector<Real> fNew, fOld;
  size_t vSize{};
  std::vector<size_t> vStarts;
  Real alpha;

public:
  explicit BroydenPotentialMixing(Real _alpha) : alpha(_alpha) {}

  void updateChargeDensity(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as) override;

  void updatePotential(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as) override;

  void prepare(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as) override;

};

class SimpleMomentDirectionMixing : public MomentMixing {
  Real alpev;
public:
  explicit SimpleMomentDirectionMixing(Real _alpha) : alpev(_alpha) {}

  void update(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as);
};

void setupMixing(MixingParameters &mix, Mixing *&mixing, int iprint);


#endif
