//
// Created by F.Moitzi on 12.12.2022.
//

#ifndef LSMS_MIXER_HPP
#define LSMS_MIXER_HPP

#include "Communication/LSMSCommunication.hpp"
#include "Main/SystemParameters.hpp"
#include "MixerType.hpp"
#include "MixingParameter.hpp"
#include "MixingVector.hpp"

namespace lsms {

/**
 * Abstract Mixer
 */
class AbstractMixer {
 public:
  virtual void mix(LSMSCommunication &comm, RealMixingVector &mix_vector) = 0;

  virtual void mix(RealMixingVector &mix_vector) = 0;

  static std::unique_ptr<AbstractMixer> generateMixer(
      unsigned int mixerType, const RealMixingVector &mixVector,
      const MixingParameter &params);

};

/**
 * No Mixer
 */

class NoMixer : public AbstractMixer {
 public:
  void mix(LSMSCommunication &comm, RealMixingVector &mix_vector) override;

  void mix(RealMixingVector &mix_vector) override;

};

/**
 * Simple Mixer
 */

class SimpleMixer : public AbstractMixer {
  Real alpha;

 public:
  explicit SimpleMixer(Real alpha);

  void mix(LSMSCommunication &comm, RealMixingVector &mix_vector) override;

  void mix(RealMixingVector &mix_vector) override;

};

/**
 * Broyden Mixer
 */
class BroydenMixer : public AbstractMixer {
 private:
  std::size_t vectorSize;        // size of vectors to be mixed
  std::size_t currentIteration;  // currently used number of iterations
  std::size_t
      iterationReset;  // number of iterations after which the Broyden mixing
  // is reset
  std::size_t maxBroydenLength;  // maximum number of iterations that are used
  std::vector<Real> vOld, F, dF;
  std::vector<Real> w, cm;
  std::vector<std::vector<Real>> u, vt;
  std::vector<int> ipiv;
  Real alpha{}, w0{};
  Matrix<Real> a, b;

  std::vector<Real> WORK;

  void save(std::vector<Real> &fOld, std::vector<Real> &fNew, Real wtmp);

  void invert(Matrix<Real> &A, int nn);

 public:
  BroydenMixer(std::size_t size,
               Real alpha,
               std::size_t maxBroyden = 10,
               std::size_t iterReset = 25,
               Real w0 = 0.01);

  void mix(LSMSCommunication &comm, RealMixingVector &mix_vector) override;

  void mix(RealMixingVector &mix_vector) override;

};

/**
 * @brief
 */
void printMixingParameters(const MixingParameterPack &mix, LSMSCommunication &comm, LSMSSystemParameters &lsms);

}  // namespace lsms

#endif  // LSMS_MIXER_HPP
