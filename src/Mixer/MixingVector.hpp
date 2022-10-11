//
// Created by F.Moitzi on 23.12.2022.
//

#ifndef LSMS_MIXINGVECTOR_HPP
#define LSMS_MIXINGVECTOR_HPP

#include <vector>

#include "Communication/LSMSCommunication.hpp"
#include "Main/SystemParameters.hpp"
#include "MixerType.hpp"
#include "Real.hpp"

namespace lsms {

/**
 * @brief Mixing Vector base class
 */
template<typename T>
class MixingVector {
 public:
  MixingVector() = default;

  std::vector<T> vec_old;
  std::vector<T> vec_new;
  T rms{};
  std::size_t size{};
  std::size_t total_size{};
  std::vector<std::size_t> starts;

  virtual void copyToVector(LSMSCommunication &comm, LSMSSystemParameters &lsms,
                            LocalTypeInfo &local) = 0;

  virtual void copyFromVector(LSMSCommunication &comm,
                              LSMSSystemParameters &lsms,
                              LocalTypeInfo &local) = 0;
};

using RealMixingVector = MixingVector<Real>;

/**
 * @brief Default mixing vector for testing purposed
 */
class DefaultMixerVector : public RealMixingVector {
 public:
  explicit DefaultMixerVector(std::size_t t);

  void copyToVector(LSMSCommunication &comm, LSMSSystemParameters &lsms,
                    LocalTypeInfo &local) override;

  void copyFromVector(LSMSCommunication &comm, LSMSSystemParameters &lsms,
                      LocalTypeInfo &local) override;
};

/**
 * @brief Potential mixing vector
 */
class PotentialMixingVector : public RealMixingVector {
 public:
  PotentialMixingVector(LSMSSystemParameters &lsms, LocalTypeInfo &local);

  void copyToVector(LSMSCommunication &comm, LSMSSystemParameters &lsms,
                    LocalTypeInfo &local) override;

  void copyFromVector(LSMSCommunication &comm, LSMSSystemParameters &lsms,
                      LocalTypeInfo &local) override;
};

/**
 * @brief Charge mixing vector
 */
class ChargeMixingVector : public RealMixingVector {

 private:

  bool no_spin_split;

 public:

  ChargeMixingVector() = default;

  ChargeMixingVector(LSMSSystemParameters &lsms, LocalTypeInfo &local, bool no_spin_split = false);

  void copyToVector(LSMSCommunication &comm, LSMSSystemParameters &lsms,
                    LocalTypeInfo &local) override;

  void copyFromVector(LSMSCommunication &comm, LSMSSystemParameters &lsms,
                      LocalTypeInfo &local) override;
};

/**
 * @brief Spin density mixing vector
 */
class SpinMixingVector : public RealMixingVector {
 public:

  SpinMixingVector() = default;

  SpinMixingVector(LSMSSystemParameters &lsms, LocalTypeInfo &local);

  void copyToVector(LSMSCommunication &comm, LSMSSystemParameters &lsms,
                    LocalTypeInfo &local) override;

  void copyFromVector(LSMSCommunication &comm, LSMSSystemParameters &lsms,
                      LocalTypeInfo &local) override;
};

}  // namespace lsms

#endif  // LSMS_MIXINGVECTOR_HPP
