//
// Created by F.Moitzi on 15.05.2022.
//

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <numeric>

#include "Mixer/Mixer.hpp"
#include "Mixer/MixerType.hpp"
#include "Mixer/MixingVector.hpp"
#include "accel_common.hpp"

extern "C" {

void vector_fun(double *x, double *f);
}

namespace mixer_tests {

template <class T>
auto squareError = [](T a, T b) {
  auto e = a - b;
  return e * e;
};

template <class T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &v) {
  for (auto ii = v.begin(); ii != v.end(); ++ii) {
    os << " " << std::setprecision(12) << std::fixed << *ii;
  }
  return os;
}

TEST(MixerTests, FixPointIterationTests) {
  LSMSCommunication comm;

  auto MAX_ITER = 10000;

  initializeCommunication(comm);

  lsms::DefaultMixerVector mixingVector(5);

  lsms::MixingParameter start_parameter{0.03, 0.01, 8, 12};

  // 0. Test creation

  std::unique_ptr<lsms::AbstractMixer> BroydenMixer =
      lsms::AbstractMixer::generateMixer(lsms::MixerType::BROYDEN_MIXER,
                                         mixingVector, start_parameter);

  mixingVector.vec_old = {0.1, 2.0, 1.0, 3.0, 0.5};

  std::vector<Real> prev_new = mixingVector.vec_old;

  int i = 1;

  for (i = 1; i <= MAX_ITER; i++) {
    // Eval function
    vector_fun(mixingVector.vec_old.data(), mixingVector.vec_new.data());

    mixingVector.rms = std::sqrt(std::transform_reduce(
        mixingVector.vec_new.begin(), mixingVector.vec_new.end(),
        mixingVector.vec_old.begin(), 0.0, std::plus<>(), squareError<Real>));

    std::cout << " --- " << i << " --- " << std::endl;
    std::cout << mixingVector.rms << std::endl;
    std::cout << mixingVector.vec_old << std::endl;

    prev_new = mixingVector.vec_new;

    std::cout << mixingVector.vec_new << std::endl;

    if (mixingVector.rms < 1.0e-12) {
      break;
    }

    // Mix
    BroydenMixer->mix(comm, mixingVector);

    std::cout << mixingVector.vec_new << std::endl;

    mixingVector.vec_old = mixingVector.vec_new;
  }

  ASSERT_TRUE(i < MAX_ITER);

  finalizeCommunication();
}

}  // namespace mixer_tests