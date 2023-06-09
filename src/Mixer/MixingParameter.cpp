//
// Created by F.Moitzi on 21.12.2022.
//

#include "MixingParameter.hpp"

namespace lsms {

MixingParameter::MixingParameter(double alpha, double w0,
                                 unsigned int max_broyden,
                                 unsigned int iter_reset
)
    : alpha{alpha}, w0{w0}, max_broyden{max_broyden}, iter_reset{iter_reset} {}
}  // namespace lsms