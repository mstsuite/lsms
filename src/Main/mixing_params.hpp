
#ifndef LSMS_MIXING_PARAMS_HPP
#define LSMS_MIXING_PARAMS_HPP

#include "Real.hpp"

struct MixingParameters {

    // Different mixing quantities and algorithms
    static const int numQuantities = 5;

    enum mixQuantity {
        no_mixing = 0, charge = 1, potential = 2, moment_magnitude = 3,
        moment_direction = 4
    };
    enum mixAlgorithm {
        noAlgorithm = 0, simple = 1, broyden = 2
    };

    // These parameters specify the which quantity(ies) is (are) being mixed and which algorithm(s) to used.
    // The correspondances of the indices are specified in mixQuantity.
    // bool values:
    // 0 : quantity is not used for mixing
    // 1 : quantity is used for mixing
    bool quantity[numQuantities];
    mixAlgorithm algorithm[numQuantities];
    Real mixingParameter[numQuantities];

};



#endif //LSMS_MIXING_PARAMS_HPP
