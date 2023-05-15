// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md



// This function calculates the coriolis acceleration vector for an object 
// in motion in a rotating frame of reference. It uses the following 
// formula: a = -2 * w * v where w is the rotational velocity of the 
// reference frame and v is the velocity of the object in the rotating frame  
#include "../include/aether.h"

std::vector<arma_cube> coriolis(std::vector<arma_cube> velocity, precision_t rotationrate) {
    std::vector<arma_cube> coriolis_vec(3);
    for (int64_t iComp = 0; iComp < 3; iComp++) {
        coriolis_vec[iComp] = -2 * rotationrate * velocity[iComp];
    }
    return coriolis_vec;
}