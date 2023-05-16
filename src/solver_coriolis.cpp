// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md



// This function calculates the coriolis acceleration vector for an object 
// in motion in a rotating frame of reference. It uses the following 
// formula: a = 2 * W x v where W is the rotational velocity of the 
// reference frame and v is the velocity of the object in the rotating frame 
// This assumes the grid is lined up with the tilt 
#include "../include/aether.h"

std::vector<arma_cube> coriolis(std::vector<arma_cube> velocity, precision_t rotation_rate) {
    std::vector<arma_cube> coriolis_vec(3);
    std::vector<arma_cube> cross;

    // Rotational velocity vector
    std::vector<int> rotational_velocity{ 0, 0, rotation_rate};

    // Cross Product of rotational velocity and velocity:
    // Latitude
    cross.push_back(velocity[1] * rotational_velocity[2] - velocity[2] * rotational_velocity[1])
    // Longitude
    cross.push_back(velocity[0] * rotational_velocity[2] - velocity[2] * rotational_velocity[0])
    // Altitude
    cross.push_back(velocity[0] * rotational_velocity[1] - velocity[1] * rotational_velocity[0])

    for (int64_t iComp; iComp < 3; iComp++) {
        coriolis_vec[iComp] = 2 * cross[iComp]
    }

    return coriolis_vec;
}