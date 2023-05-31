// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md



// This function calculates the coriolis acceleration vector for an object 
// in motion in a rotating frame of reference. It uses the coriolis terms in equations
// (14), (18), (19) in https://drive.google.com/file/d/1Q0cMzKhdd0IoXzBl3odTge9JhBwfAyk9/view?usp=share_link

#include "../include/aether.h"

std::vector<arma_cube> coriolis(std::vector<arma_cube> velocity, precision_t rotation_rate, arma_cube lat_scgc) {
    std::vector<arma_cube> coriolis_vec(3);
    // std::vector<arma_cube> cross;

    // Rotational velocity vector
    // std::vector<int> rotational_velocity{ rotation_rate * cos(geoLon_scgc) * geoLat_scgc , 0, 0};
    // std::vector<arma_cube> rotational_velocity(3);
    // rotational_velocity[0].zeros();
    // rotational_velocity[1].zeros();
    // rotational_velocity[2].fill(rotation_rate);

    coriolis_vec[0] = -2 * rotation_rate * velocity[1] % sin(lat_scgc);
    coriolis_vec[1] = 2 * rotation_rate * velocity[0] % sin(lat_scgc) - 2 * rotation_rate * velocity[2] % cos(lat_scgc);
    coriolis_vec[2] = 2 * rotation_rate * cos(lat_scgc) % velocity[1];
    // Cross Product of rotational velocity and velocity:
    // Latitude
    // cross.push_back(velocity[1] * rotational_velocity[2] - velocity[2] * rotational_velocity[1]);
    // Longitude
    // cross.push_back(velocity[0] * rotational_velocity[2] - velocity[2] * rotational_velocity[0]);
    // Altitude
    // cross.push_back(velocity[0] * rotational_velocity[1] - velocity[1] * rotational_velocity[0]);

    
    // for (int64_t iComp; iComp < 3; iComp++) {
    //     coriolis_vec[iComp] = 2 * cross[iComp];
        // coriolis_vec[iComp] = 2 * cross_product(velocity, rotational_velocity)
    //}

    return coriolis_vec;
}