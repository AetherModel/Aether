// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>

#include "aether.h"

//returns new acceleration vector:
arma_vec neutral_friction_one_cell(int64_t iLong, int64_t iLat, int64_t iAlt,
                               precision_t dt,
                               arma_vec &vels,
                               Neutrals &neutrals,
                               Report &report) {
    
    std::string function = "neutral_friction_one_cell";
    static int iFunction = -1;
    report.enter(function, iFunction);
    
    arma_mat matrix(neutrals.nSpecies, neutrals.nSpecies, fill::zeros);
    
    for (int iSpec = 0; iSpec < neutrals.nSpecies; iSpec++) {
        
        // ktom = boltzmann's constant * temperature / mass
        precision_t ktom = cKB * temperature_scgc(iLong, iLat, iAlt) / neutrals.species[iSpec].mass;
        
        for (int jSpec = 0; jSpec < neutrals.nSpecies; jSpec++) {
            
            if (iSpec == jSpec) continue;
            
            // temp_dij holds the dij binary coefficients based upon the formulation by Banks and Kokarts.
            // These coefficients demand that:
            // (1) density be in cm^-3 (hence the 1.0e-06) factor below
            // (2) Additionally, the Dij's are in cm^2/s, thus the 1.0e-04 factor
            precisiont_t temp_dij = 1e-4 * neutrals.species[iSpec].diff0[jSpec] *
                pow(neutrals.temperature_scgc(iLong, iLat, iAlt), neutrals.species[iSpec].diff_exp[jSpec]) /
                (neutrals.density_ncgc(iLong, iLat, iAlt) * 1e-6);
            
            matrix(iSpec, jSpec) = ktom * neutrals.species[jSpec].density_ncgc(iLong, iLat, iAlt) /
                (temp_dij * neutrals.density_ncgc(iLong, iLat, iAlt));
            
        } // jSpec loop
    } // iSpec loop
    
    matrix = -1 * dt * matrix;
    
    // Fill in diagonal of matrix:
    for (int iSpec = 0; iSpec < neutrals.nSpecies; iSpec++) {
        matrix(iSpec, iSpec) = 1 + dt * sum(matrix.row(iSpec));
    }
    
    
    // initialize array of each neutral species' accelerations at (iLong, iLat, iAlt):
    arma_vec accs(neutrals.nSpecies, fill::zeros);
    
    // Solve system of equations:
    arma_vec new_vels = arma::solve(matrix, vels);
    
    // put the new values into the velocity cubes:
    for (int64_t iSpec = 0; iSpec < neutrals.nSpecies; iSpec++) {
        accs(iSpec) = (new_vels(iSpec) - vels(iSpec))/dt;
    }
    
    return accs;
    
    report.exit(function);
} // neutral_friction_one_cell


    
void calc_neutral_friction(precision_t dt,
                           Grid &gGrid,
                           Neutrals &neutrals,
                           Report &report) {
            
    std::string function = "calc_neutral_friction";
    static int iFunction = -1;
    report.enter(function, iFunction);
            
    for (int64_t iAlt = 0; iAlt < gGrid.get_nZ(); iAlt++) {
        for (int64_t iLat = 0; iLat < gGrid.get_nY(); iLat++) {
            for (int64_t iLong = 0; iLong < gGrid.get_nX(); iLong++) {
                
                for (int iDir = 0; iDir < 3; iDir++) {
                    
                    arma_vec vels(neutrals.nSpecies, fill::zeros);
                    for (int64_t iSpec = 0; iSpec < neutrals.nSpecies; iSpec++) {
                        //Put the old velocities into vels:
                        vels(iSpec) = neutrals.species[iSpec].velocity_vcgc[iDir](iLong, iLat, iAlt);
                    } // for species
                    arma_vec acc = neutral_friction_one_cell(iLong, iLat, iAlt, dt, vels, neutrals neutrals, report);
                    for (int64_t iSpec = 0; iSpec < neutrals.nSpecies; iSpec++) {
                        neutrals.species[iSpec].acc_vcgc[iDir](iLong, iLat, iAlt) += acc(iSpec);
                    } // for species
                } // for direction
            } // for long
        } // for lat
    } // for alt

    
    report.exit(function);
    
    
} //calc_neutral_friction
