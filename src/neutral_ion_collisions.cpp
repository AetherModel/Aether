// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "../include/aether.h"

void Neutrals::calc_ion_collisions(Ions &ions, Grid &grid, precision_t dt, Report &report) {
    
    std::string function = "Neutrals::calc_ion_collisions";
    static int iFunction = -1;
    report.enter(function, iFunction);
    
    int64_t nX = grid.get_nX();
    int64_t nY = grid.get_nY();
    int64_t nZ = grid.get_nZ();
    
    arma_cube rho_n(nX, nY, nZ);
    arma_cube rho_i(nX, nY, nZ);
    arma_cube rho_sum(nX, nY, nZ);
    
    // Calculate acceleration due to ion drag. Based on Formula 4.124b in Ionospheres text.
    for (int64_t iNeutral = 0; iNeutral < nSpecies; iNeutral++) {
        rho_n = species[iNeutral].mass * species[iNeutral].density_scgc;
        
        for (int64_t iDir = 0; iDir < 3; iDir++) {
            for (int64_t iIon = 0; iIon < ions.nSpecies; iIon++) {
                rho_i = ion.species[iIon].mass * ion.species[iIon].density_scgc;
                rho_sum += rho_i % ions.species[iIon].nu_ion_neutral_coef[iSpecies]
                    % (ions.species[iIon].par_velocity_vcgc[iDir] + ions.species[iIon].perp_velocity_vcgc[iDir] - species[iNeutral].velocity_vcgc[iDir]);
            } // for each ion
            species[iNeutral].acc_ion_drag[iDir] = rho_sum / rho_n;
        } // for each direction
    } // for each neutral
    
    report.exit(function);
    
} // calc_ion_collisions
