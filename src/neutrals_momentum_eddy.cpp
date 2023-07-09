// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>

#include "aether.h"

void Neutrals::vertical_momentum_eddy(Grid &gGrid, Report &report, Input inputs) {
    
    std::string function = "Neutrals::vertical_momentum_eddy";
    static int iFunction = -1;
    report.enter(function, iFunction);
    
    arma_cube log_cons;
    arma_cube grad_cons;
    
    for (int64_t iSpecies = 0; iSpecies < neutrals.nSpecies; iSpecies++) {
        
        // Take the natural log of each concentration cube:
        log_cons = log(species[iSpecies].concentration_scgc);
        
        // calculate gradient:
        grad_cons = calc_gradient_alt(log_cons, gGrid);
        
        species[iSpec].acc_eddy += -grad_cons * inputs.get_eddy_coef();
        
        
    }
    
    report.exit(function);
} // neutral_friction_momentum_eddy
