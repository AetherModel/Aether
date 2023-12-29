// Copyright 2023, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md
//
// initial version - A. Ridley - May 27, 2023

#include "aether.h"

// -----------------------------------------------------------------------------
//  Set initial conditions for the neutrals.
//    Two methods implemented so far:
//      - Planet: Use fixed density values in the planet.in file and the
//                temperature profile to set the densities and temperature.
//                Densities are filled with hydrostatic solution.
//      - Msis: Use NRL MSIS to set the densities and temperatures.  If the
//              densities are not found, then set to density in planet.in
//              file and fill with hydrostatic.
// -----------------------------------------------------------------------------

//----------------------------------------------------------------------
// set_bcs - This is for setting the vertical BCs
//----------------------------------------------------------------------

bool Ions::set_bcs(Grid grid,
                   Times time,
                   Indices indices) {

    std::string function = "Ions::set_bcs";
    static int iFunction = -1;
    report.enter(function, iFunction);

    bool didWork = true;

    if (input.get_nAltsGeo() > 1) {
        didWork = set_lower_bcs(grid, time, indices);

        if (didWork)
            didWork = set_upper_bcs(grid);

        if (didWork)
            fill_electrons();
    }

    if (!didWork)
        report.error("issue with ion BCs!");

    report.exit(function);
    return didWork;
}

//----------------------------------------------------------------------
// set upper boundary conditions for the ions
//----------------------------------------------------------------------

bool Ions::set_upper_bcs(Grid grid) {

    std::string function = "Ions::set_upper_bcs";
    static int iFunction = -1;
    report.enter(function, iFunction);

    bool didWork = true;

    int64_t nAlts = grid.get_nZ();
    int64_t nX = grid.get_nX(), iX;
    int64_t nY = grid.get_nY(), iY;
    int64_t nGCs = grid.get_nGCs();
    int64_t iAlt;
    arma_mat h;
    arma_mat aveT;

    for (iAlt = nAlts - nGCs; iAlt < nAlts; iAlt++) {
        // Bulk Quantities:
        temperature_scgc.slice(iAlt) = temperature_scgc.slice(iAlt - 1);

        // For each species:
        for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
            species[iSpecies].temperature_scgc.slice(iAlt) =
                species[iSpecies].temperature_scgc.slice(iAlt - 1);

            aveT = (species[iSpecies].temperature_scgc.slice(iAlt) + 
                    electon_temperature_scgc.slice(iAlt));
            // Calculate scale height for the species:
            h = cKB * species[iSpecies].temperature_scgc.slice(iAlt) /
                (species[iSpecies].mass % abs(grid.gravity_vcgc[2]));
            // Assume each species falls of with (modified) hydrostatic:
            species[iSpecies].density_scgc.slice(iAlt) =
                species[iSpecies].density_scgc.slice(iAlt - 1) %
                exp(-grid.dalt_lower_scgc.slice(iAlt) / h);
        }
    }

    report.exit(function);
    return didWork;
}

//----------------------------------------------------------------------
// set lower boundary conditions for the ions
//----------------------------------------------------------------------

bool Ions::set_lower_bcs(Grid grid) {

    std::string function = "Ions::set_lower_bcs";
    static int iFunction = -1;
    report.enter(function, iFunction);

    bool didWork = true;

    int64_t nAlts = grid.get_nZ();
    int64_t nX = grid.get_nX(), iX;
    int64_t nY = grid.get_nY(), iY;
    int64_t nGCs = grid.get_nGCs();
    int64_t iAlt;
    arma_mat h;
    arma_mat aveT;

    for (iAlt = nGCs - 1; iAlt >= 0; iAlt--) {
        // Bulk Quantities:
        temperature_scgc.slice(iAlt) = temperature_scgc.slice(iAlt + 1);

        // For each species:
        for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
            // assign all species temperatures the bulk temperature:
            species[iSpecies].temperature_scgc.slice(iAlt) =
                temperature_scgc.slice(iAlt);
            // Assume each species falls off a bit.
            // this BC shouldn't matter, since the bottom of the code
            // should be in chemical equalibrium:
            species[iSpecies].density_scgc.slice(iAlt) =
                0.95 * species[iSpecies].density_scgc.slice(iAlt + 1);
        }
    }

    report.exit(function);
    return didWork;
}

