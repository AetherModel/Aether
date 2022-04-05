// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <math.h>
#include <iostream>
#include <fstream>

#include "aether.h"


// --------------------------------------------------------------------------
// Initialize the ion temperature - set equal to the neutral temperature
// --------------------------------------------------------------------------

void Ions::init_ion_temperature(Neutrals neutrals, Grid grid, Report &report) {

  int64_t iIon;

  for (iIon = 0; iIon < nIons; iIon++)
    species[iIon].temperature_scgc = neutrals.temperature_scgc;

  temperature_scgc = neutrals.temperature_scgc;

  return;
}


// --------------------------------------------------------------------------
// Calculate the ion temperature - Option) one temperature for all species
// --------------------------------------------------------------------------

void Ions::calc_ion_temperature(Neutrals neutrals, Grid grid, 
                                    Times time, Report &report) {

  int64_t iIon, iLon, iLat;   
  int64_t nLons = grid.get_nLons();
  int64_t nLats = grid.get_nLats();
  int64_t nAlts = grid.get_nAlts();
  
  std::string function = "Ions::calc_ion_temperature";
  static int iFunction = -1;
  report.enter(function, iFunction);

  arma_cube rhocvr23d(nLons, nLats, nAlts);
  arma_cube lambda3d(nLons, nLats, nAlts);
  arma_cube prandtl3d(nLons,nLats, nAlts);
  rhocvr23d = neutrals.rho_scgc % neutrals.Cv_scgc % grid.radius2_scgc;
  prandtl3d.zeros();
  lambda3d = (neutrals.kappa_scgc + prandtl3d) % grid.radius2_scgc;

  arma_vec temp1d(nAlts);
  arma_vec lambda1d(nAlts);
  arma_vec rhocvr21d(nAlts);
  arma_vec dalt1d(nAlts);
  arma_vec conduction1d(nAlts);

  // Get the time step size
  precision_t dt = time.get_dt();

  // --------------------------------------------------------------------------
  // Calculate heat flux (conduction) in 1D; loop over all lat,lon positions
  // --------------------------------------------------------------------------
  for (iLon = 0; iLon < nLons; iLon++) {
    for (iLat = 0; iLat < nLats; iLat++) {

      temp1d    = temperature_scgc.tube(iLon, iLat);     // current bulk ion temp
      lambda1d  = lambda3d.tube(iLon, iLat);             // lambda
      rhocvr21d = rhocvr23d.tube(iLon, iLat);            // ???
      dalt1d    = grid.dalt_lower_scgc.tube(iLon, iLat); // grid thing for solver

      conduction1d.zeros();                              // set temp variable to zero

      conduction1d = solver_conduction(temp1d, lambda1d, rhocvr21d, dt, dalt1d);

      // The conduction solver gives Tnew-Told, so divide by dt
      conduction_scgc.tube(iLon, iLat) = conduction1d / dt;

    }  // lat
  }  // lon


  // --------------------------------------------------------------------------
  // Add temperature terms together to advance ion temperature
  // --------------------------------------------------------------------------
  temperature_scgc = temperature_scgc + dt * (conduction_scgc);


  // --------------------------------------------------------------------------
  // Use the bulk ion temperature for all ion specie temperatures
  // --------------------------------------------------------------------------
  for (iIon = 0; iIon < nIons; iIon++)
    species[iIon].temperature_scgc = temperature_scgc;
   

  report.exit(function);
  return;
}
