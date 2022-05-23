// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "../include/aether.h"

// --------------------------------------------------------------------------
// Initialize the electron temperature - set equal to the neutral temperature
// --------------------------------------------------------------------------

void Ions::init_electron_temperature(Neutrals neutrals, Grid grid,
                                     Report &report) {

  electron_temperature_scgc = neutrals.temperature_scgc;

}


// --------------------------------------------------------------------------
// Calculate the electron temperature
// --------------------------------------------------------------------------

void Ions::calc_electron_temperature(Neutrals neutrals, Grid grid, 
		                     Times time, Report &report) {

  std::string function = "Ions::calc_electron_temperature";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int64_t iLon, iLat, id_O, id_N2, id_O2;
  int64_t nLons = grid.get_nLons();
  int64_t nLats = grid.get_nLats();
  int64_t nAlts = grid.get_nAlts();

  arma_vec temp1d(nAlts);
  arma_vec lambda1d(nAlts);
  arma_vec front1d(nAlts);
  arma_vec dalt1d(nAlts);
  arma_vec conduction1d(nAlts);

  // Get the time step size
  precision_t dt = time.get_dt();

  // Loop over Lon,Lat positions
  for(iLon = 0; iLon < nLons; iLon++) {
    for(iLat = 0; iLat < nLats; iLat++) {

      // ----------------------------------------------------------------------
      // Calculate heat flux (conduction) in 1D; loop over all lat,lon
      // ----------------------------------------------------------------------
      temp1d = electron_temperature_scgc.tube(iLon, iLat);
      lambda1d = 1.60217646e-19 * 100.0 * 7.7e5 
        * pow(electron_temperature_scgc.tube(iLon, iLat),(5/2))
	/ (1.0 + 3.22e4 * pow(electron_temperature_scgc.tube(iLon, iLat),2) 
	/ density_scgc.tube(iLon, iLat) %
	  (neutrals.species[iO_].density_scgc.tube(iLon, iLat) * 1.1e-16
	% (1.0 + 5.7e-4 * electron_temperature_scgc.tube(iLon, iLat)) + 
	   neutrals.species[iN2_].density_scgc.tube(iLon, iLat) * 2.82e-17 
	% sqrt(electron_temperature_scgc.tube(iLon, iLat)) 
	% (1.0 - 1.21e-4 * electron_temperature_scgc.tube(iLon, iLat)) +
	   neutrals.species[iO2_].density_scgc.tube(iLon, iLat) * 2.2e-16 
	% (1.0 + 3.6e-2 * sqrt(electron_temperature_scgc.tube(iLon, iLat))) ) );
      front1d = 2.0 / 3.0 / cKB / density_scgc.tube(iLon, iLat); 
      dalt1d = grid.dalt_lower_scgc.tube(iLon, iLat);

      conduction1d.zeros();

      conduction1d = solver_conduction(temp1d, lambda1d, front1d, dt, dalt1d);

      // The conduction solver gives Tnew-Told, so divide by dt
      electron_conduction_scgc.tube(iLon, iLat) = conduction1d / dt;

    }  // iLat
  }  // iLon

  // Add temperature terms together to advance electron temperature
  electron_temperature_scgc = electron_temperature_scgc + 
	                      dt * (electron_conduction_scgc);

  report.exit(function);
  return;
}
