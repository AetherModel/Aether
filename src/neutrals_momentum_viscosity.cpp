// Copyright 2020, the Aether Development Team (see doc/dev_team.md
// for members) Full license can be found in License.md

#include "aether.h"

// ---------------------------------------------------------------------------
// Calculate viscosity
// ---------------------------------------------------------------------------

void Neutrals::update_horizontal_velocity(Grid grid, Times time) {

  std::string function = "Neutrals::update_horizontal_velocity";
  static int iFunction = -1;
  report.enter(function, iFunction);

  precision_t dt;
  int64_t iLon, iLat, iDir;
  int64_t nLons = grid.get_nLons();
  int64_t nLats = grid.get_nLats();
  int64_t nAlts = grid.get_nAlts();
  int64_t nGCs = grid.get_nGCs();

  dt = time.get_dt();

  if (nAlts == 2 * nGCs + 1) {
    // Simply update the temperature if there is no conduction:
    velocity_vcgc[0] = velocity_vcgc[0] + dt * acc_sources_total[0];
    velocity_vcgc[1] = velocity_vcgc[1] + dt * acc_sources_total[1];
  } else {

    arma_cube rhor23d(nLons, nLats, nAlts);
    arma_cube lambda3d(nLons, nLats, nAlts);
    //arma_cube prandtl3d(nLons, nLats, nAlts);

    rhor23d = rho_scgc % grid.radius2_scgc;

    //// Need to make this eddy * rho * cv:
    //if (input.get_use_eddy_energy())
    //  prandtl3d = kappa_eddy_scgc % rho_scgc % Cv_scgc;
    //else
    //  prandtl3d.zeros();

    lambda3d = (viscosity_scgc) % grid.radius2_scgc;

    arma_vec vel1d(nAlts);
    arma_vec lambda1d(nAlts);
    arma_vec rhor21d(nAlts);
    arma_vec dalt1d(nAlts);
    arma_vec sources1d(nAlts);
    arma_vec visc1d(nAlts);

    //if (iProc == 1) {
    //  std::cout << "neutrals_energy heating source : ";
    //  display_vector(heating_sources_total.tube(11,2));
    //  std::cout << "  -> temp before: ";
    //  display_vector(temperature_scgc.tube(11,2));
    //}
    
    for (iDir = 0; iDir < 2; iDir++) {
      for (iLon = 0; iLon < nLons; iLon++) {
	for (iLat = 0; iLat < nLats; iLat++) {

	  vel1d = velocity_vcgc[iDir].tube(iLon, iLat);
	  lambda1d = lambda3d.tube(iLon, iLat);
	  rhor21d = rhor23d.tube(iLon, iLat);
	  sources1d.zeros();
	  dalt1d = grid.dalt_lower_scgc.tube(iLon, iLat);
	  visc1d.zeros();

	  visc1d = solver_conduction(vel1d,
				     lambda1d,
				     rhor21d,
				     sources1d,
				     dalt1d,
				     dt,
				     nGCs,
				     false);
	  velocity_vcgc[iDir].tube(iLon, iLat) = visc1d;

	  // Store the difference (as a rate), so we can output it later
	  // if we want:
	  //conduction_scgc.tube(iLon, iLat) =
	  //  (conduction1d - temp1d)/dt - sources1d;
	}  // lat
      }  // lon
    }

    //if (iProc == 1) {
    //  std::cout << "  -> temp after: ";
    //  display_vector(temperature_scgc.tube(11,2));
    //}
    
  } // if nAlts == 1 + 2*GCs

  report.exit(function);
  return;
}
