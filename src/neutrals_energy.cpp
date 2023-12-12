// Copyright 2020, the Aether Development Team (see doc/dev_team.md
// for members) Full license can be found in License.md

#include "aether.h"

// ---------------------------------------------------------------------------
// Calculate thermal conduction
// ---------------------------------------------------------------------------

void Neutrals::update_temperature(Grid grid, Times time) {

  std::string function = "Neutrals::calc_conduction";
  static int iFunction = -1;
  report.enter(function, iFunction);

  precision_t dt;
  int64_t iLon, iLat;
  int64_t nLons = grid.get_nLons();
  int64_t nLats = grid.get_nLats();
  int64_t nAlts = grid.get_nAlts();
  int64_t nGCs = grid.get_nGCs();

  dt = time.get_dt();

  if (nAlts == 2 * nGCs + 1) {
    // Simply update the temperature if there is no conduction:
    temperature_scgc = temperature_scgc + dt * heating_sources_total;
    conduction_scgc.zeros();
  } else {

    arma_cube rhocvr23d(nLons, nLats, nAlts);
    arma_cube lambda3d(nLons, nLats, nAlts);
    arma_cube prandtl3d(nLons, nLats, nAlts);

    rhocvr23d = rho_scgc % Cv_scgc % grid.radius2_scgc;

    // Need to make this eddy * rho * cv:
    if (input.get_use_eddy_energy())
      prandtl3d = kappa_eddy_scgc % rho_scgc % Cv_scgc;
    else
      prandtl3d.zeros();

    lambda3d = (kappa_scgc + prandtl3d) % grid.radius2_scgc;

    arma_vec temp1d(nAlts);
    arma_vec lambda1d(nAlts);
    arma_vec rhocvr21d(nAlts);
    arma_vec dalt1d(nAlts);
    arma_vec sources1d(nAlts);
    arma_vec conduction1d(nAlts);

    for (iLon = 0; iLon < nLons; iLon++) {
      for (iLat = 0; iLat < nLats; iLat++) {

        temp1d = temperature_scgc.tube(iLon, iLat);
        lambda1d = lambda3d.tube(iLon, iLat);
        rhocvr21d = rhocvr23d.tube(iLon, iLat);
        sources1d = heating_sources_total.tube(iLon, iLat);

        temp1d = temp1d + dt * sources1d;
        sources1d.zeros();


        dalt1d = grid.dalt_lower_scgc.tube(iLon, iLat);
        conduction1d.zeros();

        conduction1d = solver_conduction(temp1d,
					 lambda1d,
					 rhocvr21d,
					 sources1d,
					 dalt1d,
					 dt,
					 nGCs,
					 false);
        temperature_scgc.tube(iLon, iLat) = conduction1d;

        // Store the difference (as a rate), so we can output it later
        // if we want:
        conduction_scgc.tube(iLon, iLat) =
	  (conduction1d - temp1d)/dt - sources1d;
      }  // lat
    }  // lon

  } // if nAlts == 1 + 2*GCs

  report.exit(function);
  return;
}
