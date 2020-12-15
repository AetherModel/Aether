// (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>

#include "../include/sizes.h"

#include "../include/chemistry.h"
#include "../include/neutrals.h"
#include "../include/ions.h"
#include "../include/times.h"
#include "../include/grid.h"
#include "../include/report.h"
#include "../include/solvers.h"

void Chemistry::calc_chemistry(Neutrals &neutrals,
			       Ions &ions,
			       Times time,
			       Grid grid,
			       Report &report) {

  long nLons, nLats, nAlts, iLon, iLat, iAlt, index;
  int iSpecies;

  std::string function = "Chemistry::calc_chemistry";
  report.enter(function);

  if (grid.get_IsGeoGrid()) {
    nLons = nGeoLonsG;
    nLats = nGeoLatsG;
    nAlts = nGeoAltsG;
  } else {
    nLons = nMagLonsG;
    nLats = nMagLatsG;
    nAlts = nMagAltsG;
  }

  float neutral_density[nSpecies];
  float ion_density[nIons+1];  // add one for electron density
  float neutral_sources[nSpecies];
  float neutral_losses[nSpecies];
  float ion_sources[nIons];
  float ion_losses[nIons];

  float Tn, Te, Ti;

  float dt = time.get_dt();
  float old_density, source, loss, new_density;
  float heat_neutrals, heat_ions, heat_electrons;
  
  // ------------------------------------
  // Calculate electron densities
  // ------------------------------------

  ions.fill_electrons(grid, report);

  // Don't do chemistry in the ghostcells!

  for (iLon = 0; iLon < nLons; iLon++) {
    for (iLat = 0; iLat < nLats; iLat++) {
      for (iAlt = 0; iAlt < nAlts; iAlt++) {

	if (grid.get_IsGeoGrid()) {
	  index = ijk_geo_s3gc(iLon,iLat,iAlt);
	} else {
	  index = ijk_mag_s3gc(iLon,iLat,iAlt);
	}

	// Use the private variable sources_and_losses, so we don't
	// have to pass it around

	for (iSpecies=0; iSpecies < nSpecies; iSpecies++) {
	  neutral_density[iSpecies] = neutrals.neutrals[iSpecies].density_s3gc[index];
	  sources_and_losses.neutral_losses[iSpecies] = neutrals.neutrals[iSpecies].ionization_s3gc[index];
	  sources_and_losses.neutral_sources[iSpecies] = 0.0;
	}

	for (iSpecies=0; iSpecies < nIons; iSpecies++) {
	  ion_density[iSpecies] = ions.species[iSpecies].density_s3gc[index];
	  sources_and_losses.ion_sources[iSpecies] = ions.species[iSpecies].ionization_s3gc[index];
	  sources_and_losses.ion_losses[iSpecies] = 0.0;
	}
	ion_density[nIons] = ions.density_s3gc[index];

	Tn = neutrals.temperature_s3gc[index];
	Ti = ions.ion_temperature_s3gc[index];
	Te = ions.electron_temperature_s3gc[index];

	// If we wanted to do a higher-order solver, we would probably
	// put it starting here: (If we do that, we may want to have a
	// code that calculates the reaction rates first, then take
	// this outside of the higher-order solver, since the reaction
	// rates involve a lot of powers, which seem like there are
	// quite slow in C.)
	
	calc_chemical_sources(neutral_density,
			      ion_density,
			      Tn, Ti, Te, report);
	
	//report.enter("solver_chemistry");
	for (iSpecies=0; iSpecies < nSpecies; iSpecies++) {
	  old_density = neutral_density[iSpecies];
	  source = sources_and_losses.neutral_sources[iSpecies];
	  loss = sources_and_losses.neutral_losses[iSpecies];
	  neutrals.neutrals[iSpecies].density_s3gc[index] =
	    solver_chemistry(old_density, source, loss, dt);
	}

	for (iSpecies=0; iSpecies < nIons; iSpecies++) {
	  old_density = ion_density[iSpecies];
	  source = sources_and_losses.ion_sources[iSpecies];
	  loss = sources_and_losses.ion_losses[iSpecies];
	  ions.species[iSpecies].density_s3gc[index] =
	    solver_chemistry(old_density, source, loss, dt);
	}
	//report.exit("solver_chemistry");
	
      }
    }
  }
  
  report.exit(function);
  return;
}
