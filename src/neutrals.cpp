// (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <math.h>
#include <iostream>
#include <fstream>
#include <armadillo>

#include "../include/constants.h"
#include "../include/inputs.h"
#include "../include/file_input.h"
#include "../include/neutrals.h"
#include "../include/ions.h"
#include "../include/grid.h"
#include "../include/report.h"
#include "../include/earth.h"

using namespace arma;

// -----------------------------------------------------------------------------
//  
// -----------------------------------------------------------------------------

Neutrals::species_chars Neutrals::create_species(Grid grid) {

  long iDir, iLon, iLat, iAlt, index;
  species_chars tmp;

  long nLons = grid.get_nLons();
  long nLats = grid.get_nLats();
  long nAlts = grid.get_nAlts();
  
  long iTotal = long(nLons) * long(nLats) * long(nAlts);

  // Constants:
  tmp.DoAdvect = 0;
  tmp.thermal_cond = 0.0;
  tmp.thermal_exp = 0.0;
  tmp.lower_bc_density = -1.0;

  tmp.density_scgc.set_size(nLons, nLats, nAlts);
  tmp.chapman_scgc.set_size(nLons, nLats, nAlts);
  tmp.scale_height_scgc.set_size(nLons, nLats, nAlts);
  tmp.ionization_scgc.set_size(nLons, nLats, nAlts);

  tmp.density_scgc.ones();
  tmp.chapman_scgc.ones();
  tmp.scale_height_scgc.ones();
  tmp.ionization_scgc.zeros();

  tmp.sources_scgc.set_size(nLons, nLats, nAlts);
  tmp.sources_scgc.zeros();
  tmp.losses.set_size(nLons, nLats, nAlts);
  tmp.losses_scgc.zeros();

  tmp.density_s3gc = (float*) malloc( iTotal * sizeof(float) );
  tmp.velocity_v3gc = (float*) malloc( long(3)*iTotal * sizeof(float) );
  tmp.chapman_s3gc = (float*) malloc( iTotal * sizeof(float) );
  tmp.ionization_s3gc = (float*) malloc( iTotal * sizeof(float) );

  for (iLon = 0; iLon < nLons; iLon++) {
    for (iLat = 0; iLat < nLats; iLat++) {
      for (iAlt = 0; iAlt < nAlts; iAlt++) {
	
	index = ijk_geo_s3gc(iLon,iLat,iAlt);

	tmp.density_s3gc[index] = 1.0e-32;
	tmp.chapman_s3gc[index] = 1.0e-32;
	tmp.ionization_s3gc[index] = 1.0e-32;

	for (iDir = 0; iDir < 3; iDir++) {
	  index = ijkl_geo_v3gc(iLon,iLat,iAlt,iDir);
	  tmp.velocity_v3gc[index] = 0.0;
	}
	
      }
    }
  }
	
  return tmp;
  
}

// -----------------------------------------------------------------------------
//  
// -----------------------------------------------------------------------------

Neutrals::Neutrals(Grid grid, Inputs input, Report report) {

  int iErr;
  species_chars tmp;

  long nLons = grid.get_nLons();
  long nLats = grid.get_nLats();
  long nAlts = grid.get_nAlts();
  
  report.print(2, "Initializing Neutrals");

  long iTotal = long(nLons) * long(nLats) * long(nAlts);

  for (int iSpecies=0; iSpecies < nSpecies; iSpecies++) {
    tmp = create_species(grid);
    neutrals.push_back(tmp);
  }

  // State variables:

  density_scgc.set_size(nLons, nLats, nAlts);
  density_scgc.ones();
  temperature_scgc.set_size(nLons, nLats, nAlts);
  temperature_scgc.ones();

  // Derived quantities:

  rho_scgc.set_size(nLons, nLats, nAlts);
  rho_scgc.ones();
  mean_major_mass_scgc.set_size(nLons, nLats, nAlts);
  mean_major_mass_scgc.ones();
  pressure_scgc.set_size(nLons, nLats, nAlts);
  pressure_scgc.ones();
  sound_scgc.set_size(nLons, nLats, nAlts);
  sound_scgc.ones();

  // Heating and cooling parameters:
  Cv_scgc.set_size(nLons, nLats, nAlts);
  Cv_scgc.ones();
  gamma_scgc.set_size(nLons, nLats, nAlts);
  gamma_scgc.zeros();
  kappa_scgc.set_size(nLons, nLats, nAlts);
  kappa_scgc.zeros();

  conduction_scgc.set_size(nLons, nLats, nAlts);
  heating_euv_scgc.set_size(nLons, nLats, nAlts);
  
  density_s3gc = (float*) malloc( iTotal * sizeof(float) );
  velocity_v3gc = (float*) malloc( long(3)*iTotal * sizeof(float) );
  temperature_s3gc = (float*) malloc( iTotal * sizeof(float) );

  // Source Terms:
  heating_euv_s3gc = (float*) malloc( iTotal * sizeof(float) );
  conduction_s3gc = (float*) malloc( iTotal * sizeof(float) );

  heating_efficiency = input.get_euv_heating_eff_neutrals();
  
  initial_temperatures = NULL;
  initial_altitudes = NULL;
  
  // This gets a bunch of the species-dependent characteristics:
  iErr = read_planet_file(input, report);

  // This specifies the initial conditions for the neutrals:
  iErr = initial_conditions(grid, input, report);

}

// -----------------------------------------------------------------------------
// Read in the planet file that describes the species - only neutrals
// -----------------------------------------------------------------------------

int Neutrals::read_planet_file(Inputs input, Report report) {

  int iErr = 0;
  std::string hash;
  std::ifstream infile_ptr;

  report.print(3,"In read_planet_file for Neutrals");

  infile_ptr.open(input.get_planet_species_file());

  if (!infile_ptr.is_open()) {
    std::cout << "Could not open input file: "
	      << input.get_planet_species_file() << "!!!\n";
    iErr = 1;
  } else {

    int IsDone = 0;

    while (!IsDone) {

      hash = find_next_hash(infile_ptr);
      
      if (report.test_verbose(4))
	std::cout << "hash : -->" << hash << "<--\n";

      if (hash == "#neutrals") {

	// Read in the characteristics as CSVs:
	report.print(4,"Found #neutrals!");
	
	std::vector<std::vector<std::string>> lines = read_csv(infile_ptr);

	// I should totally redo the initialization of the species,
	// since we could just do it here, but that is for the future.

	if (lines.size()-1 != nSpecies) {
	  std::cout << "number of neutrals species defined in planet.h file : "
		    << nSpecies << "\n";
	  std::cout << "number of species defined in planet.in file : "
		    << lines.size() << "\n";
	  std::cout << "These don't match!\n";
	  iErr = 1;
	} else {

	  // assume order of rows right now:
	  // name, mass, vibration, thermal_cond, thermal_exp, advect, lower BC
	  	  
	  for (int iSpecies=0; iSpecies < nSpecies; iSpecies++) {
	    report.print(5, "setting neutral species " + lines[iSpecies+1][0]);
	    neutrals[iSpecies].cName = lines[iSpecies+1][0];
	    neutrals[iSpecies].mass = stof(lines[iSpecies+1][1])*amu;
	    neutrals[iSpecies].vibe = stof(lines[iSpecies+1][2]);
	    neutrals[iSpecies].thermal_cond = stof(lines[iSpecies+1][3]);
	    neutrals[iSpecies].thermal_exp = stof(lines[iSpecies+1][4]);
	    neutrals[iSpecies].DoAdvect = stoi(lines[iSpecies+1][5]);
	    neutrals[iSpecies].lower_bc_density = stof(lines[iSpecies+1][6]);
	  }

	}
	
      }

      if (hash == "#temperature") {

	report.print(4,"Found #temperatures!");
	
	std::vector<std::vector<std::string>> temps = read_csv(infile_ptr);

	int nTemps = temps.size()-1;
	initial_temperatures = (float*) malloc(nTemps * sizeof(float) );
	initial_altitudes = (float*) malloc(nTemps * sizeof(float) );
	for (int iTemp=0; iTemp < nTemps; iTemp++) {
	  report.print(5, "reading initial temp alt " + temps[iTemp+1][0]);
	  // convert altitudes from km to m
	  initial_altitudes[iTemp] = stof(temps[iTemp+1][0]) * 1000;
	  initial_temperatures[iTemp] = stof(temps[iTemp+1][1]);
	}
	nInitial_temps = nTemps;
      }
      
      if (infile_ptr.eof()) IsDone = 1;

    }

    infile_ptr.close();

  }

  return iErr;
  
}

// -----------------------------------------------------------------------------
//  
// -----------------------------------------------------------------------------

float Neutrals::calc_scale_height(int iSpecies,
				  long index,
				  Grid grid) {

  float g = grid.gravity_s3gc[index];
  float t = temperature_s3gc[index];
  float m = neutrals[iSpecies].mass;
  float H = boltzmanns_constant * t / m / g;

  return H;

}

// -----------------------------------------------------------------------------
//  
// -----------------------------------------------------------------------------

int Neutrals::initial_conditions(Grid grid, Inputs input, Report report) {

  int iErr = 0;
  long iDir, iLon, iLat, iAlt, index, iA, indexm;
  float alt, r, H;

  report.print(3, "Creating Neutrals initial_condition");
  
  // ---------------------------------------------------------------------
  // This section assumes we want a hydrostatic solution given the
  // temperature profile in the planet.in file.
  // ---------------------------------------------------------------------

  long nLons = grid.get_nLons();
  long nLats = grid.get_nLats();
  long nAlts = grid.get_nAlts();
  long nGCs = grid.get_nGCs();

  // Let's assume that the altitudes are not dependent on lat/lon:
  
  fvec alt1d(nAlts);
  fvec temp1d(nAlts);

  fmat H2d(nLons,nLats);

  alt1d = grid.geoAlt_scgc.tube(0,0);

  if (nInitial_temps > 0) {
    for (iAlt = 0; iAlt < nAlts; iAlt++) {
      alt = alt1d(iAlt);
      // Find temperatures:
      if (alt <= initial_altitudes[0]) {
	temp1d[iAlt] = initial_temperatures[0];
      } else {
	if (alt >= initial_altitudes[nInitial_temps-1]) {
	  temp1d[index] = initial_temperatures[nInitial_temps-1];
	} else {
	  // Linear interpolation!
	  iA = 0;
	  while (alt > initial_altitudes[iA]) iA++;
	  iA--;
	  // alt will be between iA and iA+1:
	  r = (alt - initial_altitudes[iA]) /
	    (initial_altitudes[iA+1] - initial_altitudes[iA]);
	  temp1d[iAlt] =
	    (1.0-r) * initial_temperatures[iA  ] + 
	    (    r) * initial_temperatures[iA+1];
	}
      }
    }
  } else {
    temp1d = 200.0;
  }

  // spread the 1D temperature across the globe:
  for (iLon = 0; iLon < nLons; iLon++) {
    for (iLat = 0; iLat < nLats; iLat++) {
      temperature_scgc.tube(iLon,iLat) = temp1d;
    }
  }
    
  // Set the lower boundary condition:
  for (int iSpecies=0; iSpecies < nSpecies; iSpecies++) {

    neutrals[iSpecies].
      density_scgc.slice(0).
      fill(neutrals[iSpecies].lower_bc_density);

    // Integrate with hydrostatic equilibrium up:
    for (iAlt = 1; iAlt < nAlts; iAlt++) {
      neutrals[iSpecies].scale_height_scgc.slice(iAlt) =
	boltzmanns_constant *
	temperature_scgc.slice(iAlt) /
	( neutrals[iSpecies].mass *
	  grid.gravity_scgc.slice(iAlt));

      neutrals[iSpecies].density_scgc.slice(iAlt) =
	neutrals[iSpecies].density_scgc.slice(iAlt-1) %
	exp(-grid.dalt_lower_scgc.slice(iAlt) /
	    neutrals[iSpecies].scale_height_scgc.slice(iAlt));
      
    }
  }

  // Copy into old variable:

  for (iLon = 0; iLon < nLons; iLon++) {
    for (iLat = 0; iLat < nLats; iLat++) {
      for (iAlt = 0; iAlt < nAlts; iAlt++) {
	
	index = ijk_geo_s3gc(iLon,iLat,iAlt);
	temperature_s3gc[index] = temperature_scgc(iLon,iLat,iAlt);
	for (int iSpecies=0; iSpecies < nSpecies; iSpecies++)
	  neutrals[iSpecies].density_s3gc[index] = 
	    neutrals[iSpecies].density_scgc(iLon,iLat,iAlt);
      }
    }
  }
  
  calc_mass_density(report);
  return iErr;
  
}

//----------------------------------------------------------------------
// set_bcs
//----------------------------------------------------------------------

void Neutrals::set_bcs(Report &report) {

  std::string function="Neutrals::set_bcs";
  static int iFunction = -1;
  report.enter(function, iFunction);

  long nAlts = temperature_scgc.n_slices;

  // temperature_scgc.slice(nAlts-2).fill(800.0);
  // temperature_scgc.slice(nAlts-1).fill(800.0);

  temperature_scgc.slice(nAlts-2) = temperature_scgc.slice(nAlts-3);
  temperature_scgc.slice(nAlts-1) = temperature_scgc.slice(nAlts-2);

  report.exit(function);
  
}

//----------------------------------------------------------------------
// This code takes the EUV information that was read in from the EUV
// file and tries to figure out which things are absorbtion/ionization
// cross sections.  It does this by comparing the name of the neutral
// species to the first column in the euv.csv file.  If it finds a
// match, it then checks to see if it is an absorbtion or ionization
// cross section.  If it is an ionization cs, then it tries to figure
// out which ion it is producing (the "to" column).
// ---------------------------------------------------------------------

int Neutrals::pair_euv(Euv euv, Ions ions, Report report) {

  int iErr = 0;

  if (report.test_verbose(3))
    std::cout << "Pairing EUV abs/ion cross sections... (Neutrals::pair_euv) \n";
  
  for (int iSpecies=0; iSpecies < nSpecies; iSpecies++) {

    if (report.test_verbose(5))
      std::cout << neutrals[iSpecies].cName << "\n";

    neutrals[iSpecies].iEuvAbsId_ = -1;
    neutrals[iSpecies].nEuvIonSpecies = 0;
    
    // Check each row to see if the first column "name" matches:
    for (int iEuv=0; iEuv < euv.waveinfo.size(); iEuv++) {

      if (report.test_verbose(6))
	std::cout << "  " << euv.waveinfo[iEuv].name << "\n";

      // if this matches...
      if (neutrals[iSpecies].cName == euv.waveinfo[iEuv].name) {

	// First see if we can find absorbtion:
	if (euv.waveinfo[iEuv].type == "abs") {
	  if (report.test_verbose(5)) std::cout << "  Found absorbtion\n";
	  neutrals[iSpecies].iEuvAbsId_ = iEuv;
	}
	
	// Next see if we can find ionizations:
	if (euv.waveinfo[iEuv].type == "ion") {

	  // Loop through the ions to see if names match:
	  for (int iIon=0; iIon<nIons; iIon++) {
	    if (ions.species[iIon].cName == euv.waveinfo[iEuv].to) {
	      if (report.test_verbose(5))
		std::cout << "  Found ionization!! --> "
			  << ions.species[iIon].cName << "\n";
	      neutrals[iSpecies].iEuvIonId_.push_back(iEuv);
	      neutrals[iSpecies].iEuvIonSpecies_.push_back(iIon);
	      neutrals[iSpecies].nEuvIonSpecies++;
	    }
	  }
	}
      }
      
    }
    
  }
  
  return iErr;
  
}

