// (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>
#include <fstream>

#include "../include/inputs.h"
#include "../include/file_input.h"
#include "../include/constants.h"
#include "../include/sizes.h"
#include "../include/ions.h"
#include "../include/grid.h"
#include "../include/report.h"
#include "../include/earth.h"


// -----------------------------------------------------------------------------
//  
// -----------------------------------------------------------------------------

Ions::species_chars Ions::create_species(Grid grid) {

  long iDir, iLon, iLat, iAlt, index;
  species_chars tmp;

  long nLons = grid.get_nLons();
  long nLats = grid.get_nLats();
  long nAlts = grid.get_nAlts();
  
  long iTotal = long(nLons) * long(nLats) * long(nAlts);

  // Constants:
  tmp.DoAdvect = 0;

  tmp.density_scgc.set_size(nLons, nLats, nAlts);
  tmp.density_scgc.ones();
  tmp.temperature_scgc.set_size(nLons, nLats, nAlts);
  tmp.temperature_scgc.ones();
  tmp.ionization_scgc.set_size(nLons, nLats, nAlts);
  tmp.ionization_scgc.zeros();

  tmp.density_s3gc = (float*) malloc( iTotal * sizeof(float) );
  tmp.par_velocity_v3gc = (float*) malloc( long(3)*iTotal * sizeof(float) );
  tmp.perp_velocity_v3gc = (float*) malloc( long(3)*iTotal * sizeof(float) );
  tmp.temperature_s3gc = (float*) malloc( iTotal * sizeof(float) );
  tmp.ionization_s3gc = (float*) malloc( iTotal * sizeof(float) );

  for (iLon = 0; iLon < nGeoLonsG; iLon++) {
    for (iLat = 0; iLat < nGeoLatsG; iLat++) {
      for (iAlt = 0; iAlt < nGeoAltsG; iAlt++) {
	
	index = ijk_geo_s3gc(iLon,iLat,iAlt);

	tmp.density_s3gc[index] = 1.0;
	tmp.temperature_s3gc[index] = 1.0e-32;
	tmp.ionization_s3gc[index] = 1.0e-32;

	for (iDir = 0; iDir < 3; iDir++) {
	  index = ijkl_geo_v3gc(iLon,iLat,iAlt,iDir);
	  tmp.par_velocity_v3gc[index] = 0.0;
	  tmp.perp_velocity_v3gc[index] = 0.0;
	}
	
      }
    }
  }
	
  return tmp;
  
}

// -----------------------------------------------------------------------------
//  
// -----------------------------------------------------------------------------

Ions::Ions(Grid grid, Inputs input, Report report) {

  long nLons = grid.get_nLons();
  long nLats = grid.get_nLats();
  long nAlts = grid.get_nAlts();
  
  long iTotal = long(nLons) * long(nLats) * long(nAlts);

  species_chars tmp;
  int iErr;

  report.print(2,"Initializing Ions");
  
  for (int iSpecies=0; iSpecies < nIons; iSpecies++) {
    tmp = create_species(grid);
    species.push_back(tmp);
  }

  // Create one extra species for electrons
  tmp = create_species(grid);
  species.push_back(tmp);

  // State variables:
  density_s3gc = (float*) malloc( iTotal * sizeof(float) );
  velocity_v3gc = (float*) malloc( long(3)*iTotal * sizeof(float) );
  exb_v3gc = (float*) malloc( long(3)*iTotal * sizeof(float) );
  ion_temperature_s3gc = (float*) malloc( iTotal * sizeof(float) );
  electron_temperature_s3gc = (float*) malloc( iTotal * sizeof(float) );

  density_scgc.set_size(nLons, nLats, nAlts);
  density_scgc.ones();
  ion_temperature_scgc.set_size(nLons, nLats, nAlts);
  ion_temperature_scgc.ones();
  electron_temperature_scgc.set_size(nLons, nLats, nAlts);
  electron_temperature_scgc.ones();

  tmp.sources_scgc.set_size(nLons, nLats, nAlts);
  tmp.sources_scgc.zeros();
  tmp.losses.set_size(nLons, nLats, nAlts);
  tmp.losses_scgc.zeros();
  
  // This gets a bunch of the species-dependent characteristics:
  iErr = read_planet_file(input, report);

}

// -----------------------------------------------------------------------------
// Read in the planet file that describes the species - only ions
// -----------------------------------------------------------------------------

int Ions::read_planet_file(Inputs input, Report report) {

  int iErr = 0;
  std::string hash;
  std::ifstream infile_ptr;

  report.print(3,"In read_planet_file for Ions");

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

      if (hash == "#ions") {

	// Read in the characteristics as CSVs:
	report.print(4,"Found #ions!");
	
	std::vector<std::vector<std::string>> lines = read_csv(infile_ptr);

	// I should totally redo the initialization of the species,
	// since we could just do it here, but that is for the future.

	if (lines.size()-1 != nIons) {
	  std::cout << "number of ion species (nIons) defined in planet.h file : "
		    << nIons << "\n";
	  std::cout << "number of ions defined in planet.in file : "
		    << lines.size() << "\n";
	  std::cout << "These don't match!\n";
	  iErr = 1;
	} else {

	  // assume order of rows right now:
	  // name, mass, charge, advect
	  	  
	  for (int iSpecies=0; iSpecies < nIons; iSpecies++) {
	    report.print(5, "setting ion species " + lines[iSpecies+1][0]);
	    species[iSpecies].cName = lines[iSpecies+1][0];
	    species[iSpecies].mass = stof(lines[iSpecies+1][1])*amu;
	    species[iSpecies].charge = stoi(lines[iSpecies+1][2]);
	    species[iSpecies].DoAdvect = stoi(lines[iSpecies+1][3]);
	  }

	  species[nIons].cName = "e-";
	  species[nIons].mass = mass_electron;
	  species[nIons].charge = -1;
	  species[nIons].DoAdvect = 0;
	}
	
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

void Ions::fill_electrons(Grid grid,
			  Report &report) {

  long nLons, nLats, nAlts, iLon, iLat, iAlt, index;
  int iSpecies;
  float electron_density;
  
  std::string function = "Ions::fill_electrons";
  static int iFunction = -1;
  report.enter(function, iFunction);  

  species[nIons].density_scgc.zeros();
  for (iSpecies=0; iSpecies < nIons; iSpecies++)
    species[nIons].density_scgc = 
      species[nIons].density_scgc +
      species[iIon].density_scgc; 
  density_scgc = species[nIons].density_scgc;
  
  if (grid.get_IsGeoGrid()) {
    nLons = nGeoLonsG;
    nLats = nGeoLatsG;
    nAlts = nGeoAltsG;
  } else {
    nLons = nMagLonsG;
    nLats = nMagLatsG;
    nAlts = nMagAltsG;
  }

  for (iLon = 0; iLon < nLons; iLon++) {
    for (iLat = 0; iLat < nLats; iLat++) {
      for (iAlt = 0; iAlt < nAlts; iAlt++) {

	if (grid.get_IsGeoGrid()) {
	  index = ijk_geo_s3gc(iLon,iLat,iAlt);
	} else {
	  index = ijk_mag_s3gc(iLon,iLat,iAlt);
	}

	electron_density = 0.0;

	for (iSpecies=0; iSpecies < nIons; iSpecies++)
	  electron_density = electron_density  + species[iSpecies].density_s3gc[index];

	species[nIons].density_s3gc[index] = electron_density;
	density_s3gc[index] = electron_density;

      }
    }
  }

  report.exit(function);
  return;
}
