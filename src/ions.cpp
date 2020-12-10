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

Ions::species_chars Ions::create_species() {

  long iDir, iLon, iLat, iAlt, index;
  species_chars tmp;

  long iTotal = long(nGeoLonsG) * long(nGeoLatsG) * long(nGeoAltsG);

  // Constants:
  tmp.DoAdvect = 0;
  
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

Ions::Ions(Inputs input, Report report) {

  long iTotal = long(nGeoLonsG) * long(nGeoLatsG) * long(nGeoAltsG);
  species_chars tmp;
  int iErr;

  report.print(2,"Initializing Ions");
  
  for (int iSpecies=0; iSpecies < nIons; iSpecies++) {
    tmp = create_species();
    species.push_back(tmp);
  }

  // State variables:
  density_s3gc = (float*) malloc( iTotal * sizeof(float) );
  velocity_v3gc = (float*) malloc( long(3)*iTotal * sizeof(float) );
  exb_v3gc = (float*) malloc( long(3)*iTotal * sizeof(float) );
  ion_temperature_s3gc = (float*) malloc( iTotal * sizeof(float) );
  electron_temperature_s3gc = (float*) malloc( iTotal * sizeof(float) );

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

	}
	
      }

      if (infile_ptr.eof()) IsDone = 1;

    }

    infile_ptr.close();

  }

  return iErr;
  
}

