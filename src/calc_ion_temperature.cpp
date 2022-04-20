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
// Calculate the ion temperature 
// --------------------------------------------------------------------------

void Ions::calc_ion_temperature(Neutrals neutrals, Grid grid, 
                                Times time, Inputs input, Report &report) {
  
  std::string function = "Ions::calc_ion_temperature";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int64_t iIon, iLon, iLat, nSpecs;   
  int64_t nLons = grid.get_nLons();
  int64_t nLats = grid.get_nLats();
  int64_t nAlts = grid.get_nAlts();

  arma_vec temp1d(nAlts);
  arma_vec den1d(nAlts);
  arma_vec lambda1d(nAlts);
  arma_vec front1d(nAlts);
  arma_vec dalt1d(nAlts);
  arma_vec conduction1d(nAlts);

  arma_cube tempT(nLons, nLats, nAlts);

  // Get the time step size
  precision_t dt = time.get_dt();
  std::cout << "Test ion temp 1 at Alt: 185, Lat: 47.5, Lon: 175 ---> " 
	    << temperature_scgc(20,30,20) << ", " << species[0].temperature_scgc(20,30,20) << "\n";
  //std::cout << "Test neutral temp too " << neutrals.temperature_scgc(20,30,20) << "\n";
  // Loop over all species or assume only bulk calculation
  if (input.get_do_calc_bulk_ion_temp()==true)  nSpecs=1; // First ion specie only, is O+?
  if (input.get_do_calc_bulk_ion_temp()==false) nSpecs=nIons;
  //std::cout << "Bulk ion temp flag: " << input.get_do_calc_bulk_ion_temp() 
  //	    << " so number of ions is " << nSpecs << "\n";

  // Loop over all species or assume only bulk calculation
  for (iIon = 0; iIon < nSpecs; iIon++) {

    for (iLon = 0; iLon < nLons; iLon++) {
      for (iLat = 0; iLat < nLats; iLat++) {

        // --------------------------------------------------------------------------
        // Calculate heat flux (conduction) in 1D; loop over all lat,lon positions
        // --------------------------------------------------------------------------
        temp1d   = species[iIon].temperature_scgc.tube(iLon, iLat);     // ion temperature
	
	//std::cout << "temp1d(20) " << temp1d[20] << "   for ion species " << iIon << "\n";
	//std::cout << "kb " << cKB << "   kb^2 " << pow(cKB,2) << "\n";
	//std::cout << "T^5/2 " << pow(temp1d[20],2.5) << "\n";
	
	lambda1d = 25 / 8 * pow(cKB,2) * pow(temp1d,2.5) / 
		   species[iIon].mass / species[iIon].nu_ion_ion[iIon]; // thermal conductivity
	
	//std::cout << "ion mass is: " << species[iIon].mass << "\n";
        //std::cout << "Bst value is: " << species[iIon].nu_ion_ion[iIon] << "\n";	
	std::cout << "lambda is: " << lambda1d[20] << "\n";
	std::cout << "lambda bits is: " << 25 / 8 * pow(cKB,2) * pow(temp1d[20],2.5) / 
		   species[iIon].mass / species[iIon].nu_ion_ion[iIon] << "\n\n\n";
	
	den1d    = species[iIon].density_scgc.tube(iLon,iLat);
	
	//std::cout << "density of ion specie is " << den1d[20] << "\n";
	
	front1d  = 2/3/species[iIon].density_scgc.tube(iLon,iLat)/cKB;  // front matter of the term
	
	//std::cout << "front matter is: " << front1d[20] << "\n\n\n\n";
	
	dalt1d   = grid.dalt_lower_scgc.tube(iLon, iLat);               // grid thing for solver
        
	//std::cout << "grid thing is: " << dalt1d[20] << "\n";
        
	conduction1d.zeros();                                           // set temp variable to zero

        conduction1d = solver_conduction(temp1d, lambda1d, front1d, dt, dalt1d);
        
	//std::cout << "solved conduction is: " << conduction1d[20] << "\n";
        // The conduction solver gives Tnew-Told, so divide by dt
        conduction_scgc.tube(iLon, iLat) = conduction1d / dt;
        break; /////////////////////////////////
      } // Lats
    } // Lons

    //std::cout << "lambda " << lambda1d(20) << "  front " << front1d << "  conduction1d " 
//	      << conduction1d(20) << "\n";
    // --------------------------------------------------------------------------
    // Add temperature terms together to advance ion temperature
    // As more temperature terms get coded, they are added to the parenthesis 
    // for inclusion in the advancement of the ion temperature
    // --------------------------------------------------------------------------
    if (input.get_do_calc_bulk_ion_temp()==false) {
      species[iIon].temperature_scgc = species[iIon].temperature_scgc + 
                                         dt * (conduction_scgc);
    }
  } //ions

  if (input.get_do_calc_bulk_ion_temp()==false) {
    // Use the density averaged temperature to fill the bulk temperature 
    tempT.ones();
    for (iIon = 0; iIon < nIons; iIon++)
      tempT = tempT + (species[iIon].temperature_scgc % species[iIon].density_scgc);
    temperature_scgc = tempT / density_scgc;
  }

  if (input.get_do_calc_bulk_ion_temp()==true) {
    // Add temperature terms together to advance bulk ion temperature
    temperature_scgc = temperature_scgc + dt * (conduction_scgc);

    // Use the bulk ion temperature to fill all ion specie temperatures
    for (iIon = 0; iIon < nIons; iIon++)
      species[iIon].temperature_scgc = temperature_scgc;
  }
  //std::cout << "Test ion temp 2 at Alt: 185, Lat: 47.5, Lon: 175 ---> " 
  //	    << temperature_scgc(20,30,20) << "\n";
  report.exit(function);
  return;
}
