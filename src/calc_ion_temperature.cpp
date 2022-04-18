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
                                    Times time, Report &report) {
  
  std::string function = "Ions::calc_ion_temperature";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int64_t iIon, iLon, iLat, nSpecies;   
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

  // Loop over all species or assume only bulk calculation
  if (DoCalcBulkIonTemp==true)  nSpecies=1; // First ion specie only, is O+?
  if (DoCalcBulkIonTemp==false) nSpecies=nIons;

  for (iIon = 0; iIon < nSpecies; iIon++) {

    for (iLon = 0; iLon < nLons; iLon++) {
      for (iLat = 0; iLat < nLats; iLat++) {

        // --------------------------------------------------------------------------
        // Calculate heat flux (conduction) in 1D; loop over all lat,lon positions
        // --------------------------------------------------------------------------

        // lambdas=25/8*kb^2*Ts(:,:,is).^(5/2)/ms/(Css(is)*1e-6);
        // Css is ion self collisions, diagional from the Bst table, 
        // Css is species[iIon].nu_ion_ion since the calc_ion_temperature is a function 
        // within the ion class, nu_ion_ion is a vector of length number-of-ion-species

        temp1d   = species[iIon].temperature_scgc.tube(iLon, iLat); // ion temperature
        lambda1d = 25/8*pow(cKB,2)*pow(temp1d,(5/2))/species[iIon].mass/
                   (ions.species[iIon].nu_ion_ion);         // thermal conductivity
        front1d.ones(); //... need to look up what is the front matter of the term being solved;
        dalt1d   = grid.dalt_lower_scgc.tube(iLon, iLat);   // grid thing for solver

        conduction1d.zeros();                               // set temp variable to zero

        conduction1d = solver_conduction(temp1d, lambda1d, front1d, dt, dalt1d);

        // The conduction solver gives Tnew-Told, so divide by dt
        conduction_scgc.tube(iLon, iLat) = conduction1d / dt;

      } // Lats
    } // Lons

    // --------------------------------------------------------------------------
    // Add temperature terms together to advance ion temperature
    // As more temperature terms get coded, they are added to the parenthesis 
    // for inclusion in the advancement of the ion temperature
    // --------------------------------------------------------------------------
    if (DoCalcBulkIonTemp==false) {
      species[iIon].temperature_scgc = species[iIon].temperature_scgc + 
                                         dt * (conduction_scgc);
    }
  } //ions

  if (DoCalcBulkIonTemp==false) {
    // Use the density averaged temperature to fill the bulk temperature 
    // do I need to loop this??? bleh.
    // from ions.cpp : density_scgc = species[nIons].density_scgc; It is electron density.
    temperature_scgc = (species.temperature_scgc * species.density_scgc) / density_scgc;
  }

  if (DoCalcBulkIonTemp==true) {
    // Add temperature terms together to advance bulk ion temperature
    temperature_scgc = temperature_scgc + dt * (conduction_scgc);

    // Use the bulk ion temperature to fill all ion specie temperatures
    for (iIon = 0; iIon < nIons; iIon++)
      species[iIon].temperature_scgc = temperature_scgc;
  }

  report.exit(function);
  return;
}
