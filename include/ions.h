// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_IONS_H_
#define INCLUDE_IONS_H_

#include <string>
#include <vector>

class Ions {

 public:

  // This struct contains all of the information needed for a single
  // species of ion.  We will then have a vector of these species.

  struct species_chars {
    std::string cName;
    precision_t mass;
    int charge;

    int DoAdvect;

    // Sources and Losses:

    arma_cube density_scgc;
    arma_cube par_velocity_vcgc;
    arma_cube perp_velocity_vcgc;

    arma_cube temperature_scgc;

    // Sources and Losses:

    arma_cube ionization_scgc;

    arma_cube sources_scgc;
    arma_cube losses_scgc;
  };

  // bulk quantities (states):
  arma_cube density_scgc;

  arma_cube ion_temperature_scgc;
  arma_cube electron_temperature_scgc;

  // This is the vector that will contain all of the different species:
  std::vector<species_chars> species;

  // ------------------------------
  // Functions:

  Ions(Grid grid, Inputs input, Report report);
  species_chars create_species(Grid grid);
  int read_planet_file(Inputs input, Report report);
  void fill_electrons(Report &report);

};
#endif  // INCLUDE_IONS_H_
