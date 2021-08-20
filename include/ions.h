// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_IONS_H_
#define INCLUDE_IONS_H_

#include <string>
#include <vector>

// The armadillo library is to allow the use of 3d cubes and other
// array types, with array math built in. This eliminates loops!
#include <armadillo>

#include "inputs.h"
#include "report.h"
#include "grid.h"

using namespace arma;

class Ions {

 public:

  // This struct contains all of the information needed for a single
  // species of ion.  We will then have a vector of these species.

  struct species_chars {
    std::string cName;
    float mass;
    int charge;

    int DoAdvect;

    std::vector<float> nu_ion_neutral_coef;
    std::vector<bool> nu_is_resonant;
    std::vector<float> nu_ion_ion;
    std::vector<float> nu_ion_electron;

    std::vector<float> nu_in_res_temp_min;
    std::vector<float> nu_in_res_coef1;
    std::vector<float> nu_in_res_coef2;
    std::vector<float> nu_in_res_tn_frac;
    std::vector<float> nu_in_res_ti_frac;

    std::vector<fcube> nu_ion_neutral_vcgc;
    
    // Sources and Losses:

    fcube density_scgc;
    std::vector<fcube> par_velocity_vcgc;
    std::vector<fcube> perp_velocity_vcgc;

    fcube temperature_scgc;

    // Sources and Losses:

    fcube ionization_scgc;

    fcube sources_scgc;
    fcube losses_scgc;
  };

  // bulk quantities (states):
  fcube density_scgc;
  std::vector<fcube> velocity_vcgc;
  fcube ion_temperature_scgc;
  fcube electron_temperature_scgc;

  // This is the vector that will contain all of the different species:
  std::vector<species_chars> species;

  // Electrodynamics:
  fcube potential_scgc;
  std::vector<fcube> efield_vcgc;
  std::vector<fcube> exb_vcgc;
  fmat eflux;
  fmat avee;

  // ------------------------------
  // Functions:

  Ions(Grid grid, Inputs input, Report report);
  species_chars create_species(Grid grid);
  int read_planet_file(Inputs input, Report report);
  void fill_electrons(Report &report);
  int get_species_id(std::string name, Report &report);
  void calc_efield(Grid grid, Report &report);
  void calc_exb_drift(Grid grid, Report &report);
  void calc_ion_drift(Neutrals neutrals, Grid grid, float dt, Report &report);
  void calc_ion_neutral_coll_freq(Neutrals neutrals, Report &report);
  std::vector<fcube> calc_ion_electron_pressure_gradient(int64_t iIon,
							 Grid grid,
							 Report &report);
  void calc_ion_temperature(Neutrals neutrals, Grid grid, Report &report);
};
#endif  // INCLUDE_IONS_H_
