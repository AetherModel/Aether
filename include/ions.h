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

    std::vector<precision_t> nu_ion_neutral_coef;
    std::vector<bool> nu_is_resonant;
    std::vector<precision_t> nu_ion_ion;
    std::vector<precision_t> nu_ion_electron;

    std::vector<precision_t> nu_in_res_temp_min;
    std::vector<precision_t> nu_in_res_coef1;
    std::vector<precision_t> nu_in_res_coef2;
    std::vector<precision_t> nu_in_res_tn_frac;
    std::vector<precision_t> nu_in_res_ti_frac;

    std::vector<arma_cube> nu_ion_neutral_vcgc;
    
    // Sources and Losses:

    arma_cube density_scgc;
    std::vector<arma_cube> par_velocity_vcgc;
    std::vector<arma_cube> perp_velocity_vcgc;

    arma_cube temperature_scgc;

    // Sources and Losses:

    arma_cube ionization_scgc;

    arma_cube sources_scgc;
    arma_cube losses_scgc;
  };

  // bulk quantities (states):
  arma_cube density_scgc;
  std::vector<arma_cube> velocity_vcgc;
  arma_cube ion_temperature_scgc;
  arma_cube electron_temperature_scgc;

  // This is the vector that will contain all of the different species:
  std::vector<species_chars> species;

  // Electrodynamics:
  arma_cube potential_scgc;
  std::vector<arma_cube> efield_vcgc;
  std::vector<arma_cube> exb_vcgc;
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
  void calc_ion_drift(Neutrals neutrals,
		      Grid grid,
		      precision_t dt,
		      Report &report);
  void calc_ion_neutral_coll_freq(Neutrals neutrals, Report &report);
  std::vector<arma_cube> calc_ion_electron_pressure_gradient(int64_t iIon,
							     Grid grid,
							     Report &report);
  void calc_ion_temperature(Neutrals neutrals, Grid grid, Report &report);
};
#endif  // INCLUDE_IONS_H_
