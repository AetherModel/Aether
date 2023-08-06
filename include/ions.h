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

  int64_t nSpecies = 8;
  
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
    arma_cube conduction_scgc;

    arma_cube ionization_scgc;

    arma_cube sources_scgc;
    arma_cube losses_scgc;
  };

  // bulk quantities (states):
  arma_cube density_scgc;
  std::vector<arma_cube> velocity_vcgc;
  arma_cube temperature_scgc;
  arma_cube conduction_scgc;
  arma_cube electron_temperature_scgc;

  // This is the vector that will contain all of the different species:
  std::vector<species_chars> species;

  // Electrodynamics:
  arma_cube potential_scgc;
  std::vector<arma_cube> efield_vcgc;
  std::vector<arma_cube> exb_vcgc;
  arma_mat eflux;
  arma_mat avee;

  // names and units
  const std::string density_name = "Neutral Bulk Density";
  const std::string density_unit = "/m3";

  std::vector<std::string> velocity_name;
  std::vector<std::string> par_velocity_name;
  std::vector<std::string> perp_velocity_name;
  const std::string velocity_unit = "m/s";

  const std::string temperature_name = "Temperature";
  const std::string temperature_unit = "K";

  const std::string potential_name = "Potential";
  const std::string potential_unit = "Volts";
  
  // ------------------------------
  // Functions:

  Ions(Grid grid, Planets planet);
  species_chars create_species(Grid grid);
  int read_planet_file(Planets planet);
  void init_ion_temperature(Neutrals neutrals, Grid grid);
  void fill_electrons();
  int get_species_id(std::string name);
  void calc_efield(Grid grid);
  void calc_exb_drift(Grid grid);
  void calc_ion_drift(Neutrals neutrals,
		      Grid grid,
		      precision_t dt);
  std::vector<arma_cube> calc_ion_electron_pressure_gradient(int64_t iIon,
							     Grid grid);
  void calc_ion_temperature(Neutrals neutrals, Grid grid, Times time);
  void calc_electron_temperature(Neutrals neutrals, Grid grid);
  
  bool restart_file(std::string dir, bool DoRead);
};
#endif  // INCLUDE_IONS_H_
