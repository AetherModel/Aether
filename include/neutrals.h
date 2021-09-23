// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_NEUTRALS_H_
#define INCLUDE_NEUTRALS_H_

#include <string>
#include <vector>

class Neutrals {

 public:

  // This struct contains all of the information needed for a single
  // species of neutrals.  We will then have a vector of these species.
  
  struct species_chars {
    std::string cName;
    precision_t mass;
    precision_t vibe;

    int DoAdvect;

    arma_cube density_scgc;

    std::vector<float> diff0;
    std::vector<float> diff_exp;
    std::vector<float> neutral_ion;

    precision_t thermal_cond;
    precision_t thermal_exp;

    int iEuvAbsId_;
    int nEuvIonSpecies;
    std::vector<int> iEuvIonSpecies_;
    std::vector<int> iEuvIonId_;

    // Some derived quantities:
    arma_cube chapman_scgc;
    arma_cube scale_height_scgc;

    // Sources and Losses:

    arma_cube ionization_scgc;

    arma_cube sources_scgc;
    arma_cube losses_scgc;

    // If we want a fixed lower BC:
    precision_t lower_bc_density;
  };

  // bulk quantities (states):

  arma_cube density_scgc;
  arma_cube temperature_scgc;

  arma_cube rho_scgc;
  arma_cube mean_major_mass_scgc;
  arma_cube pressure_scgc;
  arma_cube sound_scgc;

  // For heating/cooling:
  arma_cube Cv_scgc;
  arma_cube gamma_scgc;
  arma_cube kappa_scgc;

  std::vector<species_chars> species;

  precision_t max_chapman = 1.0e26;

  // Source terms:

  arma_cube conduction_scgc;
  arma_cube heating_euv_scgc;

  precision_t heating_efficiency;

  // This is an initial temperature profile, read in through the
  // planet.in file:
  float *initial_temperatures, *initial_altitudes;
  int nInitial_temps = 0;

  // names and units
  std::string density_unit = "(/m3)";
  std::string density_name = "Neutral Bulk Density";

  std::string velocity_unit = "(m/s)";
  std::vector<std::string> velocity_name;

  std::string temperature_unit = "(K)";
  std::string temperature_name = "Temperature";

  // ------------------------------
  // Functions:

  Neutrals(Grid grid, Inputs input, Report report);
  species_chars create_species(Grid grid);
  int read_planet_file(Inputs input, Report report);
  int initial_conditions(Grid grid, Inputs input, Report report);
  void fill_with_hydrostatic(Grid grid, Report report);
  void calc_mass_density(Report &report);
  void calc_specific_heat(Report &report);
  void calc_chapman(Grid grid, Report &report);
  void calc_conduction(Grid grid, Times time, Report &report);
  void add_sources(Times time, Report &report);
  void set_bcs(Report &report);
};

#endif  // INCLUDE_NEUTRALS_H_

