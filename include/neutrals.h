// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_NEUTRALS_H_
#define INCLUDE_NEUTRALS_H_

#include <string>
#include <vector>

// The armadillo library is to allow the use of 3d cubes and other
// array types, with array math built in. This eliminates loops!
#include <armadillo>

#include "grid.h"
#include "euv.h"
#include "time.h"
#include "ions.h"
#include "inputs.h"
#include "report.h"

using namespace arma;

class Neutrals {

 public:

  // This struct contains all of the information needed for a single
  // species of neutrals.  We will then have a vector of these species.
  
  struct species_chars {
    std::string cName;
    float mass;
    float vibe;

    int DoAdvect;

    fcube density_scgc;

    std::vector<float> diff0;
    std::vector<float> diff_exp;
    std::vector<float> neutral_ion;

    float thermal_cond;
    float thermal_exp;

    int iEuvAbsId_;
    int nEuvIonSpecies;
    std::vector<int> iEuvIonSpecies_;
    std::vector<int> iEuvIonId_;

    // Some derived quantities:
    fcube chapman_scgc;
    fcube scale_height_scgc;

    // Sources and Losses:

    fcube ionization_scgc;

    fcube sources_scgc;
    fcube losses_scgc;

    // If we want a fixed lower BC:
    float lower_bc_density;
  };

  // bulk quantities (states):

  fcube density_scgc;
  fcube temperature_scgc;

  fcube rho_scgc;
  fcube mean_major_mass_scgc;
  fcube pressure_scgc;
  fcube sound_scgc;

  // For heating/cooling:
  fcube Cv_scgc;
  fcube gamma_scgc;
  fcube kappa_scgc;

  std::vector<species_chars> neutrals;

  float max_chapman = 1.0e26;

  // Source terms:

  fcube conduction_scgc;
  fcube heating_euv_scgc;

  float heating_efficiency;

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
  int pair_euv(Euv euv, Ions ions, Report report);
  void calc_mass_density(Report &report);
  void calc_specific_heat(Report &report);
  void calc_chapman(Grid grid, Report &report);
  void calc_ionization_heating(Euv euv, Ions &ions, Report &report);
  void calc_conduction(Grid grid, Times time, Report &report);
  void add_sources(Times time, Report &report);
  void set_bcs(Report &report);
};

#endif  // INCLUDE_NEUTRALS_H_

