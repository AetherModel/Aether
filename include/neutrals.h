// (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef AETHER_INCLUDE_NEUTRALS_H_
#define AETHER_INCLUDE_NEUTRALS_H_

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

  struct species_chars {

    std::string cName;
    float mass;
    float vibe;

    int DoAdvect;
    
    fcube density_scgc;

    float *density_s3gc;
    float *velocity_v3gc;

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
    float *chapman_s3gc;
    fcube chapman_scgc;
    fcube scale_height_scgc;
    
    // Sources and Losses:

    float *ionization_s3gc;

    // If we want a flat lower BC:
    float lower_bc_density;
    
  };
  
  // bulk quantities (states):
  float *density_s3gc;
  float *velocity_v3gc;
  float *temperature_s3gc;

  fcube density_scgc;
  fcube temperature_scgc;
  
  float *rho_s3gc;
  float *mean_major_mass_s3gc;
  float *pressure_s3gc;
  float *sound_s3gc;
  
  // For heating/cooling:
  float *Cv_s3gc;
  float *gamma_s3gc;
  float *kappa_s3gc;

  std::vector<species_chars> neutrals;

  float max_chapman = 1.0e26;

  // Source terms:

  float *heating_euv_s3gc;
  float *conduction_s3gc;

  float heating_efficiency;
  
  // This is an initial temperature profile, read in through the
  // planet.in file:
  float *initial_temperatures, *initial_altitudes;
  int nInitial_temps=0;

  // names and units
  std::string density_unit="(/m3)";
  std::string density_name="Neutral Bulk Density";

  std::string velocity_unit="(m/s)";
  std::vector<std::string> velocity_name;

  std::string temperature_unit="(K)";
  std::string temperature_name="Temperature";
  
  // ------------------------------
  // Functions:
  
  Neutrals(Grid grid, Inputs input, Report report);
  species_chars create_species(Grid grid);
  int read_planet_file(Inputs input, Report report);
  int initial_conditions(Grid grid, Inputs input, Report report);
  float calc_scale_height(int iSpecies,
			  long index,
			  Grid grid);
  int pair_euv(Euv euv, Ions ions, Report report);
  void calc_mass_density(Report &report);
  void calc_specific_heat(Report &report);
  void calc_chapman(Grid grid, Report &report);
  void calc_ionization_heating(Euv euv, Ions &ions, Report &report);
  void calc_conduction(Grid grid, Times time, Report &report);
  void add_sources(Times time, Report &report);
  
};
  


#endif // AETHER_INCLUDE_NEUTRALS_H_

