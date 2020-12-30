// (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef AETHER_INCLUDE_IONS_H_
#define AETHER_INCLUDE_IONS_H_

#include <armadillo>

#include "inputs.h"
#include "report.h"
#include "grid.h"

using namespace arma;

class Ions {

 public:

  struct species_chars {

    std::string cName;
    float mass;
    int charge;
    
    int DoAdvect;
    
    float *density_s3gc;
    float *par_velocity_v3gc;
    float *perp_velocity_v3gc;

    float *temperature_s3gc;
    
    // Sources and Losses:

    float *ionization_s3gc;
    
    fcube density_scgc;
    fcube par_velocity_vcgc;
    fcube perp_velocity_vcgc;

    fcube temperature_scgc;
    
    // Sources and Losses:

    fcube ionization_scgc;

    fcube sources_scgc;
    fcube losses_scgc;
    
  };

  // bulk quantities (states):
  float *density_s3gc;
  fcube density_scgc;

  float *velocity_v3gc;
  float *exb_v3gc;
  float *ion_temperature_s3gc;
  float *electron_temperature_s3gc;

  fcube ion_temperature_scgc;
  fcube electron_temperature_scgc;
  
  std::vector<species_chars> species;
  
  // ------------------------------
  // Functions:
  
  Ions(Grid grid, Inputs input, Report report);
  species_chars create_species(Grid grid);
  int read_planet_file(Inputs input, Report report);
  void fill_electrons(Grid grid, Report &report);

};
#endif // AETHER_INCLUDE_NEUTRALS_H_
