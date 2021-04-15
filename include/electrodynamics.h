// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_ELECTRODYNAMICS_H_
#define INCLUDE_ELECTRODYNAMICS_H_

/**************************************************************
 * 
 *
 *
 **************************************************************/

#include <vector>

#include "aether.h"

class Electrodynamics {

 public:

  Electrodynamics(Inputs input, Report report);
  void set_time(double time);
  void set_grid(fmat lats, fmat mlts, Report report);
  void set_imf_bx(float value);
  void set_imf_by(float value);
  void set_imf_bz(float value);
  void set_sw_v(float value);
  void set_sw_n(float value);
  void set_hp(float value);
  void set_au(float value);
  void set_al(float value);
  void set_ae(float value);
  void set_kp(float value);
  
  /**************************************************************
   * 
   **/
  fmat get_potential(Report report);
  fmat get_eflux(Report report);
  fmat get_avee(Report report);
  fmat get_ion_eflux(Report report);
  fmat get_ion_avee(Report report);
  
 private:

  int iTimeInterpolationMethod;

  std::string input_file;

  /// Variables needed by user:
  double time_needed;
  fmat lats_needed;
  fmat mlts_needed;
  float imf_bx_needed;
  float imf_bx_needed;
  float imf_bx_needed;
  float sw_v_needed;
  float sw_n_needed;
  float hp_needed;
  float au_needed;
  float al_needed;
  float ae_needed;
  float kp_needed;

  int iUseGridBasedModel;
  std::string efield_model_to_use;
  std::string auroral_model_to_use;
  
  /// interpolation indices, for each one, the integer portion is
  /// the index, while the decimal part is the percentage of the way
  /// between the index and index + 1
  float time_index;
  fmat lats_index;
  fmat mlts_index;

  /// This is which index of the input_electrodynamics vector should be used:
  imat map_index;
    
  struct input_electrodynamics_struct {

    int nLats;
    int nMlts;

    fvec mlats;
    fvec mlts;

    std::vector<double> times;
    
    std::vector<fmat> energy_flux;
    std::vector<fmat> average_energy;
    std::vector<fmat> ion_energy_flux;
    std::vector<fmat> ion_average_energy;
    std::vector<fmat> potential;

  }

  std::vector<input_electrodynamics_struct> input_electrodynamics;

  
  
};

#endif // INCLUDE_ELECTRODYNAMICS_H_
