// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_INDICES_H_
#define INCLUDE_INDICES_H_

#include <vector>
#include <string>

#include "inputs.h"

struct index_file_output_struct {
  std::vector<double> times;
  int64_t nTimes;
  std::vector<std::string> var_names;
  int nVars;
  std::vector<std::vector<float>> values;
  std::vector<float> missing_values;
  std::vector<int> index_id;
};

void print_index_file_output_struct(index_file_output_struct contents);


class Indices {

 public:

  // Public Functions:

  Indices(Inputs args);

  float get_f107(double time);
  float get_f107a(double time);

  int get_f107_index_id();
  int get_f107a_index_id();
  int get_imf_bx_index_id();
  int get_imf_by_index_id();
  int get_imf_bz_index_id();
  int get_sw_vx_index_id();
  int get_sw_vy_index_id();
  int get_sw_vz_index_id();
  int get_sw_n_index_id();
  int get_sw_t_index_id();
  int get_ae_index_id();
  int get_au_index_id();
  int get_al_index_id();
  
  // This is the method for setting f107 specifically:
  void set_f107(index_file_output_struct f107_contents);

  // This is the general method for setting indices:
  void set_index(int index_id,
		 std::vector<double> time,
		 std::vector<float> values,
		 float missing);

 private:

  struct index_time_pair {
    std::vector<float> values;
    std::vector<double> times;
    std::string name;
    int64_t nValues;
  };

  std::vector<index_time_pair> all_indices_arrays;

  int f107_ = 0;
  int f107a_ = 1;
  int imf_bx_ = 2;
  int imf_by_ = 3;
  int imf_bz_ = 4;
  int sw_vx_ = 5;
  int sw_vy_ = 6;
  int sw_vz_ = 7;
  int sw_n_ = 8;
  int sw_t_ = 9;
  int ae_ = 10;
  int al_ = 11;
  int au_ = 12;
  int nIndices = 13;
  
  float get_index(double time, int index);
};


#endif // INCLUDE_INDICES_H_
