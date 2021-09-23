// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_INDICES_H_
#define INCLUDE_INDICES_H_

#include <vector>
#include <string>

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

  precision_t get_f107(double time);
  precision_t get_f107a(double time);

  int get_f107_index_id();
  int get_f107a_index_id();
  int get_imfbx_index_id();
  int get_imfby_index_id();
  int get_imfbz_index_id();
  int get_swvx_index_id();
  int get_swvy_index_id();
  int get_swvz_index_id();
  int get_swn_index_id();
  int get_swt_index_id();
  int get_ae_index_id();
  int get_au_index_id();
  int get_al_index_id();
  
  // This is the method for setting f107 specifically:
  void set_f107(index_file_output_struct f107_contents);

  // This is the general method for setting indices:
  void set_index(int index_id,
		 std::vector<double> time,
		 std::vector<float> values,
		 precision_t missing);

 private:

  struct index_time_pair {
    std::vector<float> values;
    std::vector<double> times;
    std::string name;
    int64_t nValues;
  };

  std::vector<index_time_pair> all_indices_arrays;

  const int iF107_ = 0;
  const int iF107A_ = 1;
  const int iIMFBX_ = 2;
  const int iIMFBY_ = 3;
  const int iIMFBZ_ = 4;
  const int iSWVX_ = 5;
  const int iSWVY_ = 6;
  const int iSWVZ_ = 7;
  const int iSWN_ = 8;
  const int iSWT_ = 9;
  const int iAE_ = 10;
  const int iAL_ = 11;
  const int iAU_ = 12;
  int nIndices = 13;
  
  precision_t get_index(double time, int index);
};


#endif // INCLUDE_INDICES_H_
