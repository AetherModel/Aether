// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_INDICES_H_
#define INCLUDE_INDICES_H_

#include <vector>
#include <string>

class Indices {

 public:

  // Public Functions:

  Indices(Inputs args);

  float get_f107(double time);
  float get_f107a(double time);

 private:

  struct ind_time_pair {
    float index;
    double time;
  };

  std::vector<ind_time_pair> f107;
  std::vector<ind_time_pair> f107a;

  // This is the method for setting f107 specifically:
  int set_f107(std::vector<double> time, std::vector<float> f107array);

  // This is the general method for setting indices:
  int set_index(std::vector<double> time,
                std::vector<float> indexarray,
                std::vector<ind_time_pair> &index);

  float get_index(double time, std::vector<ind_time_pair> index);
};

#endif // INCLUDE_INDICES_H_
