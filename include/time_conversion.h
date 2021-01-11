// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_TIME_CONVERSION_H_
#define INCLUDE_TIME_CONVERSION_H_

#include <vector>

double time_int_to_jday(std::vector<int> itime);
int day_of_year(int year, int month, int day);
double time_int_to_real(std::vector<int> itime);
void time_real_to_int(double timereal, std::vector<int> &itime);
int test_time_routines();
void display_itime(std::vector<int> itime);

#endif  // INCLUDE_TIME_CONVERSION_H_

