// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_TIME_CONVERSION_H_
#define INCLUDE_TIME_CONVERSION_H_

/**************************************************************
 * These are functions that convert time between different systems.
 * - real time is seconds past January 1, 1965 00 UT
 * - integer time is a vector of:
 *        (year, month, day, hour, minute, second, millisecond)
 * - jday is julian day, which is the day of year (jan 1 = 1)
 **************************************************************/

#include <vector>

/**************************************************************
  \brief convert time from integer array to day of year (jday)
  \param itime integer vector (y,m,d,h,m,s,ms)
**/
double time_int_to_jday(std::vector<int> itime);

/**************************************************************
  \brief converts year, month, day to day of year
  \param year the year
  \param month the month (1 = jan)
  \param day day of the month (1 = 1)
**/
int day_of_year(int year, int month, int day);

/**************************************************************
  \brief converts integer vector time to real time
  \param itime integer vector (y,m,d,h,m,s,ms)
**/
double time_int_to_real(std::vector<int> itime);

/**************************************************************
  \brief converts real time to an integer vector (y,m,d,h,m,s,ms)
  \param timereal seconds past January 1, 1965 00 UT
**/
std::vector<int> time_real_to_int(double timereal);

/**************************************************************
  \brief runs a test on the time conversion routines.
**/
int test_time_routines();

/**************************************************************
  \brief displays the integer time array
  \param itime integer vector (y,m,d,h,m,s,ms)
**/
void display_itime(std::vector<int> itime);

#endif  // INCLUDE_TIME_CONVERSION_H_
