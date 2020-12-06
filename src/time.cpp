// (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <math.h>
#include <vector>
#include <string>
#include <iostream>
#include <time.h>

#include "../include/times.h"
#include "../include/sizes.h"
#include "../include/time_conversion.h"

// -----------------------------------------------------------------------------
// Initialize the time variables
// -----------------------------------------------------------------------------

Times::Times() {
  
  iCurrent = {0, 0, 0, 0, 0, 0, 0};
  iStep = -1;

  year = 0;
  month = 0;
  day = 0;
  jDay = 0;
  hour = 0;
  minute = 0;
  second = 0;
  milli = 0;

  time(&sys_time_start);
}

// -----------------------------------------------------------------------------
// 
// -----------------------------------------------------------------------------

void Times::set_times(std::vector<int> itime) {

  start = time_int_to_real(itime);
  current = start;
  iStep = -1;
  dt = 0;
  // This will initiate more variables:
  increment_time();
  
}

// -----------------------------------------------------------------------------
// This function checks to see if the simulation has passed through a
// time gate By this, I mean that the user sets a dt in which to do
// something. If the simulation passes through that dt, then it will
// return a 1, else it will return a 0.  It also returns 1 if the
// current time is the start time.
// -----------------------------------------------------------------------------

int Times::check_time_gate(float dt_check) {
  int DoThing = 0;
  if (current == start) DoThing = 1;
  if ( floor((simulation - dt) / dt_check) <
       floor(simulation / dt_check)) DoThing = 1;
  return DoThing;
}

// -----------------------------------------------------------------------------
// Need to actually calculate dt here....
// -----------------------------------------------------------------------------

void Times::calc_dt() {
  dt = 5.0;
  return;
}

// -----------------------------------------------------------------------------
// 
// -----------------------------------------------------------------------------

double Times::get_current() {
  return current;
}

// -----------------------------------------------------------------------------
// 
// -----------------------------------------------------------------------------

double Times::get_end() {
  return end;
}

// -----------------------------------------------------------------------------
// 
// -----------------------------------------------------------------------------

std::string Times::get_YMD_HMS() {
  return sYMD_HMS;
}

// -----------------------------------------------------------------------------
// 
// -----------------------------------------------------------------------------

double Times::get_intermediate() {
  return intermediate;
}

// -----------------------------------------------------------------------------
// 
// -----------------------------------------------------------------------------

float Times::get_dt() {
  return dt;
}

// -----------------------------------------------------------------------------
// 
// -----------------------------------------------------------------------------

float Times::get_orbittime() {
  return orbittime;
}

  
// -----------------------------------------------------------------------------
// 
// -----------------------------------------------------------------------------

double Times::get_julian_day() {
  return julian_day;
}

  
// -----------------------------------------------------------------------------
// 
// -----------------------------------------------------------------------------

void Times::set_end_time(std::vector<int> itime) {
  end = time_int_to_real(itime);
}

  
// -----------------------------------------------------------------------------
// Increment the time by dt and increment the iteration number (iStep)
// -----------------------------------------------------------------------------

void Times::increment_time() {

  // Increment simulation time:
  simulation += dt;

  // Increment iStep (iteration number):
  iStep++;

  // Increment current time:
  current += dt;

  // Convert current time to array:
  time_real_to_int(current, iCurrent);

  // Set named variables:
  year = iCurrent[0];
  month = iCurrent[1];
  day = iCurrent[2];
  hour = iCurrent[3];
  minute = iCurrent[4];
  second = iCurrent[5];
  milli = iCurrent[6];

  char tmp[100];
  sprintf(tmp, "%04d%02d%02d_%02d%02d%02d",
          year, month, day, hour, minute, second);
  sYMD_HMS = std::string(tmp);
  sprintf(tmp, "%04d%02d%02d", year, month, day);
  sYMD = std::string(tmp);
  sprintf(tmp, "%02d%02d%02d", hour, minute, second);
  sHMS = std::string(tmp);
  
  // Calculate Julian Day (day of year):
  jDay = day_of_year(year, month, day);

  // Calculate True Julian Day:
  julian_day = time_int_to_jday(iCurrent);

  // Calculate UT (in hours):
  ut = float(iCurrent[3])    // hours
    + float(iCurrent[4])/60.0    // minutes
    + (float(iCurrent[5]) + float(iCurrent[6])/1000)/3600.0;

  // Calculate orbital parameters based on E Standish,
  // Solar System Dynamics, JPL,.
  // No constant orbital speed assumption

  float day_number = 367.0*float(iCurrent[0])
    - 7.0*(float(iCurrent[0])+(float(iCurrent[1])+9.0)/12.0)/4.0
    + 275.0 * float(iCurrent[1])/9.0
    + float(iCurrent[2])
    - 730531.5
    + ut/24.0;

  orbittime = day_number/36525.0;

  return;
}

// -----------------------------------------------------------------------------
// Display Current Time
// -----------------------------------------------------------------------------

void Times::display() {

  time(&sys_time_current);
  walltime = double(sys_time_current) - double(sys_time_start);
  //cout << "Elapsed walltime : " << walltime << "\n";

  //cout << "current time (double) : " << current << "\n";
  std::cout << "Current Time : ";
  for (int i=0;i<7;i++) std::cout << iCurrent[i] << " ";
  std::cout << "\n";
  // cout << "  (as julian day): " << julian_day << "\n";
  // cout << "  (as julian day since j2000): " << julian_day-j2000 << "\n";

  return;

}

// -----------------------------------------------------------------------------
// We have two main stop times in Aether:
//   - intermediate: to do something every once in a while, like couple
//   - end: to end the whole simulation
// -----------------------------------------------------------------------------

void Times::increment_intermediate(double dt) {

  intermediate = current + dt;
  if (intermediate > end) intermediate = end;

  return;

}
