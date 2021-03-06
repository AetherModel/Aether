// (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef AETHER_INCLUDE_TIMES_H_
#define AETHER_INCLUDE_TIMES_H_

#include <vector>
#include <string>

class Times {

public:

  Times();
  void increment_time();
  void increment_intermediate(double dt);
  void display();
  void set_times(std::vector<int> itime);
  void set_end_time(std::vector<int> itime);
  double get_current();
  double get_end();
  std::string get_YMD_HMS();
  double get_intermediate();
  float get_dt();
  float get_orbittime();
  double get_julian_day();

  int check_time_gate(float dt_check);

  void calc_dt();
  
private:
  
  // These variables are for keeping track of the time:
  double start;
  double end;
  double current;
  double intermediate;
  double simulation;
  long iStep;

  float dt;

  // Derived variables from the current time:
  std::vector<int> iCurrent;
  float ut;
  float orbittime;
  int year;
  int month;
  int day;
  int hour;
  int minute;
  int second;
  int milli;
  int jDay;
  double julian_day;
  std::string sYMD;
  std::string sHMS;
  std::string sYMD_HMS;

  // Keeping track of walltime for the run:
  time_t sys_time_start;
  time_t sys_time_current;
  float walltime;

};

#endif // AETHER_INCLUDE_TIMES_H_
