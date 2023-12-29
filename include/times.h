// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_TIMES_H_
#define INCLUDE_TIMES_H_

/**************************************************************
 * 
 * times.h:
 *
 *   Functions that are assocuated with keeping track of time in
 * Aether.  Once central concept in Aether, is that time is stored as
 * a double and is the number of seconds from January 1, 1965. This is
 * arbitrary, and is derived from code developed by UCLA in the 1990s.
 *
 **************************************************************/

#include <vector>
#include <string>

class Times {

public:

  /**************************************************************
     \brief Initialize the Times class
  **/
  Times();

  /**************************************************************
     \brief Increments the current time and derived variables

     Increments the current time, simulation time, and iStep.
   Then derives a bunch of variables from there, in order to allow
   users to use more than the time in seconds. This function should
   only be called once per iteration at the end of advance.
   **/
  void increment_time();

  /**************************************************************
     \brief Increments the intermediate/coupling time by input dt

     This increments the intermediate time (or coupling time) by
   the input time.  The intermediate time is set as the time in which
   Aether should couple with any external models (in main). When
   Aether is run in stand-alone mode, this doesn't really do much.

   \param dt amount of time to move the intermediate time forward.
   **/
  void increment_intermediate(double dt);

  /**************************************************************
     \brief Displays the current time, elapsed time, and time to complete

     This displays the current time, the elapsed wall time, estimates
   a time to completion (remaining_walltime) and displays that also.
   **/
  void display();

  /**************************************************************
     \brief Sets the start, restart, and current times.

     This sets the start time, restart time, and current time to
   the input time, and initializes iStep and dt, then calls 
   increment_time, which derives a bunch of other variables.

   \param itime year, month, day, hour, minute, second, millisecond vector
   **/
  void set_times(std::vector<int> itime);

  /**************************************************************
     \brief Sets end time

     This simply ses the end time of the simulation from the input.

   \param itime year, month, day, hour, minute, second, millisecond vector
   **/
  void set_end_time(std::vector<int> itime);

  /**************************************************************
     \brief Returns the current time in seconds since ref date
   **/
  double get_current();

  /**************************************************************
     \brief Returns the end time in seconds since ref date
   **/
  double get_end();

  /**************************************************************
     \brief Returns the current time as a string (year, month...)

     Returns the current time as a string of format YYYYMMDD_HHMMSS

     \param useSeconds if false, replace seconds with 0
   **/
  std::string get_YMD_HMS(bool useSeconds);

  /**************************************************************
     \brief Returns the intermediate time in seconds since ref date
   **/
  double get_intermediate();

  /**************************************************************
     \brief Returns dt

     dt is the delta-time between iterations in Aether.
   **/
  precision_t get_dt();

  /**************************************************************
     \brief Returns current time in specific unit for planetary calc

     Return the current time as an orbit time, to allow for calculation
   of a bunch of planetary characteristics. It's a JPL thing.
   **/
  precision_t get_orbittime();

  /**************************************************************
     \brief Returns current time in Julian Days
   **/
  double get_julian_day();

  /**************************************************************
     \brief Check to see if time just passed through dt gate.

     This function checks to see if the simulation has passed through
   a time gate. By this, I mean that the user sets a dt_check in which
   to do something. If the simulation passes through that dt_check,
   then it will return a 1, else it will return a 0.  It also returns
   1 if the current time is the start time or restart time.

   \param dt_check Repetative delta-time to do a task

   **/
  int check_time_gate(precision_t dt_check);

  /**************************************************************
     \brief Calculates the delta time in the code.
     \param dtNeutral dt from the neutrals
     \param dtIon dt from the ions
     This calculates the delta-time in Aether for the current state.
   At this moment, this is a stand-in code. In reality, we need to pass
   in the grid, neutrals, ions, and inputs for calculation.
   **/
  void calc_dt(precision_t dtNeutral, precision_t dtIon);

  /**************************************************************
     \brief Get the current time as an array
   **/
  std::vector<int> get_iCurrent();
  
  /**************************************************************
     \brief Get the current simulation time (sec since start)
   **/
  double get_simulation_time();
  
  /**********************************************************************
     \brief Read / Write restart files for time
     \param dir directory to write restart files
     \param DoRead read the restart files if true, write if false
   **/
  bool restart_file(std::string dir, bool DoRead);  

private:
  
  // -------------------------------------------------------------
  // These variables are for keeping track of the time. All in seconds
  // since reference time (except where noted).

  /// Start time of the simulation
  double start;

  /// Restart time of the simulation (if not restart = start time)
  double restart;

  /// End time of the simulation
  double end;

  /// Current time of the simulation
  double current;

  /// Time to stop and couple with other codes, if applicable
  double intermediate;

  /// Seconds since the start of the simulation
  double simulation;

  /// number of iterations into the simulation
  int64_t iStep;

  /// delta-time between current time-steps
  precision_t dt;

  // -------------------------------------------------------------
  // Derived variables from the current time:

  /// year, month, day, hour, minute, second, millisecond vector
  std::vector<int> iCurrent;

  /// Universal time in hours
  precision_t ut;
  
  /// in weird JPL units
  precision_t orbittime;

  /// current time as different integer units:
  int year;
  int month;
  int day;
  int hour;
  int minute;
  int second;
  int milli;

  /// This is day of year (and NOT real Julian Day!)
  int jDay;
  
  /// This is Julian day
  double julian_day;

  /// represented as YYMMDD
  std::string sYMD;
  
  /// represented as HHMMSS
  std::string sHMS;

  /// represented as YYYYMMDD_HHMMSS
  std::string sYMD_HMS;

  /// represented as YYYYMMDD_HHMM00
  std::string sYMD_HM0;  

  // -------------------------------------------------------------
  // Keeping track of walltime for the run:

  /// This is the system time at the start of the simulation
  time_t sys_time_start;

  /// This is the current system time of the simulation
  time_t sys_time_current;

  /// This is the difference, or the current walltime in seconds
  precision_t walltime;
};

#endif  // INCLUDE_TIMES_H_
