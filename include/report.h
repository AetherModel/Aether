// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_REPORT_H_
#define INCLUDE_REPORT_H_

/**************************************************************
 * \class Report
 *
 * \brief A system for reporting within the program

   This class allows users to manage reporting information to the
   user, including debugging information (through a verbose level) and
   timing information for parts of the code. There are functions for
   capturing when the code goes into and out of a function to allow
   for code tracing and timing.
 *
 * \author Aaron Ridley
 *
 * \date 2021/04/16
 **************************************************************/

#include <string>
#include <vector>

#include "aether.h"

class Report {

// -----------------------------------------------------------------------
// Public functions and variables
// -----------------------------------------------------------------------

public:

  // Functions:

  /**************************************************************
   \brief Initialize all of the internal variables.
   **/
  Report();

  /**************************************************************
   \brief This sets the global verbose level of the code

   The verbose level sets the amount of information presented to the
   user of the code. Some thoughts on levels:
   0 - essentially only report timing information
   1 - report only broadest level of detail
   2 - report going in and out of functions, and maybe things that are
       done once per iteration
   3-4 - Once per iteration types of outputs. More details at higher levels
   5+ - get into grid iterations. Produces a HUGE amount of info.

   This verbose level should be set only on one processor. All other
   processors should be set to -1.

   \param input global verbose level to set the code.
   **/
  void set_verbose(int input);

  /**************************************************************
   \brief How deep to go in the timing report at the end of the run
   \param input the timing depth to be set for the code
   **/
  void set_timing_depth(int input);

  /**************************************************************
   \brief limit the timing report to functions that take long time to run
   \param input minimum percent of total time to report
   **/
  void set_timing_percent(float input);

  /**************************************************************
   \brief Print message if iLevel <= verbose level of code
   \param iLevel test against verbose level of the code
   \param output_string string to output if iLevel <= verbose level
   **/
  void print(int iLevel, std::string output_string);

  /**************************************************************
   \brief Returns 1 if iLevel <= verbose level of code, 0 otherwise
   \param iLevel test against verbose level of the code
   **/
  int test_verbose(int iLevel);

  /**************************************************************
   \brief Returns the global verbosity of the code
   **/
  int get_verbose();

  /**************************************************************
   \brief Starts timer and reports when entering a function, if applicable

   This is typically placed at the start of a function. This function
   does the following:
   - Figures out the depth of the function (which is compared to verbose level)
   - Reports that it has entered function if the verbose level is high enough
   - Starts a timer
   - Records the name of the function and number of the function for exit.

   \param function_name The name of the function that is being entered
   \param iFunction first time entered, this should be -1 or something,
                    then this is altered to be the number of the function.
                    when the function is called again, this number saves
                    a string compare
   **/
  void enter(std::string function_name, int &iFunction);

  /**************************************************************
   \brief End the timer and reports exiting the function, see enter above
   \param function_name The name of the function that is being exited
   **/
  void exit(std::string function_name);

  /**************************************************************
   \brief Loop through all reported functions and report their run times
   **/
  void times();

// -----------------------------------------------------------------------
// Private functions and variables
// -----------------------------------------------------------------------
private:

  /// global verbose level of the code
  int iVerbose;
  /// the depth of the reporting for the timing at the end of the simulation
  int iTimingDepth;
  /// Only report times above the given percentage of the total run time:
  float TimingPercent;

  /// This is the information needed to be stored for each "entry" (when the
  /// enter function is called - typically a function, but could just be a
  /// section of code):
  struct item_struct {
    /// this is the name of the iten (or function)
    std::string entry;
    /// number of times this entry is called
    int nTimes;
    /// The total accumulated wall time time that this entry takes to run
    precision_t timing_total;
    /// This is the start-gate for the timer
    unsigned long long timing_start;
    /// This is the level of the function that is then compared to verbose
    int iLevel;
    /// This is a string that holds all of the function names above this one
    int iStringPosBefore;
    /// This is the function that was called just before this one, so that
    /// if can be "popped" of the queue:
    int iLastEntry;
  };

  /// A vector of entries to keep track of during the model run
  std::vector<item_struct> entries;
  /// number of entries that the code is keeping track of
  int nEntries;
  /// which entry is the current one (string)
  std::string current_entry;
  /// which entry is the current one (integer)
  int iCurrentFunction = -1;

  /// what string to use to tie together the function names
  std::string divider;
  /// the length of the divider string, just to keep track of it
  int divider_length;

  /// The current depth level of the function calls
  int iLevel;

  /// Report when leaving a function
  bool DoReportOnExit = true;
};

#endif  // INCLUDE_REPORT_H_
