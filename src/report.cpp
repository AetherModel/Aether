// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <string>
#include <iostream>
#include <chrono>

#include "aether.h"

// -----------------------------------------------------------------------
// Initialize class Report
// -----------------------------------------------------------------------

Report::Report() {
  current_entry = "";
  nEntries = 0;
  iVerbose = -2;
  divider = ">";
  divider_length = divider.length();
  // Set iLevel to -1, so that the call in main takes it to 0:
  iLevel = -1;
}

// -----------------------------------------------------------------------
// When the code enters a function, report this if the verbose level
// is high enough, and record the system start time of the entry, so that
// the run-time of the function can be recorded on exit.
// -----------------------------------------------------------------------

void Report::enter(std::string input, int &iFunction) {

  int iOldStrLen = current_entry.length();

  current_entry = current_entry + divider + input;

  int iEntry = -1;

  if (iFunction > -1)
    if (current_entry == entries[iFunction].entry)
      iEntry = iFunction;

  if (iEntry == -1) {
    for (int i = 0; i < nEntries; i++)
      if (current_entry == entries[i].entry)
        iEntry = i;
  }

  if (iEntry == -1) {
    item_struct tmp;
    tmp.entry = current_entry;
    tmp.nTimes = 0;
    tmp.timing_total = 0.0;
    tmp.iStringPosBefore = iOldStrLen;
    tmp.iLastEntry = iCurrentFunction;
    entries.push_back(tmp);
    nEntries++;
    iEntry = nEntries - 1;
    iFunction = iEntry;
  }

  // This was taken from
  // https://stackoverflow.com/questions/19555121/how-to-get-current-timestamp-in-milliseconds-since-1970-just-the-way-java-gets
  unsigned long long now = std::chrono::duration_cast<std::chrono::milliseconds>
                           (std::chrono::system_clock::now().time_since_epoch()).count();

  entries[iEntry].timing_start = now;
  iLevel++;
  entries[iEntry].iLevel = iLevel;
  iCurrentFunction = iEntry;

  // Sometimes it is good to uncomment this line and see what is happening:
  //std::cout << "iLevel : " << iLevel << " " << current_entry << "\n";
  print(iLevel, "Entering function : " + current_entry);
}

// -----------------------------------------------------------------------
// This records the exit of the exit of the function, computing the
// runtime of the function.
// This assumes that the enter and exit functions are perfectly paired,
// so that the enter function sets the iCurrentFunction variable.
// -----------------------------------------------------------------------

void Report::exit(std::string input) {

  int iEntry = -1;
  iEntry = iCurrentFunction;

  if (iEntry == -1) {
    std::cout << "Report::exit Error!!! Could not find valid entry!\n";
    std::cout << "current_entry : " << current_entry << "\n";
  } else {

    if (DoReportOnExit)
      print(iLevel, "Exiting function : " + current_entry);

    // Get current system time:
    unsigned long long now = std::chrono::duration_cast<std::chrono::milliseconds>
                             (std::chrono::system_clock::now().time_since_epoch()).count();

    // Calculate the difference in times, do get the total timing:
    entries[iEntry].timing_total = entries[iEntry].timing_total +
                                   float(now - entries[iEntry].timing_start) / 1000.0;

    // Increment the total number of times that the function has been called:
    entries[iEntry].nTimes++;

    // Pop the current function off the stack:
    current_entry = current_entry.substr(0, entries[iEntry].iStringPosBefore);
    iCurrentFunction = entries[iEntry].iLastEntry;
    iLevel--;
  }
}

// -----------------------------------------------------------------------
// Loop through all reported functions and report their run times
// -----------------------------------------------------------------------

void Report::times() {
  if (iVerbose >= 0) {
    std::cout << "Timing Summary :\n";
    float min_timing = entries[0].timing_total * TimingPercent / 100.0;

    for (int i = 0; i < nEntries; i++) {
      if (entries[i].iLevel <= iTimingDepth &&
          entries[i].timing_total >= min_timing) {
        std::cout << entries[i].entry << "\n";

        for (int j = 0; j < entries[i].iLevel; j++)
          std::cout << "  ";

        std::cout << "nTimes called : " << entries[i].nTimes << "\n";

        for (int j = 0; j < entries[i].iLevel; j++)
          std::cout << "  ";

        std::cout << "timing_total (s) : " << entries[i].timing_total << "\n";
      }
    }
  }
}

// -----------------------------------------------------------------------
// Print string if verbose level is set high enough
// -----------------------------------------------------------------------

void Report::print(int iLevel, std::string output_string) {
  if (test_verbose(iLevel))
    std::cout << output_string << "\n";
}

// -----------------------------------------------------------------------
// Test verbose level and print line starter if it is high enough
// -----------------------------------------------------------------------

int Report::test_verbose(int iLevel) {
  int iPass = 0;

  if (iLevel <= iVerbose) {
    iPass = 1;

    for (int iL = 0; iL < iLevel; iL++)
      std::cout << "=";

    std::cout << "> ";
  }

  return iPass;
}

// -----------------------------------------------------------------------
// Set the verbose level in the code.
// -----------------------------------------------------------------------

void Report::set_verbose(int input) {
  iVerbose = input;
}

// -----------------------------------------------------------------------
// Set the depth to report for timing at the end of the run
// -----------------------------------------------------------------------

void Report::set_timing_depth(int input) {
  iTimingDepth = input;
}

// -----------------------------------------------------------------------
// Set the percent to report for timing at the end of the run
// -----------------------------------------------------------------------

void Report::set_timing_percent(float input) {
  TimingPercent = input;
}

// -----------------------------------------------------------------------
// Get the verbose level
// -----------------------------------------------------------------------

int Report::get_verbose() {
  return iVerbose;
}

// -----------------------------------------------------------------------
// Report student checker
// -----------------------------------------------------------------------

void Report::student_checker_function_name(bool isStudent,
					   std::string cStudentName,
					   int iFunctionNumber,
					   std::string cFunctionName) {
  if (isStudent) {
    if (cFunctionName.length() > 1)
      print(-1, cStudentName + " found function " + cFunctionName);
    else
      std::cout << "> (" << iFunctionNumber << ")"
		<< " What function is this: " << cStudentName << "?\n";
  }
  return;
}
