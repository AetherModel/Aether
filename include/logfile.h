// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_LOGFILE_H_
#define INCLUDE_LOGFILE_H_

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
#include "aether.h"
#include <fstream>
#include <iostream>

using namespace std;

class Logfile {

public:

  Logfile(Indices indices, Inputs inputs, Report &report);

  void write_logfile(Times time,
		     Neutrals neutrals,
		     Ions ions,
		     Inputs inputs,
		     Indices indices,
		     Report report);

  void close_logfile();
private:
  vector<precision_t> lla {5,4,40};
  ofstream logfilestream;  
  bool header = false;

};

#endif  // INCLUDE_LOGFILE_H_
