// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_LOGFILE_H_
#define INCLUDE_LOGFILE_H_

/**************************************************************
 * 
 * logfile.h:
 *
 *    Write the logfile
 *
 **************************************************************/

#include "aether.h"
#include <fstream>
#include <iostream>

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

  void append_mode();

private:

  std::vector<precision_t> lla {5,4,40};

  std::ofstream logfilestream;

  bool trunc_mode = true;
};

#endif  // INCLUDE_LOGFILE_H_
