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

  Logfile(std::string logfileNameIn,
	  precision_t dtIn,
	  bool doAppend,
	  Indices indices,
	  Inputs inputs,
	  Report &report);

  bool write_logfile(Times time,
		     Neutrals neutrals,
		     Ions ions,
		     Inputs inputs,
		     Indices indices,
		     Report &report);

  void close_logfile();

private:

  std::vector<precision_t> lla {5,4,40};

  std::ofstream logfilestream;

  bool doAppend = false;
  bool isOk = true;
  std::string logfileName;
  precision_t dt;
  
};

#endif  // INCLUDE_LOGFILE_H_
