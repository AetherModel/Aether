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

  /**
   * \brief Initialize the Logfile. The logfile will output
   *    all indicies and specified neutrals and ions every dt time
   */
  Logfile(Indices &indices,
          Inputs &input,
          Report &report);
  
  /**
   * \brief Close the file stream if not append
   */
  ~Logfile();

  /**
   * \brief Check the time gate, and write the values into log file
   *    if needed
   */
  bool write_logfile(Indices &indices,
                     Neutrals &neutrals,
                     Ions &ions,
                     Grid &gGrid,
                     Times &time,
                     Report &report);

private:

  // The name of logfile
  std::string logfileName;
  // The file stream to write
  std::ofstream logfilestream;
  // The specified variables
  std::vector<std::string> species;
  // The time gate
  precision_t dt;
  // Whether append or rewrite
  bool doAppend;

  bool isOk;

  // A randomly chosen point for test
  std::vector<precision_t> lla {5,4,40};
};

#endif  // INCLUDE_LOGFILE_H_
