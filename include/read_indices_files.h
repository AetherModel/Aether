// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_READ_INDICES_H_
#define INCLUDE_READ_INDICES_H_

/**************************************************************
 * A set of routines to read in a number of different indices
 * files, including the NGDC F10.7 files, and OMNIWeb files.
 **************************************************************/

#include <vector>
#include <string>

#include "aether.h"

/**********************************************************************
  \brief Reads in all of the indices files and stores them in Indices

  This function goes through all of the input indices files and 
  reads in the files, then stores the values into the Indices class.
  At this point, it can read in the following file types:
  1. NGDC F10.7 files.
  2. OMNIWeb files that include IMF, Solar Wind, and/or AE/AU/AL.

  \param indices the indices class
**/
bool read_and_store_indices(Indices &indices);

/**********************************************************************
  \brief Read the NGDC F10.7 file format and store in the index_file struct
  \param f107_file the f10.7 file to read in
**/
index_file_output_struct read_f107_file(std::string f107_file,
					Indices indices);

/**********************************************************************
  \brief Read the OMNIWeb file format and store in the index_file struct
  \param omni_file the onmiweb file to read
  \param indices needed to get the indices index for each variable
**/
index_file_output_struct read_omni_file(std::string omni_file,
					Indices indices);

/**********************************************************************
  \brief This code compares a string to return the variable index
  \param variable_name the name of the variable to match with for OMNIWeb files
  \param indices needed to get the indices index for each variable
**/
int pair_omniweb_index(std::string, Indices index);

/**********************************************************************
  \brief This code returns the missing value for a given OMNIWeb variable
  \param variable_name the name of the variable to get the OMNIWeb missing value
**/
precision_t get_omniweb_missing_value(std::string);

#endif  // INCLUDE_READ_INDICES_H_
