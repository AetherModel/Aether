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

/**
 * The class Satellite is used to track the satellites
 * Given any time, the user can obtain the geographic location of the satellite
 * 
 * ASSUMPTION : The satellite csv layout is the same as the following
 * year   mon   day   hr    min   sec   lon      lat      alt   x   y    z    vx     vy     vz
 * (int)  (int) (int) (int) (int) (int) (degree) (degree) (km) (km) (km) (km) (km/s) (km/s) (km/s)
 */

class Satellite {

public:

  /**
   * \brief Initialize the satellite class
   *        The name of the satellite is not allowed to have any characters which can 
   *        terminate the read of a string including white space' ', endline'\n', and '\t'
   *        Different satellites must have different names (not only input file names)
   * \param csv_in The path to the satellite csv file
   * \param log_name_in The path to write the satellite log file
   * \param sat_header The first line of the log file
   * \param dt_in The time gate
   */
  Satellite(const std::string &csv_in,
            const std::string &log_name_in,
            const std::string &sat_header,
            const precision_t dt_in);

  /**
   * \brief Get the position of the satellite at given time
   * \param lons_out Vector to stroe the longitude of the location
   * \param lats_out Vector to store the latitude of the location
   * \param alts_out Vector to store the altitude of the location
   * \param time_in Time to determine location
   */
  void get_position(std::vector<precision_t> &lons_out,
                    std::vector<precision_t> &lats_out,
                    std::vector<precision_t> &alts_out,
                    const double time_in) const;

  // Get the name of the satellite
  std::string get_name() const;

  // Get the time gate of the satellite
  precision_t get_dt() const;

  // Write the content to the satellite log
  void write_log(const std::vector<int> &iCurrent,
                 const std::vector<precision_t> &variables);

  // DEBUG
  void print();

private:

  // The name of the satellite
  std::string name;
  // The time gate
  precision_t dt;
  // The file name to write log
  std::string log_name;
  // The time should always be ascending
  // The type for time is double rather than precision_t
  std::vector<double> timereals;
  // The position of the satellite at those times
  std::vector<precision_t> lons;
  std::vector<precision_t> lats;
  std::vector<precision_t> alts;
  // The rest columns are not used (i.e. from x to vz)
};

/**
 * The class to maintain the logfile
 */

class Logfile {

public:

  /**
   * \brief Initialize the Logfile.
   *    The logfile will output all indicies and specified neutrals and ions
   *    every dt time.
   */
  Logfile(Indices &indices);
  
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
                     Times &time);

private:

  // The name of logfile
  std::string logfileName;
  // The file stream to write
  std::ofstream logfilestream;
  // The specified variables
  std::vector<std::string> species;
  // The satellites as a vector
  std::vector<Satellite> satellites;
  // The time gate
  precision_t dt;
  // Whether append or rewrite
  bool doAppend;

  // A randomly chosen point for test
  std::vector<precision_t> lla {5,4,40};
};

#endif  // INCLUDE_LOGFILE_H_
