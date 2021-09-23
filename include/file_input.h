// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_FILE_INPUT_H_
#define INCLUDE_FILE_INPUT_H_

/**************************************************************
 * Functions that assist in the reading of files.
 *
 **************************************************************/

#include <string>
#include <fstream>
#include <vector>

/**************************************************************
 \brief Reads the file until it finds a #command, returns "#command"
 \param file_ptr file pointer to the open file
 **/
std::string find_next_hash(std::ifstream &file_ptr);

/**************************************************************
 \brief Reads a series of line that can be in a CSV format, returns 2D array
 \param file_ptr file pointer to the open file
 **/
std::vector<std::vector<std::string>> read_csv(std::ifstream &file_ptr);

/**************************************************************
 \brief Reads a series of line that can be in a SSV format, returns 2D array
 \param file_ptr file pointer to the open file
 \param
 **/
std::vector<std::vector<std::string>> read_ssv(std::ifstream &file_ptr);

/**************************************************************
 \brief Reads either a comma-separated time or series of lines describing time

  format is either 
  y, m, d, h, m, s, ms 
  or
  y
  m
  d
  h
  m
  s
  ms

 \param file_ptr file pointer to the open file
 \param
 **/
std::vector<int> read_itime(std::ifstream &file_ptr, std::string hash);

/**************************************************************
 \brief Converts a string to lower case, returning lower case string
 \param instring string to make lower
 **/
std::string make_lower(std::string instring);

/**************************************************************
 \brief Strips everything including and after two spaces, returns stripped str
 \param instring string to strip end off of
 **/
std::string strip_string_end(std::string instring);

/**************************************************************
 \brief Takes all of the spaces out of the sting, returns stripped string
 \param instring string to take spaces out of
 **/
std::string strip_spaces(std::string instring);

/**************************************************************
 \brief read in a string from a file
 \param file_ptr file pointer to the open file
 \param hash string to print out if there is an error
 **/
std::string read_string(std::ifstream &file_ptr, std::string hash);

/**************************************************************
 \brief read in an integer from a file
 \param file_ptr file pointer to the open file
 \param hash string to print out if there is an error
 **/
int read_int(std::ifstream &file_ptr, std::string hash);

/**************************************************************
 \brief read in a float from a file
 \param file_ptr file pointer to the open file
 \param hash string to print out if there is an error
 **/
float read_float(std::ifstream &file_ptr, std::string hash);

/**************************************************************
 \brief take a string with comma-separated-values and parse into vector
 \param line string containing comma separated values
 **/
std::vector<std::string> parse_csv_row_into_vector(std::string line);

#endif  // INCLUDE_FILE_INPUT_H_
