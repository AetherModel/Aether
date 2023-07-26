// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_OUTPUT_H_
#define INCLUDE_OUTPUT_H_

#include "aether.h"

/**************************************************************
 * \class Output
 * \brief A containing to allow storage of variables for output
 * 
 * Writing output is a multi-step process now:
 *  1. Create a container to store the variables you want to output
 *  2. Define the variables to output within the container
 *  3. Store the variables
 *  4. Write the output
 *
 * \author Aaron Ridley
 * \date 2021/10/21 
 **************************************************************/

class OutputContainer {

 public:
  
  /**********************************************************************
     \brief initialize the output container
   **/
  OutputContainer();

  /**********************************************************************
     \brief set the output type to binary
   **/
  void set_binary();

  /**********************************************************************
     \brief set output type to be netcdf
   **/
  void set_netcdf();

  /**********************************************************************
     \brief set output type to be hdf5
   **/
  void set_hdf5();

  /**********************************************************************
     \brief set the directory to output to
     \param in_dir the directory name
   **/
  void set_directory(std::string in_dir);

  /**********************************************************************
     \brief set the filename to output to
     \param in_filename the file name
   **/
  void set_filename(std::string in_filename);

  /**********************************************************************
     \brief store a variable to the list of variables to output
     \param name name of the variable to output
     \param unit units of the output
     \param value the array of the data to output
   **/
  void store_variable(std::string name,
		      std::string unit,
		      arma_cube value);

  /**********************************************************************
     \brief store a variable to the list of variables to output
     \param name name of the variable to output
     \param long_name longer name of the variable to output
     \param unit units of the output
     \param value the array of the data to output
   **/
  void store_variable(std::string name,
		      std::string long_name,
		      std::string unit,
		      arma_cube value);

  /**********************************************************************
     \brief Get an arma_cube from the Container
     \param iElement return the iElement arma_cube
   **/
  arma_cube get_element_value(int64_t iElement);

  /**********************************************************************
     \brief Get an arma_cube from the Container
     \param var_to_get variable name to return the arma_cube
   **/
  arma_cube get_element_value(std::string var_to_get);

  /**********************************************************************
     \brief Get the variable name from the Container
     \param iElement return the iElement name
   **/
  std::string get_element_name(int64_t iElement);

  /**********************************************************************
     \brief Get number of elements in output container
   **/
  int64_t get_nElements();

  /**********************************************************************
     \brief Get element number for specified name
     \param var_to_find variable to find in the list
   **/
  int64_t find_element(std::string var_to_find);

  /**********************************************************************
     \brief set the time for the output
     \param time the current time for the simulation
   **/
  void set_time(double time);

  /**********************************************************************
     \brief set the version of Aether
     \param in_version the version number
   **/
  void set_version(float in_version);

  /**********************************************************************
     \brief write a file with the information in the container
   **/
  void write();
  
  /**********************************************************************
     \brief write a json header file with the information in the container
   **/
  int write_container_header();
  
  /**********************************************************************
     \brief write a binary file with the information in the container
   **/
  int write_container_binary();
  
  /**********************************************************************
     \brief write a netcdf file with the information in the container
   **/
  int write_container_netcdf();
  
  /**********************************************************************
     \brief read from a file an load into the container
   **/
  void read();
  
  /**********************************************************************
     \brief display information contained in the container
   **/
  void display();
  
  /**********************************************************************
     \brief read a netcdf file - put the information in the container
   **/
  int read_container_netcdf();
  
  /**********************************************************************
     \brief clears the vector of variables 
   **/
  void clear_variables();
  
 private:

  /// User can set the directory for output
  std::string directory;

  /// Filenames for basic types are standard, but can be set by the user.
  /// The filename will be appended with _YYYYMMDD_HHMMSS.ext
  std::string filename;

  /// This is the container variable information
  struct var_struct {
    /// Name of the variable to output
    std::string cName;
    /// Long Name of the variable to output
    std::string cLongName;
    /// Unit of the variable to output
    std::string cUnit;
    /// 3D array of the variable to output (this can be N x M x L, or
    /// 1 x M x L, or whatever, as long as it is 3D, any of the
    /// dimensions can be 1):
    arma_cube value;
  };
  /// The vector containing all of the variable information:
  std::vector<var_struct> elements;

  /// The time of the data
  std::vector<int> itime;
  /// The version of the code / data / whatever:
  float version;

  /// The frequency of the output for this particular container:
  float dt_output;
  
  /// This is to allow the user to select different output formats
  int output_type;

  /// These are the supported output types:
  const int binary_type = 0;
  const int netcdf_type = 1;
  const int hdf5_type = 2;
  
};

/**********************************************************************
   \brief Fills output containers and outputs them for common output types

   This function takes user requested output types and fullfills them
   if it understands the types. The supported output types are:

   states - neutral/ion densities, bulk velocities, bulk
            temperatures, and electric potential
   neutrals - output all neutral densities, velocities, and temperature
   ions - ion densities, (advected) ion velocities, ion temperatures,
          electron temperature, electric potential
   bfield - magnetic coordinates, magnetic field vector

   \param neutrals all of the states for the neutrals
   \param ions all of the states for the ions
   \param grid The grid to define the neutrals on
   \param time information about the current time
   \param planet information about the planet
**/

int output(const Neutrals &neutrals,
	   const Ions &ions,
	   const Grid &grid,
	   Times time,
	   const Planets &planet);

void output_binary_3d(std::ofstream &binary,
		      arma_cube value);

#endif  // INCLUDE_OUTPUT_H_
