// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_INDICES_H_
#define INCLUDE_INDICES_H_

/**************************************************************
 * \class Indices
 *
 * \brief A class for keeping track of indices (1d vectors w/time)

 This is a class that allows users to read in and keep track of
 indices, such as IMF Bz, Kp, AE, F107, etc. These are often needed
 to drive other sub models, such as EUVAC, MSIS, Weimer, etc. Basically,
 the way you use this is to:
 - read a file (custom for each type of file)
 - get the proper index number
 - call the set_index function with the time and values array
 - call the appropriate get function with time to get index value at time

 * \author Aaron Ridley
 *
 * \date 2021/04/16 
 **************************************************************/

#include <vector>
#include <string>

/// A structure that is made available to the user to allow the
/// general reading of a file, since the information that is needed
/// from the file could all be contained within this structure

struct index_file_output_struct {

  /// number of times read in:
  int64_t nTimes;
  
  /// array of times that correspond to the values:
  std::vector<double> times;
  
  /// number of variables read in:
  int nVars;

  /// variable names as a vector of strings:
  std::vector<std::string> var_names;

  /// a 2d vector (vars vs times) of indices values:
  std::vector<std::vector<float>> values;

  /// a vector of missing values for each variable:
  std::vector<float> missing_values;

  /// The index_id returned by the call to get_XXX_index_id:
  std::vector<int> index_id;
};

/**************************************************************
 \brief Print out information in the structure for debugging
 \param contents structure containing the index file output
 **/
void print_index_file_output_struct(index_file_output_struct contents);

class Indices {

// -----------------------------------------------------------------------
// Public functions and variables
// -----------------------------------------------------------------------

 public:

  /**************************************************************
   \brief Initialize the class
   **/
  Indices();

  /**************************************************************
   \brief get the daily f107 value at the given time
   \param time time in seconds
   **/
  precision_t get_f107(double time);

  /**************************************************************
   \brief get the 81-day average f107 at the given time
   \param time time in seconds
   **/
  precision_t get_f107a(double time);

  /**************************************************************
   \brief a series of functions that return the internal index number

   In order to keep track of which index is which, the class uses 
   constants. These functions return these constants. The user doesn't
   really need to know about the constants, but they have to get the
   constant (when reading the file, for example) and then provide that 
   to the set index function. Conversely, we could create a bunch of
   set_ functions (such as the set_f107 function below). We figured
   that this minor inconvience is easier than making a bunch of set_
   functions.

   **/
  int get_f107_index_id();
  int get_f107a_index_id();
  int get_imfbx_index_id();
  int get_imfby_index_id();
  int get_imfbz_index_id();
  int get_swvx_index_id();
  int get_swvy_index_id();
  int get_swvz_index_id();
  int get_swn_index_id();
  int get_swt_index_id();
  int get_ae_index_id();
  int get_au_index_id();
  int get_al_index_id();

  /**************************************************************
   \brief Return the indices index of the variable name
   \param name the name of the variable to find the index for
   **/
  
  int lookup_index_id(std::string name);
  
  /**************************************************************
   \brief This function sets the f107, does an 81 day ave, sets f107a too
   \param f107_contents contents from the f107 file (time, f107, etc.)
   **/
  // This is the method for setting f107 specifically:
  void set_f107(index_file_output_struct f107_contents);

  /**************************************************************
   \brief set the index array into the indices class
   \param index_id which index is being checked in (bx, by, kp, ae, etc)
   \param time vector of time for each index value
   \param values vector of values for each index value
   \param missing value for missing data
   **/
  bool set_index(int index_id,
		 std::vector<double> time,
		 std::vector<float> values,
		 precision_t missing);

  /**************************************************************
   \brief set the index array into the indices class
   \param index_name which index is being checked in (bx, by, kp, ae, etc)
   \param time vector of time for each index value
   \param values vector of values for each index value
   \param missing value for missing data
   **/
  bool set_index(std::string index_name,
		 std::vector<double> timearray,
		 std::vector<float> indexarray,
		 precision_t missing);
  
  /**************************************************************
   \brief Perturbs the indices requested by user input
   **/
  bool perturb();

  /**************************************************************
   \brief Perturbs the specific indices based on the user input
   \param iIndex which index to perturb
   \param seed random seed for perturbations
   \param style characteristics of the perturbations (+/*, mean/std, const)
   \param DoReport output information if true
   **/
  void perturb_index(int iIndex, int seed, json style, bool DoReport);

  /**************************************************************
   \brief The general function that returns the index value at the time
   \param time the time in seconds that the index is requested at
   \param the index to return (i.e., one of the constants defined above)
   **/
  precision_t get_index(double time, int index);

  /**************************************************************
   * \brief Get the name of the indices at the specified index
   * \param iIndex which index to get name
   * \return The string of name if the function succeeds, empty string if iIndex is out of range
   **/
  std::string get_name(int iIndex); 

  /**************************************************************
   \brief Return the number of the indices vector
   **/
  int all_indices_array_size();

// -----------------------------------------------------------------------
// Private functions and variables
// -----------------------------------------------------------------------

private:

  /// structure containing information about the specific index:
  struct index_time_pair {

    /// the number of values in the vectors:
    int64_t nValues;

    /// a vector of values for the index:
    std::vector<precision_t> values;

    /// a vector of times for the values:
    std::vector<double> times;

    /// the name of the index as a string:
    std::string name;
  };

  /// the vector that contains all of the indices vectors:
  std::vector<index_time_pair> all_indices_arrays;

  /// constants for keeping track of indices:
  const int iF107_ = 0;
  const int iF107A_ = 1;
  const int iIMFBX_ = 2;
  const int iIMFBY_ = 3;
  const int iIMFBZ_ = 4;
  const int iSWVX_ = 5;
  const int iSWVY_ = 6;
  const int iSWVZ_ = 7;
  const int iSWN_ = 8;
  const int iSWT_ = 9;
  const int iAE_ = 10;
  const int iAL_ = 11;
  const int iAU_ = 12;
  /// number of indices that system is capable of keeping track of:
  int nIndices = 13;

  /// this will let us go back and forth between names and ids:
  json indices_lookup;

};


#endif // INCLUDE_INDICES_H_
