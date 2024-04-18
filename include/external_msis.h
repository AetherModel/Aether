// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_MSIS_H_
#define INCLUDE_MSIS_H_

/**************************************************************
 * \class Msis
 *
 * \brief create an interface to the msis model
 * 
 * MSIS is a neutral model of the atmosphere, written in
 * fortran and provided by NRL. 
 *
 * \author Aaron Ridley
 *
 * \date 2023/04/30 
 *
 **************************************************************/

class Msis {

 public:

  Msis();
  bool set_time(Times time);
  bool set_f107(precision_t f107in, precision_t f107ain);
  bool set_ap(precision_t apin);
  bool set_locations(arma_vec longitude,
		     arma_vec latitude,
		     arma_vec altitude);
  bool set_locations(arma_mat longitude,
		     arma_mat latitude,
		     arma_mat altitude);
  bool set_locations(arma_cube longitude,
		     arma_cube latitude,
		     arma_cube altitude);
  arma_vec get_vec(std::string value);
  arma_mat get_mat(std::string value);
  arma_cube get_cube(std::string value);
  bool is_valid_species(std::string value);
  bool is_ok();
  
 private:

  int iYear, iDay;
  float second;
  float f107 = 80.0, f107a = 80.0;
  float ap = 10.0;

  int nX = -1, nY = -1, nZ = -1, nVars;
  arma_cube lonDeg;
  arma_cube latDeg;
  arma_cube altKm;
  std::vector<arma_cube> msis_results;

  bool didChange = true;  
  json value_lookup;
  bool isCompiled;
  
  bool reset_interface_variable_sizes();
  bool reset_results();
};

#endif // INCLUDE_MSIS_H_
