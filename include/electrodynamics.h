// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_ELECTRODYNAMICS_H_
#define INCLUDE_ELECTRODYNAMICS_H_

/**************************************************************
 * This is the electrodynamics class for Aether.

   To use this, you have to:
   1. initialize it, which reads in the electrodynamics file (if used).
   2. set the time to use.
   3. set any indices to you need for the electrodynmics.
   4. set the magnetic latitudes (2d) you need
   5. set the magnetic local times (2d) you need
   6. call get_potential or get_eflux or ... to get the electrodynamics
      on the grid that was specified.
   2-6 can all be called as many times as you want.
 *
 **************************************************************/

#include <vector>

#include "aether.h"

class Electrodynamics {

 public:

  /**************************************************************
     \brief Initialize electrodynamics variables and routines

     This does the following:
     - initialize all variables to missing values
     - read in file if it exists

     \param input Need to pass Input class, so code can get info
                  about how user has configured things. Inside
                  input, the function uses .efield_model, .auroral_model,
                  .electrodynamics_file.

     \param report Need to pass Report class, so reporting can occur  
   **/
  Electrodynamics(Inputs input, Report &report);

  /**************************************************************
     \brief used in main.cpp to ensure electrodynamics times and
      aether input times match up. Returns false if times are
      misaligned, true if they are aligned

     \param inputStartTime input file starting time

     \param inputEndTime input file ending time
   **/

  bool check_times(double inputStartTime, double inputEndTime);
  
  /**************************************************************
     \brief used in advance.cpp to get potential, eflux, avee


     \param magLat magnetic latitude

     \param magLocalTime magnetic local time

     \param report reporting
   **/

  std::tuple<fcube, fmat, fmat> get_electrodynamics(fcube magLat, fcube magLocalTime, Report &report);

  /**************************************************************
     \brief Gets interpolation indices

     Performs 2d interpolation over search vector to get indices

     \param vals the 2d array that needs indices
     \param search The vector of values to interpolate over  
   **/

  fmat get_interpolation_indices(fmat vals, fvec search);

  /**************************************************************
     \brief Sets time needed for electrodynamics

     Internally, if there is a file read, this function:
     - finds interpolation indices in time
     - interpolates the primaary quantities (pot, avee, eflux, etc) to
       the given time, creating internal states of these quantities. This
       is done for all relavant grids.

     \param time the time requested.
     \param report Need to pass Report class, so reporting can occur  
   **/
  
  void set_time(double time, Report &report);

  /**************************************************************
     \brief Sets the current grid to request data on

     Internally, if there is a file read, this function finds the
     interpolation indices in space for each relevant grid

     \param lats a 2D matrix of magnetic latitudes to interpolate to
     \param mlts a 2D matrix of magnetic local times to interpolate to
     \param report Need to pass Report class, so reporting can occur  
   **/
  void set_grid(fmat lats, fmat mlts, Report &report);

  /**************************************************************
     \brief Set the IMF Bx for internal usage
     \param value Value to assign to IMF Bx (nT)
   **/
  void set_imf_bx(float value);

  /**************************************************************
     \brief Set the IMF By for internal usage
     \param value Value to assign to IMF By (nT)
   **/
  void set_imf_by(float value);

  /**************************************************************
     \brief  Set the IMF Bz for internal usage
     \param value Value to assign to IMF Bz (nT)
   **/
  void set_imf_bz(float value);

  /**************************************************************
     \brief Set the Solar Wind Velocity for internal usage
     \param value Value to assign to Solar wind velocity (km/s)
   **/
  void set_sw_v(float value);

  /**************************************************************
     \brief Set the Solar Wind Density for internal usage
     \param value Value to assign to Solar Wind Density (/cc)
   **/
  void set_sw_n(float value);

  /**************************************************************
     \brief Set the Hemispheric Power for internal usage
     \param value Value to assign to Hemispheric Power (GW)
     \param
   **/
  void set_hp(float value);

  /**************************************************************
     \brief Set the AU index for internal usage
     \param value Value to assign to Auroral Upper Index (nT)
   **/
  void set_au(float value);

  /**************************************************************
     \brief Set the AL index for internal usage
     \param value Value to assign to Auroral Lower Index (nT)
   **/
  void set_al(float value);

  /**************************************************************
     \brief Set the AE Index for internal usage
     \param value Value to assign to Auroral Electrojet Index (nT)
   **/
  void set_ae(float value);

  /**************************************************************
     \brief Set the Kp Index for internal usage
     \param value Value to assign to Kp index
   **/
  void set_kp(float value);
  
  /**************************************************************
     \brief Get 2D electric potential on specified grid

     This function returns the electric potential on the requested
     grid at the requested time (with the requested indices, if
     applicable)

     This function does the following:
     - creates an empty potential matrix ok, I see to return
     - Loops through the grids in priority order calling set_values
       with the potentials in the grids

     \param report Need to pass Report class, so reporting can occur  
   **/
  fcube get_potential(fcube magLat, fcube magLocalTime, Report &report);

  /**************************************************************
     \brief Get 2D electron energy flux on specified grid

     This function returns the electron energy flux on the requested
     grid at the requested time (with the requested indices, if
     applicable)

     This function does the following:
     - creates an empty eflux matrix to return
     - Loops through the grids in priority order calling set_values
       with the eflux in the grids

     \param report Need to pass Report class, so reporting can occur  
   **/
  fmat get_eflux(fcube magLat, fcube magLocalTime, Report &report);

  /**************************************************************
     \brief Get 2D electron average energy on specified grid

     This function returns the electron average energy on the requested
     grid at the requested time (with the requested indices, if
     applicable)

     This function does the following:
     - creates an empty avee matrix to return
     - Loops through the grids in priority order calling set_values
       with the avee in the grids

     \param report Need to pass Report class, so reporting can occur  
   **/
  fmat get_avee(fcube magLat, fcube magLocalTime, Report &report);

  /**************************************************************
     \brief Get 2D ion energy flux on specified grid

     This function returns the ion energy flux on the requested
     grid at the requested time (with the requested indices, if
     applicable)

     This function does the following:
     - creates an empty ion eflux matrix to return
     - Loops through the grids in priority order calling set_values
       with the ion eflux in the grids

     \param report Need to pass Report class, so reporting can occur  
   **/
  fmat get_ion_eflux(Report &report);

  /**************************************************************
     \brief Get 2D ion average energy on specified grid

     This function returns the ion average energy on the requested
     grid at the requested time (with the requested indices, if
     applicable)

     This function does the following:
     - creates an empty ion avee matrix to return
     - Loops through the grids in priority order calling set_values
       with the ion avee in the grids

     \param report Need to pass Report class, so reporting can occur  
   **/
  fmat get_ion_avee(Report &report);
  
 private:

  /// This is the interpolation method for time:
  int iTimeInterpolationMethod;
  /// For all iteration methods, code should return first value if
  /// before first time, last value if after last time.
  /// Use the previous value:
  const int iPrevious_ = 1;
  /// Use the next value:
  const int iNext_ = 2;
  // Use the closest value:
  const int iClosest_ = 3; 
  /// Interpolate:
  const int iInterp_ = 4;

  /// The input file to read that contains the grid of electrodynamics:
  std::string input_file;

  /// Variables needed by user - all set by set_ functions:

  /// This is the time needed:
  double time_needed;

  /// A 2d array of magnetic latitude needed. Can set interpolation
  /// coefficients in all of the grids when this is called:
  fmat lats_needed;

  /// A 2d array of magnetic local times needed. Can set interpolation
  /// coefficients in all of the grids when this is called:
  fmat mlts_needed;
  
  /// These are all indices that may be needed by sub-models:
  float imf_bx_needed;
  float imf_by_needed;
  float imf_bz_needed;
  float sw_v_needed;
  float sw_n_needed;
  float hp_needed;
  float au_needed;
  float al_needed;
  float ae_needed;
  float kp_needed;

  /// If we read in an electrodyanmics file, set this to 1. Otherwise,
  /// set it to 0:
  int iUseGridBasedModel;

  /// If we don't read in an electrodynamics file, then this should be
  /// set to an electric field model to use.  Need to add model types.
  std::string efield_model_to_use;

  /// If we don't read in an electrodynamics file, then this should be
  /// set to an auroral model to use.  Need to add model types.
  std::string auroral_model_to_use;
  
  /// Set the interpolation indices as a float. For each interpolation index,
  /// the integer portion is the current index, and the decimal part is the 
  /// percentage of the distance between the current index and the next
  /// index.  For example, a distance midway between index 45 and 46 
  /// would give an interpolation index of 45.5.
  /// For time, we are assuming that all grids have the same times or that
  /// there are no overlaps in time, I think.
  float time_index;

  /**************************************************************
   * input_electrodynamics_struct is a structure that contains
   * information about the electrodynamics on a magnetic latitude /
   * local time grid. It is assumed that the quantities are on a
   * regular grid, in that the lat and mlt can be described with 1D
   * arrays and the quantities can be described with 2D arrays.
   * Time then makes these matrices into vectors. The quantities
   * stored are energy flux, average energy, electric potential,
   * ion energy flux, ion average energy.
   **/
  struct input_electrodynamics_struct {

    /// Number of latitudes and magnetic local times:
    int nLats;
    int nMlts;

    /// Vectors of magnetic latitudes and magnetic local times:
    fvec mlats;
    fvec mlts;

    /// Vector of times (in seconds):
    std::vector<double> times;

    /// Vector of 2d electric potentials (in Volts):
    std::vector<fmat> potential;

    /// Potential at current time:
    fmat potential_current;
    
    /// Vector of 2d electron energy flux (in ergs/cm2/s):
    std::vector<fmat> energy_flux;
    /// Said energy flux at the current time:
    fmat energy_flux_current;
    
    /// Vector of 2d electron average energy (in keV):
    std::vector<fmat> average_energy;
    /// Average energy at current time:
    fmat average_energy_current;
    
    /// Vector of 2d ion energy flux (in ergs/cm2/s):
    std::vector<fmat> ion_energy_flux;
    /// ion energy flux at current time:
    fmat ion_energy_flux_current;
    
    /// Vector of 2d ion average energy (in keV):
    std::vector<fmat> ion_average_energy;
    /// ion average energy at current time:
    fmat ion_average_energy_current;

    /// Set to 1 if ion precipitation is included, else set to 0:
    int DoesIncludeIonPrecip;
    
    /// This sets the priority of the grid. The higher the number, the
    /// more important it is, so it should overwrite any regions of
    /// a lower priority grid. For example, you could have a global
    /// grid and a regional grid, whre the global grid is assigned a
    /// priority of 1, and the regional grid is assigned a priority of 2,
    /// so that the regions inside the regional grid will overwrite
    /// the global grid:
    int grid_priority;

    /// interpolation indices for latitudes. If the requested latitude
    /// is outside of the latitude range of the grid, then the
    /// interpolation index should be set to -1:
    fmat lats_indices;
    /// interpolation indices for mlts. If the requested mlt
    /// is outside of the mlt range of the grid, then the
    /// interpolation index should be set to -1:
    fmat mlts_indices;
    
  };
  
  /// As described above, a structure containing the grid-based
  /// values of electrodynamics as a function of time.  This is
  /// vector, because we can have nested grids, or, in theory, the
  /// grid could change as a function of time. You can then search
  /// for the apropriate grid in space and time.
  std::vector<input_electrodynamics_struct> input_electrodynamics;
  
  /// Because each grid has a priority, we need to go through them in
  /// priority order, this is the sorted indices list, so that
  /// grid_order[0] points to the input_electrodynamics with the
  /// lowest priority, grid_order[1] points to the 2nd lowest, etc.
  std::vector<int> grid_order;
  
  /// Number of input grids for electrodynamics:
  int nElectrodynamicsGrids;
  
  /**************************************************************
     \brief Reads a netcdf file that has the electrodynamics specification

     Reads a netcdf file that contains at a minimum:
     - Potential
     - Auroral energy flux
     - Auroral average energy
     May contain:
     - Ion energy flux
     - Ion average energy

     These are on a magnetic lat/mlt grid as a function of time, and
     should be put into the input_electrodynamics structure.

     \param filename
     \param report Need to pass Report class, so reporting can occur
   **/
  void read_netcdf_electrodynamics_file(std::string filename,
					Report &report);

  /**************************************************************
     \brief Takes the pot/eflux/avee/etc and interpolates to the grid

     This function takes values from the input_electrodynamics
     structure and interpolates the values onto the user specified
     grid. The tricky bit is that there can be multiple overlapping
     grids.  So, this needs to be called for each of the existing
     grids, so that the values are overwritten. To keep it
     "functional", we pass in the last round of values and those are
     moved into the output values and then the overlapping region is
     overwritten (e.g., in the get_potential function, the 
     grids need to be cycled through calling get_values with the
     potential on that grid and the interpolation indices for the grid.

     \param values_current the pot/eflux/avee/etc from
     input_electrodynamics grid
 
     \param lats_indices the interpolation indices for the current
     grid latitudes

     \param mlts_indices the interpolation indices for the current
     grid mlts

     \param values_old the output of this function for the last grid
  **/  
  
  fmat get_values(fmat matToInterpolateOn, int rows, int cols);
};

#endif // INCLUDE_ELECTRODYNAMICS_H_
