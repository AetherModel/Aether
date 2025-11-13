// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_IONS_H_
#define INCLUDE_IONS_H_

#include <string>
#include <vector>

/**************************************************************
 * \class Ions
 *
 * \brief Defines the ion states
 *
 * The Ion class defines the ion states as well as a bunch
 * of derived states and source/loss terms.
 *
 * \author Aaron Ridley
 *
 * \date 2021/03/28
 *
 **************************************************************/

class Ions {

 public:

  // This struct contains all of the information needed for a single
  // species of ion.  We will then have a vector of these species.

  int64_t nSpecies = 8;

  struct species_chars {

    /// Name of the species
    std::string cName;

    /// Mass of the species
    precision_t mass;

    /// Charge of the species
    int charge;

    /// Vibrational levels of the species
    int vibe;

    /// Advect or don't advect the species
    bool DoAdvect;

    /// Ion neutral collision frequncies (calculated):
    std::vector<arma_cube> nu_ion_neutral_vcgc;

    /// Ion - Neutral collision frequency coeffient (non-resonant):
    std::vector<precision_t> nu_ion_neutral_coef;

    /// Whether the collision frequency is resonant of not:
    std::vector<bool> nu_is_resonant;

    /// A bunch of parameters for calculating resonant collisions:
    std::vector<precision_t> nu_in_res_temp_min;
    std::vector<precision_t> nu_in_res_coef1;
    std::vector<precision_t> nu_in_res_coef2;
    std::vector<precision_t> nu_in_res_tn_frac;
    std::vector<precision_t> nu_in_res_ti_frac;

    /// Ion - Ion Collision frequencies:
    std::vector<precision_t> nu_ion_ion;

    /// Ion - Electron collision frequencies:
    std::vector<precision_t> nu_ion_electron;

    // Sources and Losses:

    /// Number density of species (/m3)
    arma_cube density_scgc;
    arma_cube newDensity_scgc;

    /// For all below:
    /// Index 0 = x/longitudinal component of velocity
    /// Index 1 = y/latitudinal
    /// Index 2 = z/altitudinal/along the field line

    /// Parallel velocity of species (m/s)
    std::vector<arma_cube> par_velocity_vcgc;

    /// Perpendicular velocity of the species (m/s)
    std::vector<arma_cube> perp_velocity_vcgc;

    /// Total velocity of the species (m/s)
    /// - this is the parallel + perpendicular together
    ///   and can be used on the spherical grid for advection
    std::vector<arma_cube> velocity_vcgc;

    /// Temperature of the given species:
    arma_cube temperature_scgc;

    /// Conduction source term:
    arma_cube conduction_scgc;

    /// Ionization source term:
    arma_cube ionization_scgc;
    /// Total chemical sources and losses:
    arma_cube sources_scgc;
    arma_cube losses_scgc;

    // Heating terms:
    /// Bulk collisional heating with neutrals and electrons (K/s)
    arma_cube heating_neutral_friction_scgc;
    arma_cube heating_electron_friction_scgc;

    /// Bulk collisional heating with neutrals and electrons (K/s)
    arma_cube heating_neutral_heat_transfer_scgc;
    arma_cube heating_electron_heat_transfer_scgc;

    /// Total heating sources
    arma_cube heating_sources_total;

    /// Specific heat (constant volume):
    arma_cube Cv_scgc;

    /// Heat Conduction:
    arma_cube lambda;

  };

  // bulk quantities (states):
  /// Bulk Density, which is the electron density (/m3):
  arma_cube density_scgc;

  /// Bulk velocity (m/s):
  std::vector<arma_cube> velocity_vcgc;

  /// Bulk Temperature (K):
  arma_cube temperature_scgc;

  /// Bulk conduction term for the ions
  arma_cube conduction_scgc;

  /// Electron temperature
  arma_cube electron_temperature_scgc;

  /// Bulk mass density of the ions
  arma_cube rho_scgc;

  /// Mean major mass of the ions:
  arma_cube mean_major_mass_scgc;

  /// cMax is the sound speed + speed in each direction:
  std::vector<arma_cube> cMax_vcgc;

  /// Bulk sound speed of the ions:
  arma_cube sound_scgc;

  /// Bulk gamma of the ions:
  arma_cube gamma_scgc;

  // This is the vector that will contain all of the different species:
  std::vector<species_chars> species;

  // Heating terms:
  /// Bulk collisional heating with neutrals and electrons (K/s)
  arma_cube heating_neutral_friction_scgc;
  arma_cube heating_electron_friction_scgc;

  /// Bulk collisional heating with neutrals and electrons (K/s)
  arma_cube heating_neutral_heat_transfer_scgc;
  arma_cube heating_electron_heat_transfer_scgc;

  /// Total heating sources
  arma_cube heating_sources_total;

  /// Specific heat (constant volume):
  arma_cube Cv_scgc;

  /// Heat Conduction (bulk):
  arma_cube lambda;

  // Electrodynamics:
  /// Electric potential:
  arma_cube potential_scgc;
  /// Electric field:
  std::vector<arma_cube> efield_vcgc;
  /// E x B drift velocity (electron velocity):
  std::vector<arma_cube> exb_vcgc;
  /// Energy flux of diffuse electron aurora (mW/m2, tbc):
  arma_mat eflux;
  /// Average energy of diffuse electron aurora (keV, tbc):
  arma_mat avee;

  // Some variables that we are going to use in calc_ion_v:
  std::vector<arma_cube> gravity_vcgc;
  std::vector<arma_cube> wind_acc;
  std::vector<arma_cube> total_acc;
  std::vector<arma_cube> efield_acc;
  std::vector<arma_cube> a_par;
  std::vector<arma_cube> a_perp;
  std::vector<arma_cube> a_x_b;
  std::vector<arma_cube> grad_Pi_plus_Pe;
  arma_cube rho, nuin, nuin_sum, Nie, sum_rho;
  arma_cube top, bottom;


  /// Number of species to advect:
  int nSpeciesAdvect;

  /// IDs of species to advect:
  std::vector<int> species_to_advect;

  // names and units
  const std::string density_name = "Neutral Bulk Density";
  const std::string density_unit = "/m3";

  std::vector<std::string> velocity_name;
  std::vector<std::string> par_velocity_name;
  std::vector<std::string> perp_velocity_name;
  const std::string velocity_unit = "m/s";

  const std::string temperature_name = "Temperature";
  const std::string temperature_unit = "K";

  const std::string potential_name = "Potential";
  const std::string potential_unit = "Volts";

  // --------------------------------------------------------------------
  // Functions:

  /**********************************************************************
     \brief Initialize the ions
     \param grid The grid to define the ions on
     \param planet contains information about the species to simulate
   **/
  Ions(Grid &grid, Planets planet);

  /**********************************************************************
     \brief Creates the variables within the species_chars structure
     \param grid The grid to define the ions on
   **/
  species_chars create_species(Grid &grid);

  /**********************************************************************
     \brief
     \param planet contains information about the species to simulate
   **/
  int read_planet_file(Planets planet);

  /**********************************************************************
     \brief Initialize the ion temperature (to the neutral temperature)
     \param neutrals the neutral class to grab the temperature from
     \param grid The grid that the ions are defined on
   **/
  void init_ion_temperature(Neutrals neutrals, Grid &grid);

  /**********************************************************************
     \brief Sets the floor of the ion densities, just in case!
     \param none
   **/
  void set_floor();

  /**********************************************************************
     \brief Sum the ion component densities to calculate the electron density
     \param none
   **/
  void fill_electrons();

  /**********************************************************************
     \brief Calculates the sound speed of the ions
     \param none
   **/
  void calc_sound_speed();

  /**********************************************************************
     \brief Calculates cMax of the ions
     \param none
   **/
  void calc_cMax();

  /**********************************************************************
     \brief Calculate the individual and bulk specific heats
   **/
  void calc_specific_heat();

  /**********************************************************************
     \brief Calculate the individual and bulk thermal conductivities
   **/
  void calc_lambda();

  /**********************************************************************
     \brief Sets the boundary conditions of the ions
     \param grid The grid that the ions are defined on
     \param time The time class to get dt and the current time
     \param indices The indices class to get different indices that may be needed
   **/
  bool set_bcs(Grid &grid, Times time, Indices indices);

  /**********************************************************************
     \brief Sets the upper boundary conditions for the ions
     \param grid The grid that the ions are defined on
   **/
  bool set_upper_bcs(Grid &grid);

  /**********************************************************************
     \brief Sets the lower boundary condition for the ions
     \param grid The grid that the ions are defined on
     \param time The time class to get dt and the current time
     \param indices The indices class to get different indices that may be needed
   **/
  bool set_lower_bcs(Grid &grid, Times time, Indices indices);

  /**********************************************************************
     \brief Advect the ions along the 3rd dimension (could be altitude)
     \param grid The grid that the ions are defined on
     \param time The time class to get dt and the current time
   **/
  bool advect_vertical(Grid &grid, Times time);

  /**********************************************************************
     \brief Get the ID of the ion species with the given name
     \param name a string that describes the species
   **/
  int get_species_id(const std::string &name)const;

  /**********************************************************************
     \brief Calculates the electric field
     \param grid The grid that the ions are defined on
   **/
  void calc_efield(Grid &grid);

  /**********************************************************************
     \brief Calculates the E x B drift
     \param grid The grid that the ions are defined on
   **/
  void calc_exb_drift(Grid &grid);

  /**********************************************************************
     \brief Calculate the ion drift
     \param neutrals these are needed for the collision terms
     \param grid The grid that the ions are defined on
     \param dt the delta-t for the current time
   **/
  void calc_ion_drift(Neutrals &neutrals,
		      Grid &grid,
		      precision_t dt);
  
  /**********************************************************************
     \brief Calculate the ion + electron pressure gradient
     \param iIon which ion to act upon
     \param grid this is the grid to solve the equation on
   **/
  std::vector<arma_cube> calc_ion_electron_pressure_gradient(int64_t iIon,
                                                             Grid grid);

  /**********************************************************************
     \brief Calculates the ion temperature(s) on the given grid
     \param neutrals these are needed for the collision terms
     \param grid this is the grid to solve the equation on
     \param time the time class to know dt
   **/
  void calc_ion_temperature(const Neutrals &neutrals, Grid &grid, Times time);

  /**********************************************************************
     \brief Calculates the electron temperature on the given grid
     \param neutrals these are needed for the collision terms
     \param grid this is the grid to solve the equation on
     \param time the time class to know dt
   **/
  void calc_electron_temperature(Neutrals neutrals, Grid &grid, Times time);

  /**********************************************************************
  /// @brief Calculate epsilon
  /// @details intermediate variable used in photoelectron & ionization heating
  /// From (Smithro & Solomon, 2008).
  /// @param neutrals
  /// @return epsilon
  **/
  arma_cube calc_epsilon(Neutrals &neutrals);

  /**********************************************************************
    \brief Calculates photoelectron heating
    \details Based on (Swartz & Nisbet, 1972) & (Smithro & Solomon, 2008)
            Uses equations 9-12 from (Zhu & Ridley, 2016)
            https://doi.org/10.1016/j.jastp.2016.01.005
    \param epsilon
    \return Qphe
  **/
  arma_cube calc_photoelectron_heating(arma_cube epsilon);

  /**********************************************************************
    \brief Calculates auroral heating
    \details NOTE: in GITM this is solved separately for ion precipitation & auroral
        ionization. In Aether these are both in ions.species[iIon].ionization_scgc...
    \param epsilon
    \return Qaurora
  **/
  arma_cube calc_ionization_heating(arma_cube epsilon);

  /**********************************************************************
    \brief Calculates electron-ion (elastic) collisional heating
    \details From Schunk and Nagy 2009, and Bei-Chen Zhang and Y. Kamide 2003
    - This differs slightly from the GITM implementation, which assumes several ion species are present.
      Instead, here we use each ion species for the sum.
    - electon-ion collision frequency (from Schunk and Nagy 2009) = 5.45E-5
    - This is capable of handling BOTH the bulk & individual ion temperatures
    \return vector<Qeicp, Qeicm, Qeic_v>
  **/
  std::vector<arma_cube> calc_electron_ion_collisions();

  /**********************************************************************
    \brief Calculates electron-neutral elastic collisional heating
    \details From Schunk and Nagy 2009
    \param neutrals
    \return vector<Qencp, Qencm, Qenc_v>
  **/
  std::vector<arma_cube> calc_electron_neutral_elastic_collisions(
    Neutrals &neutrals);

  /**********************************************************************
    \brief Calculates the electron-neutral inelastic collisional heating
    \details From Schunk and Nagy 2009 pages 277, 282.
    This includes N2, O2 rotation, fine structure, O(1D) exitation & vibration, N2 vibration.
    See equation 15 from (Zhu, Ridley, Deng, 2016) https://doi.org/10.1016/j.jastp.2016.01.005
    \param neutrals
    \return vector<Qencp, Qencm, Qenc_v>
  **/
  std::vector<arma_cube> calc_electron_neutral_inelastic_collisions(
    Neutrals &neutrals);

  /**********************************************************************
    \brief Calculate the thermoelectric current (same at all altitudes)
    \details Use eq. 6 of https://doi.org/10.1016/j.jastp.2016.01.005
    - Since we do not know e- parallel velocity, the dipole needs to do it this way too.
    \param grid
    \return arma_mat JParaAlt
  **/
  arma_mat calc_thermoelectric_current(Grid &grid);

  /**********************************************************************
     \brief Check all of the variables for nonfinites, such as nans
     \param none
   **/
  bool check_for_nonfinites(std::string location);

  /**********************************************************************
     \brief Run through a test of an arma_cube to see if it contains nans
     \param variable this is the variable to check
   **/
  void nan_test(std::string variable);

  /**********************************************************************
     \brief Read or write to the restart file
     \param dir the directory to read or write from/to
     \param DoRead whether to read (true) or write (false)
   **/
  bool restart_file(std::string dir, std::string cGridtype, bool DoRead);

  /**********************************************************************
     \brief Exchange messages between processors
     \param grid The grid to define the ions on
   **/
  // bool exchange(Grid &grid);
  bool exchange_old(Grid &grid);

  /**********************************************************************
     \brief Vertical advection solver - Rusanov
     \param grid The grid to define the neutrals on
     \param time contains information about the current time
   **/
  void solver_vertical_rusanov(Grid &grid, Times time);

};
#endif  // INCLUDE_IONS_H_
