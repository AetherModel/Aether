// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_NEUTRALS_H_
#define INCLUDE_NEUTRALS_H_

/**************************************************************
 * \class Neutrals
 *
 * \brief Defines the neutral states
 * 
 * The Neutrals class defines the neutrals states as well as a bunch
 * of derived states and source/loss terms.  The initial temperature
 * structure as well as the lower boundary densities can be set
 * through the planet input file.
 *
 * \author Aaron Ridley
 *
 * \date 2021/03/28 
 *
 **************************************************************/

#include <string>
#include <vector>

#include "aether.h"

class Neutrals {

 public:

  /// This struct contains all of the information needed for a single
  /// species of neutrals.  We will then have a vector of these species.

  struct species_chars {

    /// Name of the species
    std::string cName;

    /// Mass of the species (kg)
    float mass;

    /// Vibrations of species (for calculation specific heat)
    float vibe;

    /// Advect this species? (1 = yes, 0 = no)
    int DoAdvect;

    /// Number density of species (/m3)
    fcube density_scgc;

    /// Diffusion through other neutral species:
    std::vector<float> diff0;
    std::vector<float> diff_exp;

    /// Neutral - Ion collision frequency coefficients
    std::vector<float> neutral_ion;

    /// Thermal conduction coefficients:
    float thermal_cond;
    float thermal_exp;

    /// Which row in the EUV CSV file is for absorption:
    int iEuvAbsId_;
    /// How many rows in the EUV CSV file are for ionization of this species?
    int nEuvIonSpecies;
    /// Which row in the EUV CSV file if for the particular ionization?
    std::vector<int> iEuvIonId_;
    /// Which ion species results from the ionization?
    std::vector<int> iEuvIonSpecies_;

    // --------------------------------------------------
    // Some derived quantities:

    /// Chapman Integrals for the species for EUV calculation (/m2)
    fcube chapman_scgc;

    /// Scale height for the species (m)
    fcube scale_height_scgc;

    // --------------------------------------------------
    // Sources and Losses:

    /// How much of this species is lost to ionization (/m3/s)
    fcube ionization_scgc;

    /// Chemistry source rate (/m3/s)
    fcube sources_scgc;
    
    /// Chemistry loss rate (/m3/s)
    fcube losses_scgc;

    /// If we want a fixed lower BC:
    float lower_bc_density;
  };

  /// bulk number density (/m3)
  fcube density_scgc;

  /// bunk temperature (K)
  fcube temperature_scgc;

  /// bulk mass density (kg/m3)
  fcube rho_scgc;

  /// mean major mass (kg)
  fcube mean_major_mass_scgc;

  /// mean pressure (Pa)
  fcube pressure_scgc;

  /// speed of sound (m/s)
  fcube sound_scgc;

  /// Specific heat (constant volume):
  fcube Cv_scgc;

  /// Bulk Gamma:
  fcube gamma_scgc;

  /// Bulk thermal heat conduction:
  fcube kappa_scgc;

  /// Vector of all species-specific items:
  std::vector<species_chars> species;

  /// Maximum Chapman integral (will give nearly infinite tau in EUV)
  float max_chapman = 1.0e26;

  // Source terms:

  /// Bulk neutral thermal conduction temperature change rate (K/s)
  fcube conduction_scgc;

  /// Bulk neutral EUV heating temperatuare change (K/s)
  fcube heating_euv_scgc;

  /// Nuetral gas direct absorption heating efficiency (~5%)
  float heating_efficiency;

  /// Initial temperature profile, read in through the planet.in file:
  float *initial_temperatures, *initial_altitudes;
  int nInitial_temps = 0;

  // names and units
  std::string density_name = "Neutral Bulk Density";
  std::string density_unit = "(/m3)";

  std::vector<std::string> velocity_name;
  std::string velocity_unit = "(m/s)";

  std::string temperature_name = "Temperature";
  std::string temperature_unit = "(K)";

  // --------------------------------------------------------------------
  // Functions:

  /**********************************************************************
     \brief Initialize the neutrals
     \param grid The grid to define the neutrals on
     \param input info about how user has configured things
     \param report allow reporting to occur
   **/
  Neutrals(Grid grid, Inputs input, Report report);

  /**********************************************************************
     \brief Creates the variables within the species_chars structure
     \param grid The grid to define the neutrals on
     \param input info about how user has configured things
     \param report allow reporting to occur
   **/
  species_chars create_species(Grid grid);

  /**********************************************************************
     \brief Read in the planet-specific file

     This file specifies the species to model, their masses, 
     diffusion coefficients and all of the other things needed
     for specifying the neutrals.

     \param input info about how user has configured things
     \param report allow reporting to occur
   **/
  int read_planet_file(Inputs input, Report report);

  /**********************************************************************
     \brief Sets the initial conditions of the neutrals
     \param grid The grid to define the neutrals on
     \param input info about how user has configured things
     \param report allow reporting to occur
   **/
  int initial_conditions(Grid grid, Inputs input, Report report);

  /**********************************************************************
     \brief temporary function to set neutral densities with in the model

     This function integrates the species densities from the bottom
     of the model using a hydrostatic approximation and the bulk
     temperature.  It is temporary until we get a vertical solver.

     \param grid The grid to define the neutrals on
     \param report allow reporting to occur
   **/
  void fill_with_hydrostatic(Grid grid, Report report);

  /**********************************************************************
     \brief Calculate the bulk mass density from individual species densities
     \param report allow reporting to occur
   **/
  void calc_mass_density(Report &report);

  /**********************************************************************
     \brief Calculate the bulk specific heat from individual species
     \param report allow reporting to occur
   **/
  void calc_specific_heat(Report &report);

  /**********************************************************************
     \brief Calculate the chapman integrals for the individual species
     \param grid The grid to define the neutrals on
     \param report allow reporting to occur
   **/
  void calc_chapman(Grid grid, Report &report);

  /**********************************************************************
     \brief Calculate the neutral bulk vertical thermal conduction
     \param grid The grid to define the neutrals on
     \param time The times within the model (dt is needed)
     \param report allow reporting to occur
   **/
  void calc_conduction(Grid grid, Times time, Report &report);

  /**********************************************************************
     \brief Add all of the neutral source terms to each of the equations
     \param time The times within the model (dt is needed)
     \param report allow reporting to occur
   **/
  void add_sources(Times time, Report &report);

  /**********************************************************************
     \brief Set boundary conditions for the neutrals
     \param report allow reporting to occur
   **/
  void set_bcs(Report &report);
};

#endif  // INCLUDE_NEUTRALS_H_

