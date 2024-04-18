// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef AETHER_INCLUDE_PLANETS_H_
#define AETHER_INCLUDE_PLANETS_H_

/**************************************************************
 * A class for keeping track of all of the planetary characteristics
 *
 *
 **************************************************************/

#include <string>
#include <vector>

class Planets {

// -----------------------------------------------------------------------
// Public functions and variables
// -----------------------------------------------------------------------

public:

  // --------------------------------------------------------------------
  // Functions:

  /**********************************************************************
     \brief Initialize the Planet class
   **/
  Planets();

  /**********************************************************************
     \brief Returns the distance from the star to the planet being modeled
     \param time needed for getting the current time in different formats
   **/
  precision_t get_star_to_planet_dist(Times time);

  /**********************************************************************
     \brief Returns the orbit angle of the planet around the star
     \param time needed for getting the current time in different formats
   **/
  precision_t get_orbit_angle(Times time);

  /**********************************************************************
     \brief Returns the declination angle of the planet (zero = equinox)
     \param time needed for getting the current time in different formats
     \param
   **/
  precision_t get_declination(Times time);

  /**********************************************************************
     \brief Gets mu, which is planet mass times gravitational constant
   **/
  precision_t get_mu();

  /**********************************************************************
     \brief Gets J2
   **/
  precision_t get_J2();

  /**********************************************************************
     \brief Returns radius of the planet, which can be a function of latitude

     currently, this ignores the latitude, but should be implemented.

     \param latitude the latitude to get the radius at.
   **/
  precision_t get_radius(precision_t latitude);

  /**********************************************************************
     \brief Returns the longitude offset to convert from longitude to local time

     this function returns a value that when added to the longitude (in
     radians) will provide the local time (in radians)

     \param time needed for getting the current time in different formats
   **/
  precision_t get_longitude_offset(Times time);

  /**********************************************************************
     \brief Returns the sin of the planet's declination angle
     \param time needed for getting the current time in different formats
   **/
  precision_t get_sin_dec(Times time);

  /**********************************************************************
     \brief Returns the cos of the planet's declination angle
     \param time needed for getting the current time in different formats
   **/
  precision_t get_cos_dec(Times time);

  /**********************************************************************
     \brief Returns the location of the center of the dipole (in meters)
   **/
  std::vector<float> get_dipole_center();

  /**********************************************************************
     \brief Returns the rotation angle of the dipole in longitude (radians)
   **/
  precision_t get_dipole_rotation();

  /**********************************************************************
     \brief Returns the tilt angle of the dipole (co-latitude) (radians)
   **/
  precision_t get_dipole_tilt();

  /**********************************************************************
     \brief Returns the strength of the dipole at the surface (in nT)
   **/
  precision_t get_dipole_strength();

  /**********************************************************************
     \brief Returns omega (rotation rate) of the planet 
   **/
  precision_t get_omega();

  /**********************************************************************
     \brief returns neutrals json for neutral density BCs
   **/
  json get_neutrals();
  
  /**********************************************************************
     \brief returns neutral temperature json for temperature ICs
   **/
  json get_temperatures();
  
  /**********************************************************************
     \brief returns ions json for ion density characteristics
   **/
  json get_ions();
  
  /**********************************************************************
     \brief Check to see if internal state of class is ok
   **/
  
  bool is_ok();

// -----------------------------------------------------------------------
// Private functions and variables
// -----------------------------------------------------------------------
  
 private:

  /// A structure to describe the planetary characteristics for each planet
  /// Many of these are constants and read in from a file, and
  /// Some are updated (for the planet being used) during the run
  struct planet_chars {

    // ---------------------------------------
    // These are set by an input file:

    /// Planet's name:
    std::string name;

    /// Semi Major Axis of the planet around the star
    precision_t semimajoraxis;
    /// Eccentricity of the planet around the star
    precision_t eccentricity;
    /// Inclination of the planet around the star
    precision_t inclination;
    /// Mean longitude of the planet around the star
    precision_t meanlongitude;
    /// Perihelion Longitude of the planet around the star
    precision_t perihelionlongitude;
    /// Node Longitude of the planet around the star
    precision_t nodelongitude;

    /// Rate of Change of the Semi Major Axis of the planet around the star
    precision_t rates_semimajoraxis;
    /// Rate of Change of the Eccentricity of the planet around the star
    precision_t rates_eccentricity;
    /// Rate of Change of the Inclination of the planet around the star
    precision_t rates_inclination;
    /// Rate of Change of the Mean longitude of the planet around the star
    precision_t rates_meanlongitude;
    /// Rate of Change of the Perihelion Longitude of the planet around the star
    precision_t rates_perihelionlongitude;
    /// Rate of Change of the Node Longitude of the planet around the star
    precision_t rates_nodelongitude;

    /// Tilt of the planet with respect to the ecliptic plane
    precision_t planet_tilt;

    /// rotation_period, omega, and length of day are all related to each other
    /// Rotation period of the planet (time to rotate 360 deg)
    precision_t rotation_period;
    /// Rotation rate (1/period)
    precision_t omega;
    /// Length of a day - time to rotate from noon to noon
    precision_t length_of_day;

    /// Length of a year for the planet
    double length_of_year;
    /// Longitude of midnight at January 1, 2000 (to convert to local time)
    precision_t longitude_jb2000;

    /// planet mass
    precision_t mass;
    /// mu = planet mass * gravitational constant
    precision_t mu;
    /// Mean radius of the equator
    precision_t equator_radius;
    /// Radius at pole
    precision_t polar_radius;
    /// mean radius
    precision_t radius;
    /// Difference between polar and equator radius
    precision_t delta_radius;

    /// Dipole strength (nT)
    precision_t dipole_strength;
    /// Dipole rotation around in longitude (radians)
    precision_t dipole_rotation;
    /// Dipole tilt from the rotation axis (co-latitude of pole) (radians)
    precision_t dipole_tilt;
    /// Offset of the dipole center from the geographic center of planet
    std::vector<float> dipole_center{ 0.0, 0.0, 0.0 };

    //J2 value
    precision_t J2;

    // ---------------------------------------
    // These are updated during the run:

    /// The declination, sin declination, and cos declination of the planet
    precision_t declination, sin_dec, cos_dec;
    /// Angle to add to longitude to give local time
    precision_t longitude_offset;

    /// Distance between the star and planet
    precision_t star_planet_distance;
    /// Angle of the planet around the star
    precision_t orbit_angle;
    /// L sub s (Fractional day of the planet's year):
    precision_t ls;

    /// How often to update these variables
    double update_time;

  };

  /// All of the planets in the file are stored here
  std::vector<planet_chars> planets;

  /// The chose planet is stored here
  planet_chars planet;

  /// Information about the neutrals of the planet
  json neutrals;

  /// Information about the initial temperature of the planet
  json temperatures;
  
  /// Information about the ions of the planet
  json ions;

  /// An internal variable to hold the state of the class
  bool IsOk;

  // --------------------------------------------------------------------
  // Functions:

  /**********************************************************************
     \brief Sets the planetary characterists for the chosen planet

     This function just copies over the chosen planet's characteristics
     from the vector into the stand-alone structure.
   **/
  bool set_planet();

  /**********************************************************************
     \brief Updates the planetary characteristics that depend on time

     time is needed in a variety of formats, since many of the
     calculations require different time formats.

     \param time needed for getting the current time in different formats
   **/
  int update(Times time);

  /**********************************************************************
     \brief Reads in the planetary characteristics and stores them
   **/
  bool read_file();

  /**********************************************************************
     \brief Reads in the planets neutral, ion, temperature characteristics
   **/
  bool read_planet_specific_file();

};

#endif  // INCLUDE_PLANETS_H_
