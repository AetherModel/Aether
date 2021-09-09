// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef AETHER_INCLUDE_PLANETS_H_
#define AETHER_INCLUDE_PLANETS_H_

#include <string>
#include <vector>

class Planets {

 public:
  Planets(Inputs args, Report report);
  precision_t get_star_to_planet_dist(Times time);
  precision_t get_orbit_angle(Times time);
  precision_t get_declination(Times time);
  precision_t get_mu();
  precision_t get_radius(precision_t latitude);
  precision_t get_longitude_offset(Times time);
  precision_t get_sin_dec(Times time);
  precision_t get_cos_dec(Times time);

  std::vector<float> get_dipole_center();
  precision_t get_dipole_rotation();
  precision_t get_dipole_tilt();
  precision_t get_dipole_strength();
  
 private:

  int set_planet(Inputs args, Report report);
  int update(Times time);

  struct planet_chars {

    // These are set by an input file:

    std::string name;

    precision_t semimajoraxis;
    precision_t eccentricity;
    precision_t inclination;
    precision_t meanlongitude;
    precision_t perihelionlongitude;
    precision_t nodelongitude;

    precision_t rates_semimajoraxis;
    precision_t rates_eccentricity;
    precision_t rates_inclination;
    precision_t rates_meanlongitude;
    precision_t rates_perihelionlongitude;
    precision_t rates_nodelongitude;

    precision_t planet_tilt;
    // rotation_period, omega, and length of day are all related to each other
    precision_t rotation_period;
    precision_t omega;
    precision_t length_of_day;
    double length_of_year;
    precision_t longitude_jb2000;
    precision_t longitude_offset;

    precision_t mass;
    precision_t mu;
    precision_t equator_radius;
    precision_t polar_radius;
    precision_t radius;

    precision_t dipole_strength;
    precision_t dipole_rotation;
    precision_t dipole_tilt;
    std::vector<float> dipole_center{ 0.0, 0.0, 0.0 };

    // These are updated during the run:

    precision_t declination, sin_dec, cos_dec;

    precision_t star_planet_distance;
    precision_t orbit_angle;
    precision_t ls;

    double update_time;
  };

  std::vector<planet_chars> planets;

  planet_chars planet;

  int read_file(Inputs args, Report report);
};

#endif  // INCLUDE_PLANETS_H_
