// (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef AETHER_INCLUDE_PLANETS_H_
#define AETHER_INCLUDE_PLANETS_H_

#include <string>
#include <vector>

#include "inputs.h"
#include "report.h"

class Planets {

 public:
  Planets(Inputs args, Report report);
  float get_star_to_planet_dist(Times time);
  float get_orbit_angle(Times time);
  float get_declination(Times time);
  float get_mu();
  float get_radius(float latitude);
  float get_longitude_offset(Times time);
  float get_sin_dec(Times time);
  float get_cos_dec(Times time);

 private:

  int set_planet(Inputs args, Report report);
  int update(Times time);

  struct planet_chars {

    // These are set by an input file:

    std::string name;

    float semimajoraxis;
    float eccentricity;
    float inclination;
    float meanlongitude;
    float perihelionlongitude;
    float nodelongitude;

    float rates_semimajoraxis;
    float rates_eccentricity;
    float rates_inclination;
    float rates_meanlongitude;
    float rates_perihelionlongitude;
    float rates_nodelongitude;

    float planet_tilt;
    float rotation_period;
    float omega; // rotation frequency
    float length_of_day; // this and rotation period are related to each other
    double length_of_year;
    float longitude_jb2000;
    float longitude_offset;

    float mass;
    float mu;
    float equator_radius;
    float polar_radius;
    float radius;

    // These are updated during the run:

    float declination, sin_dec, cos_dec;

    float star_planet_distance;
    float orbit_angle;
    float ls;

    double update_time;

  };

  std::vector<planet_chars> planets;

  planet_chars planet;

  int read_file(Inputs args, Report report);

};

#endif // AETHER_INCLUDE_PLANETS_H_
