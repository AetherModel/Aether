
#include <math.h>
#include <vector>
#include <string>
#include <iostream>
#include <sstream> // std::stringstream
#include <fstream>

#include "../include/constants.h"
#include "../include/planets.h"
#include "../include/times.h"
#include "../include/report.h"
#include "../include/file_input.h"

// -----------------------------------------------------------------------------
// Constructor (initiaze the class):
// -----------------------------------------------------------------------------

Planets::Planets(Inputs args, Report report) {
  int iErr = 0;

  iErr = read_file(args, report);
  if (iErr == 0) iErr = set_planet(args, report);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------

float Planets::get_longitude_offset(Times time) {
  int iErr = update(time);
  return planet.longitude_offset;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------

float Planets::get_sin_dec(Times time) {
  int iErr = update(time);
  return planet.sin_dec;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------

float Planets::get_cos_dec(Times time) {
  int iErr = update(time);
  return planet.cos_dec;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------

float Planets::get_radius(float latitude) {
  // Should modify this to allow an oblate spheriod, but not now.
  return planet.radius;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------

float Planets::get_mu() {
  return planet.mu;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------

float Planets::get_star_to_planet_dist(Times time) {

  int iErr = update(time);
  return planet.star_planet_distance;

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------

float Planets::get_orbit_angle(Times time) {

  int iErr = update(time);
  return planet.orbit_angle;

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------

float Planets::get_declination(Times time) {

  int iErr = update(time);
  return planet.declination;;

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------

int Planets::update(Times time) {

  int iErr = 0;
  
  // Update planetary stuff once a minute:
  if (time.get_current() - planet.update_time > 60.0) {

    planet.update_time = time.get_current();

    float sma = planet.semimajoraxis +
      planet.rates_semimajoraxis * time.get_orbittime();
    float ecc = planet.eccentricity +
      planet.rates_eccentricity * time.get_orbittime();
    float inc = planet.inclination +
      planet.rates_inclination * time.get_orbittime();
    float meanlon = planet.meanlongitude +
      planet.rates_meanlongitude * time.get_orbittime();
    float node_long = planet.nodelongitude +
      planet.rates_nodelongitude * time.get_orbittime();
    float perilon = planet.perihelionlongitude +
      planet.rates_perihelionlongitude * time.get_orbittime();

    // Compute argument of perihelion and mean anomaly
    float arg_peri = perilon - node_long;

    // computation of M for Jupiter and out is supposed to be modified
    // by additional terms for the time interval 3000BC -
    // 30000AD... This probably doesn't matter.
    float meananomaly = fmod(meanlon - perilon + 360.0, 360.0);
    if (meananomaly > 180.) meananomaly = meananomaly - 360.0;

    // Need to solve Kepler's equation by iterating
    float d_ecc_anomaly = 10000.0;
    float tol = 1.0e-6;
    float ecc_deg = ecc * rtod;
    float ecc_anomaly = meananomaly+ecc_deg*sin(meananomaly*dtor);
    int i = 0;
    float dm = 0.0;
    while (abs(d_ecc_anomaly) > tol && i < 100) {
      dm = meananomaly - (ecc_anomaly - ecc_deg*sin(ecc_anomaly*dtor));
      d_ecc_anomaly = dm/(1-ecc*cos(ecc_anomaly*dtor));
      ecc_anomaly = ecc_anomaly+d_ecc_anomaly;
      i++;
    }

    // Get heliocentric coordinates, TrueAnomaly and sunplanetdistance
    float x_heliocentric = sma*(cos(ecc_anomaly*dtor)-ecc);
    float y_heliocentric = sma*sqrt(1-ecc*ecc)*sin(ecc_anomaly*dtor);
    float z_heliocentric = 0.0;

    float true_anomaly = atan2(y_heliocentric,x_heliocentric)*rtod;

    planet.star_planet_distance = sqrt(x_heliocentric * x_heliocentric +
				       y_heliocentric * y_heliocentric);

    // convert to J2000 coordinates with x-axis aligned with vernal equinox so
    // we can get solar longitude in the coorect system.  We don't need z.
    float x_ecl =
      x_heliocentric *
      (cos(arg_peri*dtor)*cos(node_long*dtor) -
       sin(arg_peri*dtor)*sin(node_long*dtor)*cos(inc*dtor)) +
      y_heliocentric *
      (-sin(arg_peri*dtor)*cos(node_long*dtor) -
       cos(arg_peri*dtor)*sin(node_long*dtor)*cos(inc*dtor));

    float y_ecl =
      x_heliocentric *
      (cos(arg_peri*dtor)*sin(node_long*dtor) +
       sin(arg_peri*dtor)*cos(node_long*dtor)*cos(inc*dtor)) +
      y_heliocentric *
      (-sin(arg_peri*dtor)*sin(node_long*dtor) +
       cos(arg_peri*dtor)*cos(node_long*dtor)*cos(inc*dtor));

    // Calculate orbit angle, aka Ls. In this CS, need the angle from -x axis.
    planet.orbit_angle = atan(y_ecl/x_ecl);
    if (x_ecl > 0) planet.orbit_angle = planet.orbit_angle+pi;
    if (x_ecl < 0 && y_ecl > 0) planet.orbit_angle = planet.orbit_angle + 2*pi;

    planet.declination = atan(tan(planet.planet_tilt*dtor)*sin(planet.orbit_angle));

    planet.sin_dec = sin(planet.declination);
    planet.cos_dec = cos(planet.declination);

    // This is take from :
    //   https://www.aa.quae.nl/en/reken/zonpositie.html (section 7)
    //   Basically, it identifies the longitude that is pointed along
    //   the sun-planet line at 00 UT on Jan 1 2000.  The rotation
    //   that is used on the webpage is to calculate the sidereal
    //   time, though, which we want the local time.  So, we use the
    //   length of a day instead of the rotation rate.  Ok, this is
    //   idea is correct, but the angles reported on this website at
    //   not correct.  They may be providing the angle for RAAN or
    //   something like that. It is not defining the longitude for
    //   midnight or noon. Ugh. (Fixed Earth, but no other planet.)
    double nEarthDaysSince2000 = (time.get_julian_day() - j2000);
    double deg_per_day = 360.0 / (planet.length_of_day/seconds_per_day);
    double deg_since_2000 = deg_per_day * nEarthDaysSince2000;
    double rotations = deg_since_2000/360.0;
    double left_over = (rotations - int(rotations)) * 360.0;
    // put into radians, so it is consistent with the rest of the code:
    planet.longitude_offset =
      fmod(planet.longitude_jb2000 + left_over + 180.0,360.0) * dtor;

  }

  return iErr;

}

int Planets::set_planet(Inputs args, Report report) {

  int iErr = 0;
  int IsFound = 0;

  for (int i=0; i<planets.size(); i++) {
    if (planets[i].name == args.get_planet()) {
      IsFound = 1;
      planet.name = planets[i].name;
      if (report.test_verbose(2)) {
	std::cout << "Planet set to : " << planet.name << "\n";
      }
      planet.semimajoraxis = planets[i].semimajoraxis;
      planet.eccentricity = planets[i].eccentricity;
      planet.inclination = planets[i].inclination;
      planet.meanlongitude = planets[i].meanlongitude;
      planet.perihelionlongitude = planets[i].perihelionlongitude;
      planet.nodelongitude = planets[i].nodelongitude;
      planet.rates_semimajoraxis = planets[i].rates_semimajoraxis;
      planet.rates_eccentricity = planets[i].rates_eccentricity;
      planet.rates_inclination = planets[i].rates_inclination;
      planet.rates_meanlongitude = planets[i].rates_meanlongitude;
      planet.rates_perihelionlongitude = planets[i].rates_perihelionlongitude;
      planet.rates_nodelongitude = planets[i].rates_nodelongitude;
      planet.planet_tilt = planets[i].planet_tilt;

      double mu = solar_mass * gravitational_constant;
      double a = planet.semimajoraxis * au_to_m;
      double loy = 2 * pi * sqrt(a*a*a/mu);
      if (abs(loy - planets[i].length_of_year)/loy > 0.01) {
	std::cout << "Hmmm.... For some reason, the calculated length of year is different\n";
	std::cout << "Read in : " << planet.length_of_year/seconds_per_day << " (days) \n";
	std::cout << "Calculated from SMA : " << loy/seconds_per_day << " (days)\n";
	std::cout << "Trusting SMA version!\n";
	iErr = 1;
      }
      planet.length_of_year = loy;
      planet.length_of_day = planets[i].length_of_day;
      planet.longitude_jb2000 = planets[i].longitude_jb2000;

      float rotrate = 2*pi/planet.length_of_day*(1.0+1.0/(loy/seconds_per_day));
      planet.omega = rotrate; // frequency (rad/s)
      planet.rotation_period = 2*pi/rotrate; // (seconds)

      planet.mass = planets[i].mass;
      planet.mu = planets[i].mass * gravitational_constant;
      planet.equator_radius = planets[i].equator_radius * 1000.0; // km -> m
      planet.polar_radius = planets[i].polar_radius * 1000.0; // km -> m
      // Looking at Earth and Saturn, it seems like the Volumetric mean radius
      // is at roughly 47 deg (cos(47)=0.68)
      // Obviously an approximation...
      planet.radius = 0.68 * planet.equator_radius + 0.32 * planet.polar_radius;
      if (report.test_verbose(2))
	std::cout << "Planet Radius set to : "
		  << planet.radius/1000.0 << " (km)\n";

      planet.update_time = -1e32;

      break;
    }
  }

  if (!IsFound) {
    std::cout << "Can't file planet " << args.get_planet()
	      << " in planet file information!\n";
    iErr = 1;
  }

  return iErr;

}

int Planets::read_file(Inputs args, Report report) {

  planet_chars tmp;
  std::string line, col;
  std::ifstream myFile;
  int iErr = 0;

  report.print(1, "Reading planetary file : " + args.get_planetary_file());

  myFile.open(args.get_planetary_file());

  if (!myFile.is_open()) {
    std::cout << "Could not open planetary file : "
	      << args.get_planetary_file() << "\n";
    iErr = 1;
  } else {

    if (myFile.good()) {

      // Header line can be ignored:
      getline(myFile,line);
      getline(myFile,line);

      while (getline(myFile,line)) {

	report.print(5, line);
	std::stringstream ss(line);

	getline(ss, tmp.name, ',');
	tmp.name = make_lower(tmp.name);
	getline(ss, col, ',');
	if (col.size() > 1) {
	  tmp.semimajoraxis = stof(col);
	  getline(ss, col, ',');
	  tmp.eccentricity = stof(col);
	  getline(ss, col, ',');
	  tmp.inclination = stof(col);
	  getline(ss, col, ',');
	  tmp.meanlongitude = stof(col);
	  getline(ss, col, ',');
	  tmp.perihelionlongitude = stof(col);
	  getline(ss, col, ',');
	  tmp.nodelongitude = stof(col);
	  getline(ss, col, ',');
	  tmp.rates_semimajoraxis = stof(col);
	  getline(ss, col, ',');
	  tmp.rates_eccentricity = stof(col);
	  getline(ss, col, ',');
	  tmp.rates_inclination = stof(col);
	  getline(ss, col, ',');
	  tmp.rates_meanlongitude = stof(col);
	  getline(ss, col, ',');
	  tmp.rates_perihelionlongitude = stof(col);
	  getline(ss, col, ',');
	  tmp.rates_nodelongitude = stof(col);
	  // Length of Day (hours)
	  getline(ss, col, ',');
	  tmp.length_of_day = stof(col)*seconds_per_hour;
	  // Length of year (days)
	  getline(ss, col, ',');
	  tmp.length_of_year = stof(col)*seconds_per_day;
	  // Length of year (days)
	  getline(ss, col, ',');
	  tmp.longitude_jb2000 = stof(col);
	  // Mass of planet
	  getline(ss, col, ',');
	  tmp.mass = stof(col);
	  // Equatorial Radius
	  getline(ss, col, ',');
	  tmp.equator_radius = stof(col);
	  // Polar Radius
	  getline(ss, col, ',');
	  tmp.polar_radius = stof(col);
	  // Planetary Tilt Angle
	  getline(ss, col, ',');
	  tmp.planet_tilt = stof(col);

	  planets.push_back(tmp);

	}

      }

    }

    myFile.close();

  }

  return iErr;

}
