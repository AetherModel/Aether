// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "aether.h"

// -----------------------------------------------------------------------------
// This is the head bfield function that determines which bfield
// to get, then gets it.
// -----------------------------------------------------------------------------

bfield_info_type get_bfield(precision_t lon,
                            precision_t lat,
                            precision_t alt,
                            bool DoDebug,
                            Planets planet) {

  std::string function = "get_bfield";
  static int iFunction = -1;

  if (DoDebug)
    report.enter(function, iFunction);

  if (lat > cPI / 2) {
    lat = cTWOPI - lat;
    lon = lon + cPI;

    if (lon > cTWOPI)
      lon = lon - cTWOPI;
  }

  if (lat < -cPI / 2) {
    lat = -cTWOPI + lat;
    lon = lon + cPI;

    if (lon > cTWOPI)
      lon = lon - cTWOPI;
  }

  bfield_info_type bfield_info;

  if (input.get_bfield_type() == "none") {
    bfield_info.b[0] = 0.0;
    bfield_info.b[1] = 0.0;
    bfield_info.b[2] = 0.0;
    bfield_info.lat = lat;
    bfield_info.lon = lon;
  } else if (input.get_bfield_type() == "dipole")
    bfield_info = get_dipole(lon, lat, alt, DoDebug, planet);

  if (DoDebug)
    report.exit(function);

  return bfield_info;
}

// -----------------------------------------------------------------------------
// This function finds the magnetic pole in either the north or south.
// Return values in radians!
// -----------------------------------------------------------------------------

arma_vec get_magnetic_pole(int IsNorth,
                           Planets planet) {

  arma_vec lonlat(2, fill::zeros);

  if (input.get_bfield_type() == "none") {
    // No magnetic field, so set location to the geo pole:
    lonlat(0) = 0.0;

    if (IsNorth)
      lonlat(1) = cPI / 2.0;

    else
      lonlat(1) = -cPI / 2.0;
  } else if (input.get_bfield_type() == "dipole") {
    // This is an approximation right now, due to the fact that the
    // pole on some planets (including Earth), have an offset, so the
    // pole location should be altitude dependent.  For many planets,
    // this is a very small error, but should be fixed at some point.

    precision_t magnetic_pole_rotation = planet.get_dipole_rotation();
    precision_t magnetic_pole_tilt = planet.get_dipole_tilt();

    if (IsNorth) {
      lonlat(0) = magnetic_pole_rotation;
      lonlat(1) = cPI / 2 - magnetic_pole_tilt;
    } else {
      lonlat(0) = cTWOPI - magnetic_pole_rotation;
      lonlat(1) = -cPI / 2 + magnetic_pole_tilt;
    }
  }

  return lonlat;
}
