// (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "../include/inputs.h"
#include "../include/planets.h"
#include "../include/report.h"
#include "../include/bfield.h"
#include "../include/constants.h"

bfield_info_type get_bfield(float lon,
			    float lat,
			    float alt,
			    Planets planet,
			    Inputs input,
			    Report &report) {

  std::string function = "get_bfield";
  report.enter(function);

  if (lat > pi/2) {
    lat = twopi - lat;
    lon = lon + pi;
    if (lon > twopi) lon = lon - twopi;
  }
  if (lat < -pi/2) {
    lat = -twopi + lat;
    lon = lon + pi;
    if (lon > twopi) lon = lon - twopi;
  }

  bfield_info_type bfield_info;
  
  if (input.get_bfield_type() == "none") {
    bfield_info.b[0] = 0.0;
    bfield_info.b[1] = 0.0;
    bfield_info.b[2] = 0.0;
    bfield_info.lat = lat;
    bfield_info.lon = lon;
  } else if (input.get_bfield_type() == "dipole") {
    bfield_info = get_dipole(lon, lat, alt, planet, input, report);
  }
  
  report.exit(function);
  return bfield_info;
  
}

