// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>

#include "../include/aether.h"

void Grid::init_geo_grid(Planets planet, Inputs input, Report &report) {

  std::string function = "Grid::init_geo_grid";
  static int iFunction = -1;
  report.enter(function, iFunction);

  // This is just an example:

  Inputs::grid_input_struct grid_input = input.get_grid_inputs();

  int64_t iLon, iLat, iAlt;

  // Longitudes:
  // - Make a 1d vector
  // - copy it into the 3d cube
  fvec lon1d(nLons);
  float dlon = (grid_input.lon_max - grid_input.lon_min) / (nLons-2*nGCs);
  std::cout << "dlon : " << dlon/dtor << "\n";
  for (iLon=0; iLon < nLons; iLon++)
    lon1d(iLon) = grid_input.lon_min + (iLon-nGCs+0.5) * dlon;

  for (iLat=0; iLat < nLats; iLat++) {
    for (iAlt=0; iAlt < nAlts; iAlt++) {
      geoLon_scgc.subcube(0, iLat, iAlt, nLons-1, iLat, iAlt) = lon1d;
    }
  }

  // Latitudes:
  // - Make a 1d vector
  // - copy it into the 3d cube
  fvec lat1d(nLats);
  float dlat = (grid_input.lat_max - grid_input.lat_min) / (nLats-2*nGCs);
  std::cout << "dlat : " << dlat/dtor << "\n";
  for (iLat=0; iLat < nLats; iLat++)
    lat1d(iLat) = grid_input.lat_min + (iLat-nGCs+0.5) * dlat;

  for (iLon=0; iLon < nLons; iLon++) {
    for (iAlt=0; iAlt < nAlts; iAlt++) {
      geoLat_scgc.subcube(iLon, 0, iAlt, iLon, nLats-1, iAlt) = lat1d;
    }
  }

  fvec alt1d(nAlts);

  if (grid_input.IsUniformAlt) {
    for (iAlt = 0; iAlt < nAlts; iAlt++)
      alt1d(iAlt) =
        grid_input.alt_min +
        (iAlt-nGeoGhosts) * grid_input.dalt;
  }
  for (iLon = 0; iLon < nLons; iLon++) {
    for (iLat = 0; iLat < nLats; iLat++) {
      geoAlt_scgc.tube(iLon, iLat) = alt1d;
    }
  }

  IsGeoGrid = 1;

  // Calculate the radius, etc:

  fill_grid_radius(planet, report);

  fill_grid_bfield(planet, input, report);

  report.exit(function);
}
