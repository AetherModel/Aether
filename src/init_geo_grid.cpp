// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>

#include "aether.h"

// ----------------------------------------------------------------------
// Create a geographic grid
//    - if restarting, read in the grid
//    - if not restarting, initialize the grid
// ----------------------------------------------------------------------

void Grid::create_simple_lat_lon_alt_grid(Inputs input, Report &report) {

  std::string function = "Grid::create_simple_lat_lon_alt_grid";
  static int iFunction = -1;
  report.enter(function, iFunction);

  // This is just an example:

  Inputs::grid_input_struct grid_input = input.get_grid_inputs();

  int64_t iLon, iLat, iAlt;

  // Longitudes:
  // - Make a 1d vector
  // - copy it into the 3d cube
  arma_vec lon1d(nLons);
  precision_t dlon = (grid_input.lon_max - grid_input.lon_min) /
                     (nLons - 2 * nGCs);

  for (iLon = 0; iLon < nLons; iLon++)
    lon1d(iLon) = grid_input.lon_min + (iLon - nGCs + 0.5) * dlon;

  for (iLat = 0; iLat < nLats; iLat++) {
    for (iAlt = 0; iAlt < nAlts; iAlt++)
      geoLon_scgc.subcube(0, iLat, iAlt, nLons - 1, iLat, iAlt) = lon1d;
  }

  // Latitudes:
  // - Make a 1d vector
  // - copy it into the 3d cube
  arma_vec lat1d(nLats);
  precision_t dlat = (grid_input.lat_max - grid_input.lat_min) /
                     (nLats - 2 * nGCs);

  for (iLat = 0; iLat < nLats; iLat++)
    lat1d(iLat) = grid_input.lat_min + (iLat - nGCs + 0.5) * dlat;

  for (iLon = 0; iLon < nLons; iLon++) {
    for (iAlt = 0; iAlt < nAlts; iAlt++)
      geoLat_scgc.subcube(iLon, 0, iAlt, iLon, nLats - 1, iAlt) = lat1d;
  }

  arma_vec alt1d(nAlts);

  if (grid_input.IsUniformAlt) {
    for (iAlt = 0; iAlt < nAlts; iAlt++)
      alt1d(iAlt) =
        grid_input.alt_min +
        (iAlt - nGeoGhosts) * grid_input.dalt;
  }

  for (iLon = 0; iLon < nLons; iLon++) {
    for (iLat = 0; iLat < nLats; iLat++)
      geoAlt_scgc.tube(iLon, iLat) = alt1d;
  }

  report.exit(function);
}


// ----------------------------------------------------------------------
// Initialize the geographic grid.  At the moment, this is a simple
// Lon/Lat/Alt grid.  The grid structure is general enough that each
// of the lon, lat, and alt can be a function of the other variables.
// ----------------------------------------------------------------------

void Grid::init_geo_grid(Planets planet, Inputs input, Report &report) {

  std::string function = "Grid::init_geo_grid";
  static int iFunction = -1;
  report.enter(function, iFunction);
  bool DidWork = true;

  IsGeoGrid = 1;

  if (input.get_do_restart()) {
    report.print(1, "Restarting! Reading grid files!");
    DidWork = read_restart(input.get_restartin_dir());
  } else {
    create_simple_lat_lon_alt_grid(input, report);
    DidWork = write_restart(input.get_restartout_dir());
  }

  if (report.test_verbose(3)) {
    std::cout << "init_geo_grid testing\n";
    std::cout << "   DidWork (reading/writing restart files): "
              << DidWork << "\n";
    report_grid_boundaries();
  }

  // Calculate the radius, etc:
  fill_grid_radius(planet, report);

  // Calculate magnetic field and magnetic coordinates:
  fill_grid_bfield(planet, input, report);

  report.exit(function);
}
