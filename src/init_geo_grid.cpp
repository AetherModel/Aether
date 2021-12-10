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

  IsLatLonGrid = true;
  
  // This is just an example:

  Inputs::grid_input_struct grid_input = input.get_grid_inputs();

  int64_t iLon, iLat, iAlt;

  // Figure out dlon and lon0:
  int64_t nLonsTotal = (nLons - 2 * nGCs) * input.get_nBlocksLonGeo();
  int64_t iBlockLon = iGrid % input.get_nBlocksLonGeo();
  precision_t dlon = (grid_input.lon_max - grid_input.lon_min) / nLonsTotal;
  precision_t lon0 = grid_input.lon_min + iBlockLon * (nLons - 2 * nGCs) * dlon;
  arma_vec lon1d(nLons);

  // Longitudes:
  // - Make a 1d vector
  // - copy it into the 3d cube
  for (iLon = 0; iLon < nLons; iLon++)
    lon1d(iLon) = lon0 + (iLon - nGCs + 0.5) * dlon;

  for (iLat = 0; iLat < nLats; iLat++) {
    for (iAlt = 0; iAlt < nAlts; iAlt++)
      geoLon_scgc.subcube(0, iLat, iAlt, nLons - 1, iLat, iAlt) = lon1d;
  }

  // Figure out dlat and lat0:
  int64_t nLatsTotal = (nLats - 2 * nGCs) * input.get_nBlocksLatGeo();
  int64_t iBlockLat = iGrid / input.get_nBlocksLonGeo();
  precision_t dlat = (grid_input.lat_max - grid_input.lat_min) / nLatsTotal;
  precision_t lat0 = grid_input.lat_min + iBlockLat * (nLats - 2 * nGCs) * dlat;
  arma_vec lat1d(nLats);

  // Latitudes:
  // - Make a 1d vector
  // - copy it into the 3d cube
  for (iLat = 0; iLat < nLats; iLat++)
    lat1d(iLat) = lat0 + (iLat - nGCs + 0.5) * dlat;

  if (report.test_verbose(2)) {
    std::cout << "Grid Initialization :\n";
    std::cout << "  lon0 : " << lon0 * cRtoD
	      << "; lat0 : " << lat0 * cRtoD << "\n";
  }
  
  for (iLon = 0; iLon < nLons; iLon++) {
    for (iAlt = 0; iAlt < nAlts; iAlt++)
      geoLat_scgc.subcube(iLon, 0, iAlt, iLon, nLats - 1, iAlt) = lat1d;
  }

  arma_vec alt1d(nAlts);

  if (grid_input.IsUniformAlt) {
    for (iAlt = 0; iAlt < nAlts; iAlt++)
      alt1d(iAlt) = grid_input.alt_min + (iAlt - nGeoGhosts) * grid_input.dalt;
  }

  for (iLon = 0; iLon < nLons; iLon++) {
    for (iLat = 0; iLat < nLats; iLat++)
      geoAlt_scgc.tube(iLon, iLat) = alt1d;
  }

  // --------------------------------------------------------------
  // Figure out message passing connectivity information

  int nCols = input.get_nBlocksLonGeo();
  int iRow = iProc / nCols;
  int iLeftMost = iRow * nCols;
  int iRightMost = iLeftMost + nCols-1;
  int iCol = iProc % nCols;

  // Determine y-minus (connectivity to south):
  
  // Figure out if the grid is touching the south pole:
  if (compare(lat0, -cPI/2)) {
    DoesTouchSouthPole = true;
    // Touching the south pole -> just shift by 180 degrees:
    iProcYm = iLeftMost + (iCol + nCols/2) % nCols;
  } else {
    DoesTouchSouthPole = false;
    if (iBlockLat == 0)
      // Not touching the south pole, but no blocks to the south:
      iProcYm = -1;
    else
      // Not touching the south pole, and there is a block to the south:
      iProcYm = iProc - input.get_nBlocksLonGeo();
  }

  // Determine y-plus (connectivity to north):

  // Figure out if the grid is touching the north pole:
  precision_t latmax = lat0 + (nLats - 2 * nGCs) * dlat;
  if (compare(latmax, cPI/2)) {
    DoesTouchNorthPole = true;
    // Touching the north pole -> just shift by 180 degrees:
    iProcYp = iLeftMost + (iCol + nCols/2) % nCols;
  } else {
    DoesTouchNorthPole = false;
    if (iBlockLat == input.get_nBlocksLatGeo()-1)
      // Not touching the north pole, but no blocks to the north:
      iProcYp = -1;
    else
      // Not touching the north pole, and there is a block to the north:
      iProcYp = iProc + input.get_nBlocksLonGeo();
  }

  // Determine x-minus (connectivity to west):
  
  if (iCol > 0) {
    iProcXm = iProc - 1;
  } else {
    // Figure out whether grid is touching East/West:
    if (compare(grid_input.lon_min, 0) &&
	compare(grid_input.lon_max, cTWOPI)) {
      // wrapped completely around the Earth
      if (iCol == 0)
	// wrap around to the right-most point:
	iProcXm = iRightMost;
    } else {
      // doesn't wrap around the earth, so need fixed boundaries
      if (iBlockLon == 0)
	// Fixed BCs to the west:
	iProcXm = -1;
    }
  }
      
  // Determine x-plus (connectivity to east):
  
  if (iCol < nCols-1) {
    iProcXp = iProc + 1;
  } else {
    // Figure out whether grid is touching East/West:
    if (compare(grid_input.lon_min, 0) &&
	compare(grid_input.lon_max, cTWOPI)) {
      // wrapped completely around the Earth
      if (iBlockLon == nCols-1)
	// wrap around to the left-most point:
	iProcXp = iLeftMost;
    } else {
      // doesn't wrap around the earth, so need fixed boundaries
      if (iBlockLon == nCols-1)
	// Fixed BCs to the east:
	iProcXp = -1;
    }
  }

  if (report.test_verbose(0)) 
    std::cout << "connectivity : "
	      << "  iProc : " << iProc
	      << "  isnorth : " << DoesTouchNorthPole
	      << "  issouth : " << DoesTouchSouthPole
	      << "  iProcYm : " << iProcYm
	      << "  iProcYp : " << iProcYp
	      << "  iProcXm : " << iProcXm
	      << "  iProcXp : " << iProcXp
	      << "\n";
  
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
