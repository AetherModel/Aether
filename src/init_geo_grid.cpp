// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>

#include "aether.h"
#include <math.h>

// ----------------------------------------------------------------------
// Create a geographic grid
//    - if restarting, read in the grid
//    - if not restarting, initialize the grid
// ----------------------------------------------------------------------

void Grid::create_cubesphere_grid(Quadtree quadtree,
                                  Inputs input,
                                  Report &report) {

  std::string function = "Grid::create_cubesphere_grid";
  static int iFunction = -1;
  report.enter(function, iFunction);

  IsLatLonGrid = false;

  arma_vec lower_left_norm = quadtree.get_vect("LL");
  arma_vec middle_norm = quadtree.get_vect("MID");
  arma_vec size_right_norm = quadtree.get_vect("SR");
  arma_vec size_up_norm = quadtree.get_vect("SU");

  arma_vec down_norm = middle_norm - size_up_norm;
  arma_vec up_norm = middle_norm + size_up_norm;
  arma_vec left_norm = middle_norm - size_right_norm;
  arma_vec right_norm = middle_norm + size_right_norm;

  iProcYm = quadtree.find_point(down_norm);
  iProcYp = quadtree.find_point(up_norm);
  iProcXm = quadtree.find_point(left_norm);
  iProcXp = quadtree.find_point(right_norm);

  iRoot = quadtree.find_root(middle_norm);
  iRootYm = quadtree.find_root(down_norm);
  iRootYp = quadtree.find_root(up_norm);
  iRootXm = quadtree.find_root(left_norm);
  iRootXp = quadtree.find_root(right_norm);

  // These should be the exact edge of the face.
  // The from and to processors should get these in the same place,
  // so they can be used to match which processor to send / receive info
  edge_Xp = middle_norm + size_right_norm / 2.0;
  edge_Xm = middle_norm - size_right_norm / 2.0;
  edge_Yp = middle_norm + size_up_norm / 2.0;
  edge_Ym = middle_norm - size_up_norm / 2.0;

  if (report.test_verbose(2))
    std::cout << "connectivity : "
              << "  iProc : " << iProc << "\n"
              << "  isnorth : " << DoesTouchNorthPole << "\n"
              << "  issouth : " << DoesTouchSouthPole << "\n"
              << "  iProcXp : " << iProcXp << "\n"
              << "  iProcYp : " << iProcYp << "\n"
              << "  iProcXm : " << iProcXm << "\n"
              << "  iProcYm : " << iProcYm << "\n"
              << "  iRoot   : " << iRoot << "\n"
              << "  iRootXp : " << iRootXp << "\n"
              << "  iRootYp : " << iRootYp << "\n"
              << "  iRootXm : " << iRootXm << "\n"
              << "  iRootYm : " << iRootYm << "\n";

  int64_t iProcSelf = quadtree.find_point(middle_norm);

  arma_vec dr(3), du(3), ll(3);
  double xn, yn, zn, rn;
  double xp, yp, zp, rp, latp, lonp;

  double a = sqrt(3);

  dr = size_right_norm * a / (nLons - 2 * nGCs);
  du = size_up_norm * a / (nLats - 2 * nGCs);
  ll = lower_left_norm * a;

  arma_mat lat2d(nLons, nLats), lon2d(nLons, nLats);

  for (int iDU = 0; iDU < nY; iDU++) {
    for (int iLR = 0; iLR < nX; iLR++) {

      double iD = iDU - nGCs + 0.5;
      double iL = iLR - nGCs + 0.5;

      xn = ll(0) + dr(0) * iL + du(0) * iD;
      yn = ll(1) + dr(1) * iL + du(1) * iD;
      zn = ll(2) + dr(2) * iL + du(2) * iD;
      rn = sqrt(xn * xn + yn * yn + zn * zn);

      xp = xn / rn;
      yp = yn / rn;
      zp = zn / rn;
      rp = sqrt(xp * xp + yp * yp + zp * zp);

      latp = asin(zp / rp);
      lonp = atan2(yp, xp) + 3 * cPI / 4;

      if (lonp > cTWOPI)
        lonp = lonp - cTWOPI;

      if (lonp < 0.0)
        lonp = lonp + cTWOPI;

      lat2d(iLR, iDU) = latp;
      lon2d(iLR, iDU) = lonp;

    }
  }

  int64_t iAlt, iLon, iLat;

  for (iAlt = 0; iAlt < nAlts; iAlt++) {
    geoLon_scgc.slice(iAlt) = lon2d;
    geoLat_scgc.slice(iAlt) = lat2d;
  }

  arma_vec alt1d(nAlts);

  Inputs::grid_input_struct grid_input = input.get_grid_inputs();

  if (grid_input.IsUniformAlt) {
    for (iAlt = 0; iAlt < nAlts; iAlt++)
      alt1d(iAlt) = grid_input.alt_min + (iAlt - nGeoGhosts) * grid_input.dalt;
  }

  for (iLon = 0; iLon < nLons; iLon++) {
    for (iLat = 0; iLat < nLats; iLat++)
      geoAlt_scgc.tube(iLon, iLat) = alt1d;
  }

}

// ----------------------------------------------------------------------
// Create a geographic grid
//    - if restarting, read in the grid
//    - if not restarting, initialize the grid
// ----------------------------------------------------------------------

void Grid::create_simple_lat_lon_alt_grid(Quadtree quadtree,
                                          Inputs input,
                                          Report &report) {

  std::string function = "Grid::create_simple_lat_lon_alt_grid";
  static int iFunction = -1;
  report.enter(function, iFunction);

  IsLatLonGrid = true;

  arma_vec lower_left_norm = quadtree.get_vect("LL");
  arma_vec middle_norm = quadtree.get_vect("MID");
  arma_vec size_right_norm = quadtree.get_vect("SR");
  arma_vec size_up_norm = quadtree.get_vect("SU");

  std::cout << "middle : " << iProc << " "
            << middle_norm(0) << " "
            << middle_norm(1) << " "
            << middle_norm(2) << "\n";

  arma_vec down_norm = middle_norm - size_up_norm;
  arma_vec up_norm = middle_norm + size_up_norm;
  arma_vec left_norm = middle_norm - size_right_norm;
  arma_vec right_norm = middle_norm + size_right_norm;

  // The first component could wrap around:
  right_norm(0) = fmod(right_norm(0), quadtree.limit_high(0));
  left_norm(0) = fmod((left_norm(0) + quadtree.limit_high(0)),
                      quadtree.limit_high(0));

  // These should be the exact edge of the face.
  // The from and to processors should get these in the same place,
  // so they can be used to match which processor to send / receive info
  edge_Xp = middle_norm + size_right_norm / 2.0;
  // wrap in longitude:
  edge_Xp(0) = fmod(edge_Xp(0), quadtree.limit_high(0));
  edge_Xm = middle_norm - size_right_norm / 2.0;
  edge_Yp = middle_norm + size_up_norm / 2.0;
  edge_Ym = middle_norm - size_up_norm / 2.0;

  iProcYm = quadtree.find_point(down_norm);
  iProcYp = quadtree.find_point(up_norm);
  iProcXm = quadtree.find_point(left_norm);
  iProcXp = quadtree.find_point(right_norm);

  iRoot = quadtree.find_root(middle_norm);
  iRootYm = quadtree.find_root(down_norm);
  iRootYp = quadtree.find_root(up_norm);
  iRootXm = quadtree.find_root(left_norm);
  iRootXp = quadtree.find_root(right_norm);

  // Check if touching South Pole:
  if (lower_left_norm(1) == quadtree.limit_low(1)) {
    DoesTouchSouthPole = true;

    // edges need to be adjusted to deal with longitudes, since the
    // pole will 180deg different for the from and to processors
    if (edge_Ym(0) < 1.0)
      edge_Ym(0) += 0.5;
    else
      edge_Ym(0) -= 0.5;
  }

  // Check if touching North Pole:
  if (lower_left_norm(1) + size_up_norm(1) == quadtree.limit_high(1)) {
    DoesTouchNorthPole = true;

    // edge need to be adjusted to deal with longitudes, since the
    // pole will 180deg different for the from and to processors
    if (edge_Yp(0) < 1.0)
      edge_Yp(0) += 0.5;
    else
      edge_Yp(0) -= 0.5;
  }

  Inputs::grid_input_struct grid_input = input.get_grid_inputs();

  int64_t iLon, iLat, iAlt;

  precision_t dlon = size_right_norm(0) * cPI / (nLons - 2 * nGCs);
  precision_t lon0 = lower_left_norm(0) * cPI;
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

  precision_t dlat = size_up_norm(1) * cPI / (nLats - 2 * nGCs);
  precision_t lat0 = lower_left_norm(1) * cPI;
  arma_vec lat1d(nLats);

  // Latitudes:
  // - Make a 1d vector
  // - copy it into the 3d cube
  for (iLat = 0; iLat < nLats; iLat++)
    lat1d(iLat) = lat0 + (iLat - nGCs + 0.5) * dlat;

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

  if (report.test_verbose(2))
    std::cout << "connectivity : "
              << "  iProc : " << iProc << "\n"
              << "  isnorth : " << DoesTouchNorthPole << "\n"
              << "  issouth : " << DoesTouchSouthPole << "\n"
              << "  iProcYm : " << iProcYm << "\n"
              << "  iProcYp : " << iProcYp << "\n"
              << "  iProcXm : " << iProcXm << "\n"
              << "  iProcXp : " << iProcXp << "\n";

  report.exit(function);
}


// ----------------------------------------------------------------------
// Initialize the geographic grid.  At the moment, this is a simple
// Lon/Lat/Alt grid.  The grid structure is general enough that each
// of the lon, lat, and alt can be a function of the other variables.
// ----------------------------------------------------------------------

bool Grid::init_geo_grid(Quadtree quadtree,
                         Planets planet,
                         Inputs input,
                         Report &report) {

  std::string function = "Grid::init_geo_grid";
  static int iFunction = -1;
  report.enter(function, iFunction);
  bool DidWork = true;

  IsGeoGrid = 1;

  // Logic here is flawed:
  //   - the "create" functions both make the grid and build connections
  //   - the restart simply builds the grid.
  // Need to separate the two.

  if (input.get_do_restart()) {
    report.print(1, "Restarting! Reading grid files!");
    DidWork = read_restart(input.get_restartin_dir());
  } else {
    if (input.get_is_cubesphere())
      create_cubesphere_grid(quadtree, input, report);
    else
      create_simple_lat_lon_alt_grid(quadtree, input, report);

    DidWork = write_restart(input.get_restartout_dir());
  }

  // Calculate the radius, etc:
  fill_grid_radius(planet, report);

  // Calculate magnetic field and magnetic coordinates:
  fill_grid_bfield(planet, input, report);

  report.exit(function);
  return DidWork;
}
