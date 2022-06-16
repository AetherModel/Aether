// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>

#include "aether.h"
#include <math.h>

// ----------------------------------------------------------------------
// This function takes the normalized coordinates and makes latitude
// and longitude arrays from them.  It can do this for the corners or
// edges, depending on the offset.
// ----------------------------------------------------------------------

void fill_cubesphere_lat_lon_from_norms(Quadtree quadtree,
                                        arma_vec dr,
                                        arma_vec du,
                                        arma_vec ll,
                                        int64_t nGCs,
                                        precision_t left_off,
                                        precision_t down_off,
                                        arma_mat &lat2d,
                                        arma_mat &lon2d) {

  int64_t nX = lat2d.n_rows;
  int64_t nY = lat2d.n_cols;

  double xn, yn, zn, rn;
  double xp, yp, zp, rp, latp, lonp;

  double a = sqrt(3);

  arma_vec xyz, xyzn, xyz_wrapped;

  for (int iDU = 0; iDU < nY; iDU++) {
    for (int iLR = 0; iLR < nX; iLR++) {

      double iD = iDU - nGCs + down_off;
      double iL = iLR - nGCs + left_off;

      xyz = ll + dr * iL + du * iD;
      xyz_wrapped = quadtree.wrap_point_cubesphere(xyz) * a;
      xyzn = normalise(xyz_wrapped);
      xp = xyzn(0);
      yp = xyzn(1);
      zp = xyzn(2);

      latp = asin(zp);
      lonp = atan2(yp, xp) + 3 * cPI / 4;

      if (lonp > cTWOPI)
        lonp = lonp - cTWOPI;

      if (lonp < 0.0)
        lonp = lonp + cTWOPI;

      lat2d(iLR, iDU) = latp;
      lon2d(iLR, iDU) = lonp;
    }
  }
}


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

  dr = size_right_norm / (nLons - 2 * nGCs);
  du = size_up_norm / (nLats - 2 * nGCs);
  ll = lower_left_norm;

  int64_t iAlt, iLon, iLat;

  // ---------------------------------------------
  // Cell Centers
  // ---------------------------------------------
  arma_mat lat2d(nLons, nLats);
  arma_mat lon2d(nLons, nLats);
  fill_cubesphere_lat_lon_from_norms(quadtree, dr, du, ll, nGCs, 0.5, 0.5,
                                     lat2d, lon2d);
  for (iAlt = 0; iAlt < nAlts; iAlt++) {
    geoLon_scgc.slice(iAlt) = lon2d;
    geoLat_scgc.slice(iAlt) = lat2d;
  }

  // ---------------------------------------------
  // Left Sides - edges on left side (no offset left)
  // ---------------------------------------------
  arma_mat lat2d_left(nLons + 1, nLats);
  arma_mat lon2d_left(nLons + 1, nLats);
  fill_cubesphere_lat_lon_from_norms(quadtree, dr, du, ll, nGCs, 0.0, 0.5,
                                     lat2d_left, lon2d_left);
  for (iAlt = 0; iAlt < nAlts; iAlt++) {
    geoLon_Left.slice(iAlt) = lon2d_left;
    geoLat_Left.slice(iAlt) = lat2d_left;
  }

  // ---------------------------------------------
  // Down Sides - edges on down side (no offset down)
  // ---------------------------------------------
  arma_mat lat2d_down(nLons, nLats + 1);
  arma_mat lon2d_down(nLons, nLats + 1);
  fill_cubesphere_lat_lon_from_norms(quadtree, dr, du, ll, nGCs, 0.5, 0.0,
                                     lat2d_down, lon2d_down);
  for (iAlt = 0; iAlt < nAlts; iAlt++) {
    geoLon_Down.slice(iAlt) = lon2d_down;
    geoLat_Down.slice(iAlt) = lat2d_down;
  }

  // ---------------------------------------------
  // Corners (lower left) - no offsets
  // ---------------------------------------------
  arma_mat lat2d_corner(nLons + 1, nLats + 1);
  arma_mat lon2d_corner(nLons + 1, nLats + 1);
  fill_cubesphere_lat_lon_from_norms(quadtree, dr, du, ll, nGCs, 0.0, 0.0,
                                     lat2d_corner, lon2d_corner);
  for (iAlt = 0; iAlt < nAlts + 1; iAlt++) {
    geoLon_Corner.slice(iAlt) = lon2d_corner;
    geoLat_Corner.slice(iAlt) = lat2d_corner;
  }

  arma_vec alt1d(nAlts);

  Inputs::grid_input_struct grid_input = input.get_grid_inputs();

  if (grid_input.IsUniformAlt) {
    for (iAlt = 0; iAlt < nAlts; iAlt++)
      alt1d(iAlt) = grid_input.alt_min + (iAlt - nGeoGhosts) * grid_input.dalt;
  }

  arma_vec alt1d_below = calc_bin_edges(alt1d);

  for (iLon = 0; iLon < nLons; iLon++) {
    for (iLat = 0; iLat < nLats; iLat++) {
      geoAlt_scgc.tube(iLon, iLat) = alt1d;
      geoAlt_Below.tube(iLon, iLat) = alt1d_below;
    }
  }

  for (iLon = 0; iLon < nLons + 1; iLon++) {
    for (iLat = 0; iLat < nLats + 1; iLat++)
      geoAlt_Corner.tube(iLon, iLat) = alt1d_below;
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

  // Get some coordinates and sizes in normalized coordinates:
  arma_vec lower_left_norm = quadtree.get_vect("LL");
  arma_vec middle_norm = quadtree.get_vect("MID");
  arma_vec size_right_norm = quadtree.get_vect("SR");
  arma_vec size_up_norm = quadtree.get_vect("SU");

  // Move to the next block in 4 directions:
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

  // ---------------------------------------------
  // Left Sides - edges on left side (no offset left)
  // ---------------------------------------------
  arma_mat lat2d_left(nLons + 1, nLats);
  arma_mat lon2d_left(nLons + 1, nLats);
  for (iLat = 0; iLat < nLats; iLat++) {
    for (iLon = 0; iLon < nLons + 1; iLon++) {
      lat2d_left(iLon, iLat) = lat0 + (iLat - nGCs + 0.5) * dlat;
      lon2d_left(iLon, iLat) = lon0 + (iLon - nGCs) * dlon;
    }
  }
  for (iAlt = 0; iAlt < nAlts; iAlt++) {
    geoLon_Left.slice(iAlt) = lon2d_left;
    geoLat_Left.slice(iAlt) = lat2d_left;
  }

  // ---------------------------------------------
  // Down Sides - edges on down side (no offset lat)
  // ---------------------------------------------
  arma_mat lat2d_down(nLons, nLats + 1);
  arma_mat lon2d_down(nLons, nLats + 1);
  for (iLat = 0; iLat < nLats + 1; iLat++) {
    for (iLon = 0; iLon < nLons; iLon++) {
      lat2d_down(iLon, iLat) = lat0 + (iLat - nGCs) * dlat;
      lon2d_down(iLon, iLat) = lon0 + (iLon - nGCs + 0.5) * dlon;
    }
  }
  for (iAlt = 0; iAlt < nAlts; iAlt++) {
    geoLon_down.slice(iAlt) = lon2d_down;
    geoLat_down.slice(iAlt) = lat2d_down;
  }

  // ---------------------------------------------
  // Corner Sides - corner (no offset lat or lon)
  // ---------------------------------------------
  arma_mat lat2d_corner(nLons + 1, nLats + 1);
  arma_mat lon2d_corner(nLons + 1, nLats + 1);
  for (iLat = 0; iLat < nLats + 1; iLat++) {
    for (iLon = 0; iLon < nLons + 1; iLon++) {
      lat2d_corner(iLon, iLat) = lat0 + (iLat - nGCs) * dlat;
      lon2d_corner(iLon, iLat) = lon0 + (iLon - nGCs) * dlon;
    }
  }
  for (iAlt = 0; iAlt < nAlts + 1; iAlt++) {
    geoLon_corner.slice(iAlt) = lon2d_corner;
    geoLat_corner.slice(iAlt) = lat2d_corner;
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
