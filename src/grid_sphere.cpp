// Copyright 2025, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "aether.h"


// ----------------------------------------------------------------------
// Create connectivity between the nodes for message passing for sphere
// ----------------------------------------------------------------------

void Grid::create_sphere_connection(Quadtree quadtree) {

  std::string function = "Grid::create_sphere_connection";
  static int iFunction = -1;
  report.enter(function, iFunction);

  IsLatLonGrid = true;

  // Get some coordinates and sizes in normalized coordinates:
  arma_vec lower_left_norm = quadtree.get_vect("LL");
  arma_vec middle_norm = quadtree.get_vect("MID");
  arma_vec size_right_norm = quadtree.get_vect("SR");
  arma_vec size_up_norm = quadtree.get_vect("SU");

  // Move to the next block in 4 directions:
  arma_vec down_norm = middle_norm - 0.51 * size_up_norm;
  arma_vec up_norm = middle_norm + 0.51 * size_up_norm;
  arma_vec left_norm = middle_norm - 0.51 * size_right_norm;
  arma_vec right_norm = middle_norm + 0.51 * size_right_norm;

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

  iProcYm = quadtree.find_point(down_norm) + iMember * nGrids;
  iProcYp = quadtree.find_point(up_norm) + iMember * nGrids;
  iProcXm = quadtree.find_point(left_norm) + iMember * nGrids;
  iProcXp = quadtree.find_point(right_norm) + iMember * nGrids;

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
  return;
}

// ----------------------------------------------------------------------
// Create a spherical grid with lon/lat/alt coordinates
// ----------------------------------------------------------------------

void Grid::create_sphere_grid(Quadtree quadtree) {

  std::string function = "Grid::create_simple_lat_lon_alt_grid";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int64_t iLon, iLat, iAlt;

  // Get some coordinates and sizes in normalized coordinates:
  arma_vec lower_left_norm = quadtree.get_vect("LL");
  arma_vec size_right_norm = quadtree.get_vect("SR");
  arma_vec size_up_norm = quadtree.get_vect("SU");

  precision_t dlon = size_right_norm(0) * cPI / (nLons - 2 * nGCs);
  precision_t lon0 = lower_left_norm(0) * cPI;
  arma_vec lon1d(nLons);

  // if we are not doing anything in the lon direction, then set dlon to
  // something reasonable:
  if (!HasXdim)
    dlon = 1.0 * cDtoR;

  // Longitudes:
  // - Make a 1d vector
  // - copy it into the 3d cube
  for (iLon = 0; iLon < nLons; iLon++)
    lon1d(iLon) = lon0 + (iLon - nGCs + 0.5) * dlon;

  if (report.test_verbose(1)) {
    std::cout << function << ": " << lon0 << " " << dlon << "\n";
    display_vector("in function " + function + " lon1d : ", lon1d * cRtoD);
  }

  for (iLat = 0; iLat < nLats; iLat++) {
    for (iAlt = 0; iAlt < nAlts; iAlt++) {
      geoLon_scgc.subcube(0, iLat, iAlt, nLons - 1, iLat, iAlt) = lon1d;
      i_center_scgc.subcube(0, iLat, iAlt, nLons - 1, iLat, iAlt) = lon1d;
    }
  }

  precision_t dlat = size_up_norm(1) * cPI / (nLats - 2 * nGCs);
  precision_t lat0 = lower_left_norm(1) * cPI;
  arma_vec lat1d(nLats);

  // if we are not doing anything in the lat direction, then set dlat to
  // something reasonable:
  if (!HasYdim)
    dlat = 1.0 * cDtoR;

  // Latitudes:
  // - Make a 1d vector
  // - copy it into the 3d cube
  for (iLat = 0; iLat < nLats; iLat++)
    lat1d(iLat) = lat0 + (iLat - nGCs + 0.5) * dlat;

  if (report.test_verbose(1)) {
    std::cout << function << ": " << lat0 << " " << dlat << "\n";

    display_vector("in function " + function + " lat1d : ", lat1d * cRtoD);
  }

  for (iLon = 0; iLon < nLons; iLon++) {
    for (iAlt = 0; iAlt < nAlts; iAlt++) {
      geoLat_scgc.subcube(iLon, 0, iAlt, iLon, nLats - 1, iAlt) = lat1d;
      j_center_scgc.subcube(iLon, 0, iAlt, iLon, nLats - 1, iAlt) = lat1d;
    }
  }

  arma_cube cos_lat = cos(geoLat_scgc);
  cos_lat.clamp(0.0001, 1.0);

  y_Center = geoLat_scgc;
  x_Center = geoLon_scgc % cos_lat;
  cell_area = dlat * dlon * cos_lat;

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
    i_edge_scgc.slice(iAlt) = lon2d_left;
  }

  arma_cube cos_lat_L = cos(geoLat_Left);
  cos_lat_L.clamp(0.0001, 1.0);

  x_Left = geoLon_Left % cos_lat_L;
  dy_Left.set_size(nLons, nLats, nAlts);
  dy_Left.fill(dlat);

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
    geoLon_Down.slice(iAlt) = lon2d_down;
    geoLat_Down.slice(iAlt) = lat2d_down;
    j_edge_scgc.slice(iAlt) = lat2d_down;
  }

  arma_cube cos_lat_D = cos(geoLat_Down);
  cos_lat_D.clamp(0.0001, 1.0);

  y_Down = geoLat_Down;
  dx_Down = dlon * cos_lat_D;

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
    geoLon_Corner.slice(iAlt) = lon2d_corner;
    geoLat_Corner.slice(iAlt) = lat2d_corner;
    i_corner_scgc.slice(iAlt) = lon2d_corner;
    j_corner_scgc.slice(iAlt) = lat2d_corner;
  }

  report.exit(function);
  return;
}
