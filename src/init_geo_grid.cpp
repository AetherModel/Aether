// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>

#include "aether.h"
#include <math.h>

// ----------------------------------------------------------------------
// Create connectivity between the nodes for message passing for cubesphere
// ----------------------------------------------------------------------

void Grid::create_cubesphere_connection(Quadtree quadtree,
                                        Inputs input,
                                        Report &report) {

  std::string function = "Grid::create_cubesphere_connection";
  static int iFunction = -1;
  report.enter(function, iFunction);

  IsLatLonGrid = false;

  arma_vec lower_left_norm = quadtree.get_vect("LL");
  arma_vec middle_norm = quadtree.get_vect("MID");
  arma_vec size_right_norm = quadtree.get_vect("SR");
  arma_vec size_up_norm = quadtree.get_vect("SU");

  // These points go off the edge to the next block in each direction:
  arma_vec down_norm = middle_norm - 0.51 * size_up_norm;
  arma_vec up_norm = middle_norm + 0.51 * size_up_norm;
  arma_vec left_norm = middle_norm - 0.51 * size_right_norm;
  arma_vec right_norm = middle_norm + 0.51 * size_right_norm;

  // Find those points in the quadtree to figure out which processor
  // they are on
  iProcYm = quadtree.find_point(down_norm) + iMember * nGrids;
  iProcYp = quadtree.find_point(up_norm) + iMember * nGrids;
  iProcXm = quadtree.find_point(left_norm) + iMember * nGrids;
  iProcXp = quadtree.find_point(right_norm) + iMember * nGrids;

  // Need to know which side the current block is on and which side each
  // of the blocks in the different directions is on.  Need this so we can
  // know how to unpack the variables after the message pass.
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

  report.exit(function);
  return;
}

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

  // Loop through each point and derive the coordinate
  for (int iDU = 0; iDU < nY; iDU++) {
    for (int iLR = 0; iLR < nX; iLR++) {

      // the offsets are so we can find cell centers, edges, and corners
      double iD = iDU - nGCs + down_off;
      double iL = iLR - nGCs + left_off;

      // This is the normalized coordinate:
      xyz = ll + dr * iL + du * iD;
      // Ghost cells could be off the edge, so wrap to other face:
      xyz_wrapped = quadtree.wrap_point_cubesphere(xyz) * a;
      // Normalize the coordinate to a unit vector:
      xyzn = normalise(xyz_wrapped);
      xp = xyzn(0);
      yp = xyzn(1);
      zp = xyzn(2);

      // Derive lat and lon from unit vector:
      latp = asin(zp);
      // offset for lon is to put the left edge of face 0 at 0 longitude:
      lonp = atan2(yp, xp) + 3 * cPI / 4;

      if (lonp > cTWOPI)
        lonp = lonp - cTWOPI;

      if (lonp < 0.0)
        lonp = lonp + cTWOPI;

      lat2d(iLR, iDU) = latp;
      lon2d(iLR, iDU) = lonp;
    }
  }

  return;
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

  arma_vec dr(3), du(3), ll(3);
  double xn, yn, zn, rn;
  double xp, yp, zp, rp, latp, lonp;

  double a = sqrt(3);

  arma_vec lower_left_norm = quadtree.get_vect("LL");
  arma_vec size_right_norm = quadtree.get_vect("SR");
  arma_vec size_up_norm = quadtree.get_vect("SU");

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

  report.exit(function);
  return;
}

// ----------------------------------------------------------------------
// Create connectivity between the nodes for message passing for sphere
// ----------------------------------------------------------------------

void Grid::create_sphere_connection(Quadtree quadtree,
                                    Inputs input,
                                    Report &report) {

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

void Grid::create_sphere_grid(Quadtree quadtree,
                              Inputs input,
                              Report &report) {

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
    geoLon_Down.slice(iAlt) = lon2d_down;
    geoLat_Down.slice(iAlt) = lat2d_down;
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
    geoLon_Corner.slice(iAlt) = lon2d_corner;
    geoLat_Corner.slice(iAlt) = lat2d_corner;
  }

  report.exit(function);
  return;
}

// ----------------------------------------------------------------------
// Create a spherical grid with lon/lat/alt coordinates
// ----------------------------------------------------------------------

void Grid::create_altitudes(Planets planet, Inputs input, Report &report) {

  std::string function = "Grid::create_altitudes";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int64_t iLon, iLat, iAlt;

  arma_vec alt1d(nAlts);

  Inputs::grid_input_struct grid_input = input.get_grid_inputs();

  if (grid_input.IsUniformAlt) {
    for (iAlt = 0; iAlt < nAlts; iAlt++)
      alt1d(iAlt) = grid_input.alt_min + (iAlt - nGeoGhosts) * grid_input.dalt;
  } else {

    json neutrals = planet.get_neutrals();
    json temperatures = planet.get_temperatures();
    std::vector<double> input_alt;
    std::vector<double> input_temp;

    for (int i = 0; i < temperatures["alt"].size(); i++) {
      input_alt.push_back(double(temperatures["alt"][i]) * 1000.0);
      input_temp.push_back(temperatures["temp"][i]);
    }

    precision_t scale_height, temperature, gravity, radius, mass, density;
    int64_t nSp = neutrals["name"].size();
    arma_vec densities(nSp);
    arma_vec masses(nSp);
    arma_vec h(nSp);

    int64_t iSp;

    report.print(1, "Making non-uniform altitude grid!");

    if (grid_input.dalt > 0.5) {
      if (report.test_verbose(0)) {
        std::cout << "-----------------------------------------------------\n";
        std::cout << "WARNING: dAlt is set to > 0.5, with non-uniform grid!\n";
        std::cout << "   dAlt = " << grid_input.dalt << "\n";
        std::cout << "-----------------------------------------------------\n";
      }
    }

    double alt = grid_input.alt_min;
    radius = planet.get_radius(0.0,input) + alt;
    precision_t mu = planet.get_mu();
    gravity = mu / (radius * radius);

    temperature = interpolate_1d(alt, input_alt, input_temp);

    mass = 0.0;
    density = 0.0;

    for (iSp = 0; iSp < nSp; iSp++) {
      masses(iSp) = double(neutrals["mass"][iSp]) * cAMU;
      densities[iSp] = neutrals["BC"][iSp];
      h(iSp) = cKB * temperature / (masses(iSp) * gravity);
      mass = mass + masses(iSp) * densities[iSp];
      density = density + densities[iSp];
    }

    // convert mass density into mass:
    mass = mass / density;
    scale_height = cKB * temperature / (mass * gravity);

    precision_t dalt = scale_height * grid_input.dalt;
    precision_t dAltLimiter = dalt * 10.0;

    // Fills bottom ghost cells with constant dAlt
    // Fills bottom cell with actual desired bottom altitude
    for (iAlt = 0; iAlt <= nGeoGhosts; iAlt++) {
      alt1d(iAlt) = grid_input.alt_min + (iAlt - nGeoGhosts) * dalt;

      if (report.test_verbose(1))
        std::cout << "iAlt : " << iAlt
                  << " Altitude : " << alt1d(iAlt) / 1000.0
                  << " (km)\n";
    }

    for (iAlt = nGeoGhosts + 1; iAlt < nAlts; iAlt++) {

      alt = alt1d(iAlt - 1);
      temperature = interpolate_1d(alt, input_alt, input_temp);
      radius = planet.get_radius(0.0, input) + alt;
      gravity = mu / (radius * radius);

      mass = 0.0;
      density = 0.0;

      for (iSp = 0; iSp < nSp; iSp++) {
        mass = mass + masses(iSp) * densities[iSp];
        density = density + densities[iSp];
      }

      // convert mass density into mass:
      mass = mass / density;
      scale_height = cKB * temperature / (mass * gravity);

      dalt = scale_height * grid_input.dalt;

      if (dalt > dAltLimiter)
        dalt = dAltLimiter;

      alt1d(iAlt) = alt + dalt;

      h = cKB * temperature / (masses * gravity);
      densities = densities % exp(-dalt / h);

      if (report.test_verbose(1))
        std::cout << "iAlt : " << iAlt
                  << " Altitude : " << alt1d(iAlt) / 1000.0
                  << " (km)\n";
    }
  }

  // This takes cell centers and calculates the edges:
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

  report.exit(function);
  return;
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

  if (input.get_is_cubesphere())
    create_cubesphere_connection(quadtree, input, report);
  else
    create_sphere_connection(quadtree, input, report);

  if (input.get_do_restart()) {
    report.print(1, "Restarting! Reading grid files!");
    DidWork = read_restart(input.get_restartin_dir());
  } else {
    if (input.get_is_cubesphere())
      create_cubesphere_grid(quadtree, input, report);
    else
      create_sphere_grid(quadtree, input, report);

    MPI_Barrier(aether_comm);
    create_altitudes(planet, input, report);

    DidWork = write_restart(input.get_restartout_dir());
  }

  // Calculate the radius (for spherical or non-spherical)
  fill_grid_radius(planet, input, report);
  // Calculate grid spacing
  calc_grid_spacing(planet, report);
  //calculate radial unit vector (for spherical or oblate planet)
  calc_rad_unit(planet, input, report);
  // Calculate gravity (including J2 term, if desired)
  calc_gravity(planet, input, report);

  // Calculate magnetic field and magnetic coordinates:
  fill_grid_bfield(planet, input, report);

  // Throw a little message for students:
  report.student_checker_function_name(input.get_is_student(),
				       input.get_student_name(),
				       4, "");
  
  report.exit(function);
  return DidWork;
}
