// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "aether.h"

// Hepler variables / function begins.
// These are only used inside this cpp file and neither declared
// nor visible in any other file

// The size of a 2*2*2 arma cube
const arma::SizeCube unit_cube_size = arma::size(2, 2, 2);

// --------------------------------------------------------------------------
// Return the first index of two vectors on which they have different values
// --------------------------------------------------------------------------

int64_t first_diff_index(const arma_vec &a, const arma_vec &b) {
  int64_t i;

  for (i = 0; i < std::min(a.n_rows, b.n_rows); ++i) {
    if (std::abs(a[i] - b[i]) > cSmall)
      return i;
  }

  return i;
}

// --------------------------------------------------------------------------
// Assign any point on the surface of a cube a number within [0,5]
// --------------------------------------------------------------------------

int64_t get_cube_surface_number(precision_t x_in,
                                precision_t y_in,
                                precision_t z_in) {
  // The assigned number mainly follows from the iProc
  // i.e. 0 for left, 1 for front, 2 for right, 3 for back, 4 for below and 5 for top
  // The edge condition is a purely random choice
  // i.e. there are 8 corners and 6 surface, no perfect assignment
  if (z_in == 1)
    return 5;

  else if (y_in == -1 && x_in != 1)
    return 0;

  else if (x_in == 1 && y_in != 1)
    return 1;

  else if (y_in == 1 && x_in != -1)
    return 2;

  else if (x_in == -1 && y_in != -1)
    return 3;

  else if (z_in == -1)
    return 4;

  else {
    // The point is not on any of 6 surfaces of the a cube
    return -1;
  }
}

int64_t get_cube_surface_number(const arma_vec &point_in) {
  if (point_in.n_rows != 3) {
    // The input doesn't represent a point
    return -1;
  } else {
    return get_cube_surface_number(point_in[0],
                                   point_in[1],
                                   point_in[2]);
  }
}

// Helper variables / function ends. The following are all member
// functions of Grid class

// --------------------------------------------------------------------------
// Return the index of the last element that has altitude smaller than
// or equal to the input
// --------------------------------------------------------------------------

uint64_t Grid::search_altitude(const precision_t alt_in) const {
  // Copy from std::upper_bound. Can't directly use it
  // mainly because geoAlt_scgc(0, 0, *) can't be formed as an iterator
  uint64_t first, last, len;
  first = nGCs;
  last = nAlts - nGCs;
  len = last - first;

  while (len > 0) {
    uint64_t half = len >> 1;
    uint64_t mid = first + half;

    if (geoAlt_scgc(0, 0, mid) > alt_in)
      len = half;

    else {
      first = mid + 1;
      len = len - half - 1;
    }
  }

  return first - 1;
}

// --------------------------------------------------------------------------
// Return the index of the last element that has a value smaller than
//   or equal to the input
// - Optional argument (nGCs=0) since we cannot see grid info.
// --------------------------------------------------------------------------

// this replaces the above

uint64_t bisect_search_array(precision_t val_in, arma_vec ref_arr,
                             int64_t nGCs = 0) {
  uint64_t first, last, len;
  first = nGCs;
  last = ref_arr.size();
  len = last - first;

  while (len > 0) {
    uint64_t half = len >> 1;
    uint64_t mid = first + half;

    if (ref_arr(mid) > val_in)
      len = half;

    else {
      first = mid + 1;
      len = len - half - 1;
    }
  }

  return first - 1;
}

// --------------------------------------------------------------------------
// Get the range of a spherical grid
// --------------------------------------------------------------------------

void Grid::get_sphere_grid_range(struct sphere_range &sr) const {
  // Retrieve the range and delta of longitude, latitude and altitude
  sr.lon_min = geoLon_Corner(nGCs, nGCs, nGCs);
  sr.lon_max = geoLon_Corner(nLons - nGCs, nLats - nGCs, nAlts - nGCs);
  sr.lat_min = geoLat_Corner(nGCs, nGCs, nGCs);
  sr.lat_max = geoLat_Corner(nLons - nGCs, nLats - nGCs, nAlts - nGCs);
  // See init_geo_grid.cpp:443. The geoAlt_scgc doesn't add the coefficient 0.5
  // Use geoAlt_scgc instead of geoAlt_Corner
  sr.alt_min = geoAlt_scgc(nGCs, nGCs, nGCs);
  sr.alt_max = geoAlt_scgc(nLons - nGCs, nLats - nGCs, nAlts - nGCs);

  sr.dLon = geoLon_Corner(1, 0, 0) - geoLon_Corner(0, 0, 0);
  sr.dLat = geoLat_Corner(0, 1, 0) - geoLat_Corner(0, 0, 0);
}

// --------------------------------------------------------------------------
// Get the range of a cubesphere grid
// --------------------------------------------------------------------------

void Grid::get_cubesphere_grid_range(struct cubesphere_range &cr) const {
  // Get the location of the lower left corner, one step for row and one step for column
  arma_vec corner = sphere_to_cube(geoLon_Corner(nGCs, nGCs, nGCs),
                                   geoLat_Corner(nGCs, nGCs, nGCs));
  arma_vec step_row = sphere_to_cube(geoLon_Corner(nGCs + 1, nGCs, nGCs),
                                     geoLat_Corner(nGCs + 1, nGCs, nGCs));
  arma_vec step_col = sphere_to_cube(geoLon_Corner(nGCs, nGCs + 1, nGCs),
                                     geoLat_Corner(nGCs, nGCs + 1, nGCs));

  // Determine which axis the row expands along
  cr.row_direction = first_diff_index(corner, step_row);
  // Get the row_min and delta row;
  cr.row_min = corner[cr.row_direction];
  cr.drow = step_row[cr.row_direction] - cr.row_min;
  // Do the same for column
  cr.col_direction = first_diff_index(corner, step_col);
  cr.col_min = corner[cr.col_direction];
  cr.dcol = step_col[cr.col_direction] - cr.col_min;

  // Get the surface number
  arma_vec away_from_edge = corner;
  away_from_edge[cr.row_direction] += cr.drow;
  away_from_edge[cr.col_direction] += cr.dcol;
  cr.surface_number = get_cube_surface_number(away_from_edge);

  // Get the range of altitude, use geoAlt_scgc because the coefficient 0.5
  // is not added. See init_geo_grid.cpp:443 for detail
  cr.alt_min = geoAlt_scgc(nGCs, nGCs, nGCs);
  cr.alt_max = geoAlt_scgc(nLons - nGCs, nLats - nGCs, nAlts - nGCs);

  // Inclusion/exclusion part begins
  // The default settings are left-hand inclusive and right-hand exclusive
  // The surface number 0,1,2,3 always follow this rule
  cr.row_min_exclusive = false;
  cr.row_max_exclusive = true;
  cr.col_min_exclusive = false;
  cr.col_max_exclusive = true;

  if (cr.surface_number == 4) {
    // The bottom surface excludes all of its 4 edges, so when the min
    // equals -1, we need to turn exclusive to be true
    if (cr.row_min == -1)
      cr.row_min_exclusive = true;

    if (cr.col_min == -1)
      cr.col_min_exclusive = true;
  } else if (cr.surface_number == 5) {
    // The top surface includes all of its 4 edges, so when the max
    // equals 1 for row or -1 for col, we need to turn exclusive to be false
    if (cr.row_min + cr.drow * (nLons - 2 * nGCs) == 1)
      cr.row_max_exclusive = false;

    if (cr.col_min + cr.dcol * (nLons - 2 * nGCs) == -1)
      cr.col_max_exclusive = false;
  }
}


// --------------------------------------------------------------------------
// Get the range of a Dipole grid
// --------------------------------------------------------------------------

void Grid::get_dipole_grid_range(struct dipole_range &dr) const {
  // Retrieve the range and delta of longitude, latitude and altitude
  // ** Note the max/min are magnetic coordinates.  **
  dr.lon_min = i_corner_scgc(nGCs, nGCs, nGCs);
  dr.lon_max = i_corner_scgc(nLons - nGCs, nLats - nGCs, nAlts - nGCs);

  dr.lat_min = j_corner_scgc(nGCs, nGCs, nGCs);
  dr.lat_max = j_corner_scgc(nLons - nGCs, nLats - nGCs, nAlts - nGCs);

  // magAlt and geoAlt are the same, doesn't matter which we use:
  dr.alt_min = k_corner_scgc(nGCs, nGCs, nGCs);
  dr.alt_max = k_corner_scgc(nLons - nGCs, nLats - nGCs, nAlts - nGCs);

  // MagLon steps are uniform:
  dr.dLon = magLon_Corner(1, 0, 0) - magLon_Corner(0, 0, 0);
}


// --------------------------------------------------------------------------
// Set interpolation coefficients helper function for spherical grid
// Almost the copy of interp_sphere_linear_helper
// --------------------------------------------------------------------------


struct interp_coef_t Grid::get_interp_coef_sphere(const sphere_range &sr,
						  const precision_t lon_in,
						  const precision_t lat_in,
						  const precision_t alt_in) {

  // WARNING: IF WE ARE DEALING WITH LESS THAN THE WHOLE EARTH, THEN ALL THE POINTS WITH
  // LONGITUDE = geo_grid_input.lon_max = settings["GeoGrid"]["MaxLon"]
  // OR LATITUDE = geo_grid_input.lat_max = settings["GeoGrid"]["MaxLat"]
  // ARE EXCLUDED.
  // TO FIX IT, EACH GRID SHOULD BE ABLE TO ACCESS THE MaxLon and MaxLat

  // The structure which will be put into the interp_coefs.
  // Initialize in_grid to be false
  struct interp_coef_t coef;
  coef.in_grid = false;

  // Determine whether the point is inside this grid
  // Treat north pole specially because latitude is inclusive for
  //   both -cPI/2 and cPI/2
  // Don't check for altitude here!
  if (lon_in < sr.lon_min || lon_in >= sr.lon_max || lat_in < sr.lat_min
      || lat_in > sr.lat_max || (lat_in == sr.lat_max && sr.lat_max != cPI / 2)) {
    return coef;
  }

  // This point is in the grid!
  coef.in_grid = true;

  // ASSUMPTION: LONGITUDE AND LATITUDE ARE LINEARLY SPACED, nGCs >= 1
  // For the cell containing it, directly calculate its x and y index
  // Find its z index using binary search

  // The number of dLon between the innermost ghost cell and the given point
  coef.rRow = (lon_in - sr.lon_min) / sr.dLon + 0.5;
  // Take the integer part
  coef.iRow = static_cast<uint64_t>(coef.rRow);
  // Calculate the fractional part, which is the ratio for Longitude
  coef.rRow -= coef.iRow;
  // The actual x-axis index of the bottom-left of the cube used for interpolation
  coef.iRow += nGCs - 1;
  // Do the same for the Latitude
  coef.rCol = (lat_in - sr.lat_min) / sr.dLat + 0.5;
  coef.iCol = static_cast<uint64_t>(coef.rCol);
  coef.rCol -= coef.iCol;
  coef.iCol += nGCs - 1;

  // The altitude may not be linearly spaced, so use binary search to find
  // the first element smaller than or equal to the altitude of the give point
  // Implemented in search_altitude

  if (alt_in < sr.alt_min) {
    coef.iAlt = nGCs;
    coef.rAlt = alt_in - sr.alt_min;
    coef.below_grid = true;
    coef.above_grid = false;
  } else {
    if (alt_in > sr.alt_max) {
      coef.iAlt = nAlts - nGCs;
      coef.rAlt = alt_in - sr.alt_max;
      coef.below_grid = false;
      coef.above_grid = true;
    } else {
      coef.iAlt = bisect_search_array(alt_in,
				      geoAlt_scgc.tube(coef.iRow, coef.iCol),
				      nGCs);
      coef.rAlt =
	(alt_in - geoAlt_scgc(coef.iRow, coef.iCol, coef.iAlt))
	/ (geoAlt_scgc(coef.iRow, coef.iCol, coef.iAlt + 1) -
	   geoAlt_scgc(coef.iRow, coef.iCol, coef.iAlt));
      coef.below_grid = false;
      coef.above_grid = false;
    }
  }
  return coef;
}

// --------------------------------------------------------------------------
// Set interpolation coefficients helper function for cubesphere grid
// Almost the copy of interp_cubesphere_linear_helper
// --------------------------------------------------------------------------

struct interp_coef_t Grid::get_interp_coef_cubesphere(const cubesphere_range &cr,
                                      const precision_t lon_in,
                                      const precision_t lat_in,
                                      const precision_t alt_in) {
  // ASSUMPTION: THE SURFACES OF THE CUBE IS LINEARLY SPACED
  // I.E. init_geo_grid.cpp:106-137 WILL NEVER BE CHANGED

  // The structure which will be put into the interp_coefs. Initialize in_grid to be false
  struct interp_coef_t coef;
  coef.in_grid = false;

  // Find the projection point onto the cube and its surface number
  arma_vec point_in = sphere_to_cube(lon_in, lat_in);
  int64_t surface_in = get_cube_surface_number(point_in);

  // Determine whether the projection point is on the surface of the grid
  if (surface_in != cr.surface_number) {
    return coef;
  }

  // Calculate the theoretical fractional row index and column index
  precision_t row_frac_index, col_frac_index, row_in, col_in;
  row_in = point_in(cr.row_direction);
  col_in = point_in(cr.col_direction);
  row_frac_index = (row_in - cr.row_min) / cr.drow;
  col_frac_index = (col_in - cr.col_min) / cr.dcol;

  // Determine whether the projection point is out of range
  int64_t row_index_max, col_index_max;
  row_index_max = nLons - 2 * nGCs;
  col_index_max = nLats - 2 * nGCs;

  if (row_frac_index < 0 || (row_frac_index == 0 && cr.row_min_exclusive)
      || col_frac_index < 0 || (col_frac_index == 0 && cr.col_min_exclusive)
      || row_frac_index > row_index_max || (row_frac_index == row_index_max &&
                                            cr.row_max_exclusive)
      || col_frac_index > col_index_max || (col_frac_index == col_index_max &&
                                            cr.col_max_exclusive)) {
    return coef;
  }

  // This point is in the grid!
  coef.in_grid = true;

  // Get the real integer index and the interpolation coefficient
  uint64_t row_index, col_index, alt_index;
  precision_t rRow, rCol, rAlt;
  // Add 0.5 because the data we have is at the center of the cell rather than corner of the cell
  row_frac_index += 0.5;
  // Take the integer part
  coef.iRow = static_cast<uint64_t>(row_frac_index);
  // Calculate the fractional part, which is the coefficient
  coef.rRow = row_frac_index - coef.iRow;
  // The actual index considering the ghost cells
  coef.iRow += nGCs - 1;
  // Do the same for the column
  col_frac_index += 0.5;
  coef.iCol = static_cast<uint64_t>(col_frac_index);
  coef.rCol = col_frac_index - coef.iCol;
  coef.iCol += nGCs - 1;


  // The altitude may not be linearly spaced, so use binary search to find
  // the first element smaller than or equal to the altitude of the give point
  // Implemented in search_altitude

  if (alt_in < cr.alt_min) {
    coef.iAlt = nGCs;
    coef.rAlt = alt_in - cr.alt_min;
    coef.below_grid = true;
    coef.above_grid = false;
  } else {
    if (alt_in > cr.alt_max) {
      coef.iAlt = nAlts - nGCs;
      coef.rAlt = alt_in - cr.alt_max;
      coef.below_grid = false;
      coef.above_grid = true;
    } else {
      coef.iAlt = bisect_search_array(alt_in,
				      geoAlt_scgc.tube(coef.iRow, coef.iCol),
				      nGCs);
      coef.rAlt =
	(alt_in - geoAlt_scgc(coef.iRow, coef.iCol, coef.iAlt))
	/ (geoAlt_scgc(coef.iRow, coef.iCol, coef.iAlt + 1) -
	   geoAlt_scgc(coef.iRow, coef.iCol, coef.iAlt));
      coef.below_grid = false;
      coef.above_grid = false;
    }
  }
  return coef;
}


struct interp_coef_t Grid::get_interp_coef_dipole(const dipole_range &dr,
						  const precision_t lon_in,
						  const precision_t lat_in,
						  const precision_t alt_in) {

  // The structure which will be put into the interp_coefs. Initialize
  // in_grid to be false
  struct interp_coef_t coef;
  coef.in_grid = false;

  // Determine whether the point is inside this grid Treat north pole
  // specially because latitude is inclusive for both -cPI/2 and cPI/2
  if (lon_in < dr.lon_min ||
      lon_in >= dr.lon_max ||
      lat_in < dr.lat_min  ||
      lat_in > dr.lat_max ||
      (lat_in == dr.lat_max && dr.lat_max != cPI / 2)
      || alt_in < dr.alt_min || alt_in > dr.alt_max) {
    return coef;
  }

  // Put the coefficient into the vector
  coef.in_grid = true;
  
  // ASSUMPTION: LONGITUDE IS LINEARLY SPACED, nGCs >= 1
  // For the cell containing it, directly calculate its x index
  // Find y & z indices using a bisecting search

  // The number of dLon between the innermost ghost cell and the given point
  coef.rRow = (lon_in - dr.lon_min) / dr.dLon + 0.5;
  // Take the integer part
  coef.iRow = static_cast<uint64_t>(coef.rRow);
  // Calculate the fractional part, which is the ratio for Longitude
  coef.rRow -= coef.iRow;
  // The actual x-axis index of the bottom-left of the cube used for
  // interpolation
  coef.iRow += nGCs - 1;

  // Different from the sphere, latitude & altitude are not evenly spaced.
  // Use the bisect search function for both.

  // Lat needs to be done a little different because it could be increasing or
  // decreasing (depending on the hemisphere we're in). Take the absolute value!
  coef.iCol = bisect_search_array(abs(lat_in),
                                  abs(j_center_scgc.tube(coef.iRow, coef.iCol)), nGCs);

  // Use binary search to find the index for altitude
  if (alt_in < dr.alt_min) {
    coef.iAlt = nGCs;
    coef.rAlt = 0.0;
    coef.below_grid = true;
    coef.above_grid = false;
  } else {
    if (alt_in > dr.alt_max) {
      coef.iAlt = nAlts - nGCs;
      coef.rAlt = 0.0;
      coef.below_grid = false;
      coef.above_grid = true;
    } else {
      // Use binary search to find the index for altitude (handles
      // oblate planets)

      // need alt index to find lat coef
      coef.iAlt = bisect_search_array(alt_in,
				      k_center_scgc.tube(coef.iRow, coef.iCol),
				      nGCs);
      // then we can do the ratios:
      coef.rCol =
	(lat_in - magLat_scgc(coef.iRow, coef.iCol, coef.iAlt))
	/ (magLat_scgc(coef.iRow, coef.iCol + 1, coef.iAlt)
	   - magLat_scgc(coef.iRow, coef.iCol, coef.iAlt));
      coef.rAlt =
	(alt_in - geoAlt_scgc(coef.iRow, coef.iCol, coef.iAlt))
	/ (geoAlt_scgc(coef.iRow, coef.iCol, coef.iAlt + 1) -
	   geoAlt_scgc(coef.iRow, coef.iCol, coef.iAlt));
      coef.below_grid = false;
      coef.above_grid = false;
    }
  }

  return coef;
}


// --------------------------------------------------------------------------
// Set the interpolation coefficients
// --------------------------------------------------------------------------

bool Grid::set_interpolation_coefs(const std::vector<precision_t> &i_coords,
                                   const std::vector<precision_t> &j_coords,
                                   const std::vector<precision_t> &k_coords,
                                   bool areLocsGeo,// geo or mag?
                                   bool areLocsIJK // Are locs in 'native' coords?
                                  ) {
  /*
  Inputs:
    i_coord: longitude, either geo or mag (depends if areLocsGeo)
    j_coord:
      - Latitude if geographic
      - Invariant latitude if magnetic (areLocsGeo = false)
      - L-shell / dipole 'p': if magnetic and (areLocsIJK = false)
    k_coord:
      - Altitude/radius if NOT areLocsIJK
      - distance along field line, or dipole 'q' if areLocsIJK
  */

  std::string function = "Grid::set_interpolation_coefs";
  static int iFunction = -1;
  report.enter(function, iFunction);

  report.print(1, "interpolation gridtype : " + gridType);
  
  struct interp_coef_t coef;

  // If the size of Lons, Lats and Alts are not the same, return false
  if (i_coords.size() != j_coords.size() ||
      j_coords.size() != k_coords.size()) {
    report.error("Length of i,j,k vectors do not match!");
    return false;
  }

  // Clear the previous interpolation coefficients
  interp_coefs.clear();

  // ---------------------------------------------------
  // Cubesphere
  if (IsCubeSphereGrid) {
    // Calculate the range of the grid
    struct cubesphere_range cr;
    get_cubesphere_grid_range(cr);

    // Calculate the index and coefficients for each point
    for (size_t i = 0; i < i_coords.size(); ++i) {
      coef = get_interp_coef_cubesphere(cr,
					i_coords[i],
					j_coords[i],
					k_coords[i]);
      interp_coefs.push_back(coef);
    }

  }

  if (iGridShape_ == iSphere_) {
    report.print(1, "interpolation grid is sphere");

    struct sphere_range sr;
    get_sphere_grid_range(sr);

    // Calculate the index and coefficients for each point
    for (size_t i = 0; i < i_coords.size(); ++i) {
      coef = get_interp_coef_sphere(sr,
				    i_coords[i],
				    j_coords[i],
				    k_coords[i]);
      interp_coefs.push_back(coef);
    }
  }

  if (iGridShape_ == iDipole_) { // IsDipole
    report.print(1, "interpolation grid is dipole");

    // Calculate the range of the grid
    struct dipole_range dr;
    get_dipole_grid_range(dr);

    Planets planet;

    // make holders for dipole coordinates.
    int64_t iLoc, nPts = i_coords.size();
    std::vector<precision_t> mlon(nPts), p_coord(nPts), q_coord(nPts), dipijk(3);

    // these are the magnetic coordinates. A temporary step!  this is
    // a vector of cubes with shape (nPts, 1, 1) - avoids having to
    // overload things
    std::vector<arma_cube> magCoords;

    if (areLocsGeo) {
      arma_cube cubeCoord;
      cubeCoord = vec2cube(i_coords);
      magCoords.push_back(cubeCoord);
      cubeCoord = vec2cube(j_coords);
      magCoords.push_back(cubeCoord);
      cubeCoord = vec2cube(k_coords);
      magCoords.push_back(cubeCoord);

      magCoords = geo_to_mag(magCoords[0], magCoords[1], magCoords[2], planet);
      // for (iLoc = 0; iLoc < nPts; iLoc++) {
      //   magCoords = geo_to_mag(i_coords[iLoc], j_coords[iLoc], k_coords[iLoc], planet);
      //   mlon[iLoc] = magCoords[0];
      //   p_coord[iLoc] = magCoords[1];
      //   q_coord[iLoc] = magCoords[2];
    }

    else {
      magCoords = {vec2cube(i_coords), vec2cube(j_coords), vec2cube(k_coords)};
    }

    // std::vector<precision_t> dipcoords = geo_to_mag(i_coords[0], j_coords[0], k_coords[0], planet);
    std::vector<precision_t> dipCoords;

    if (!areLocsIJK) {
      std::vector<precision_t> planet_radii(nPts);

      for (iLoc = 0; iLoc < nPts; iLoc++) {
        // Convert from mag->dipole coordinates.
        if (areLocsGeo)  // we were given the geo-latitude
          planet_radii[iLoc] = planet.get_radius(j_coords[iLoc]);

        else {
          // equatorial radius :(
          planet_radii[iLoc] = planet.get_radius(0.0);
        }

        dipCoords = mag_to_ijk(i_coords[iLoc], j_coords[iLoc], k_coords[iLoc],
                               planet_radii[iLoc]);
        mlon[iLoc] = dipCoords[0];
        p_coord[iLoc] = dipCoords[1];
        q_coord[iLoc] = dipCoords[2];
      }
    } else {
      // just save the values
      for (iLoc = 0; iLoc < nPts; iLoc++) {
        mlon[iLoc] = i_coords[iLoc];
        p_coord[iLoc] = j_coords[iLoc];
        q_coord[iLoc] = k_coords[iLoc];
      }
    }

    // Calculate the index and coefficients for each point
    for (size_t i = 0; i < i_coords.size(); ++i) {
      coef = get_interp_coef_dipole(dr, mlon[i], p_coord[i], q_coord[i]);
      interp_coefs.push_back(coef);
    }
  }

  report.exit(function);
  return true;
}

// --------------------------------------------------------------------------
// Set the interpolation coefficients 
//    (v2 - return a list of interpolation coefficients)
// --------------------------------------------------------------------------

std::vector<struct interp_coef_t> Grid::get_interpolation_coefs(
                                    const std::vector<precision_t> &Lons,
                                    const std::vector<precision_t> &Lats,
                                    const std::vector<precision_t> &Alts) {

  int64_t nPts = Lons.size(), iPt;
  std::vector<struct interp_coef_t> listOfCoefs;
  struct interp_coef_t singleCoef;
  bool isBad = false;

  // If this is not a geo grid, return false
  if (!IsGeoGrid)
    isBad = true;

  // If the size of Lons, Lats and Alts are not the same, return false
  if (Lons.size() != Lats.size() || Lats.size() != Alts.size())
    isBad = true;

  if (isBad) {
    for (iPt = 0; iPt < nPts; ++iPt) {
      // Put the coefficient into the vector
      singleCoef.in_grid = false;
      listOfCoefs.push_back(singleCoef);
    }
    return listOfCoefs;
  }

  // Handle according to whether it is cubesphere or not
  if (IsCubeSphereGrid) {
    // Calculate the range of the grid
    struct cubesphere_range cr;
    get_cubesphere_grid_range(cr);

    // Calculate the index and coefficients for each point
    for (iPt = 0; iPt < nPts; ++iPt) {
      singleCoef = get_interp_coef_cubesphere(cr, Lons[iPt], Lats[iPt], Alts[iPt]);
      listOfCoefs.push_back(singleCoef);
    }
  } else {
    // Calculate the range of the grid
    struct sphere_range sr;
    get_sphere_grid_range(sr);

    // Calculate the index and coefficients for each point
    for (iPt = 0; iPt < nPts; ++iPt) {
      singleCoef = get_interp_coef_sphere(sr, Lons[iPt], Lats[iPt], Alts[iPt]);
      listOfCoefs.push_back(singleCoef);
    }      
  }

  return listOfCoefs;
}

// --------------------------------------------------------------------------
// Do the interpolation based on the coefficients stored in interp_coefs
// --------------------------------------------------------------------------

std::vector<precision_t> Grid::get_interpolation_values(
  const arma_cube &data) const {
  std::vector<precision_t> ans;

  // If the size of data is not the same as the size of grid, return an empty vector
  if (data.n_rows != nLons || data.n_cols != nLats || data.n_slices != nAlts)
    return ans;

  for (auto &it : interp_coefs) {
    // Do interpolation if in_grid = true. Push cNinf otherwise
    if (it.in_grid) {
      ans.push_back(interpolate_unit_cube(
                      data.subcube(it.iRow, it.iCol, it.iAlt, unit_cube_size),
                      it.rRow,
                      it.rCol,
                      it.rAlt
                    ));
      // Add std::cout if needed here
      // std::cout << "iProc = " << iProc << " interpolates the point successfully\n";
    } else
      ans.push_back(cNinf);
  }

  return ans;
}

// --------------------------------------------------------------------------
// Do the interpolation based on the coefficients passed in
// --------------------------------------------------------------------------

std::vector<precision_t> Grid::get_interpolation_values(arma_cube data,
                                                        std::vector<struct interp_coef_t> coefArray ) {
  std::vector<precision_t> ans;

  // If the size of data is not the same as the size of grid, return an empty vector
  if (data.n_rows != nLons || data.n_cols != nLats || data.n_slices != nAlts)
    return ans;

  for (auto &it : coefArray) {
    // Do interpolation if in_grid = true. Push cNinf otherwise
    if (it.in_grid) {
      ans.push_back(interpolate_unit_cube(
                      data.subcube(it.iRow, it.iCol, it.iAlt, unit_cube_size),
                      it.rRow,
                      it.rCol,
                      it.rAlt
                    ));
      // Add std::cout if needed here
      // std::cout << "iProc = " << iProc << " interpolates the point successfully\n";
    } else
      ans.push_back(cNinf);
  }

  return ans;
}

