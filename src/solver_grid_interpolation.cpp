// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "aether.h"

// Hepler varialbes / function begins. These are only used inside this cpp file and neither declared nor visible in any other file

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
// Project a point described by lon and lat to a point on a surface of the 2-2-2 cube
// --------------------------------------------------------------------------

arma_vec sphere_to_cube(precision_t lon_in, precision_t lat_in) {
  // See init_geo_grid.cpp:126. The offset for lon is subtracted
  lon_in = lon_in - 3 * cPI / 4;

  // Transfer polar coordinate to cartesian coordinate
  precision_t xy_temp;
  arma_vec ans(3);
  ans[2] = sin(lat_in);
  xy_temp = cos(lat_in);
  ans[1] = xy_temp * sin(lon_in);
  ans[0] = xy_temp * cos(lon_in);

  // Project this point onto the surface of cube
  precision_t coef = 1.0 / std::max({std::abs(ans[0]), std::abs(ans[1]), std::abs(ans[2])});
  ans *= coef;

  // Round the number if it is close to 1 or -1, otherwise the == and != operator
  // won't behave as expected because of the accuracy problem of floating point numbers
  for (int64_t i = 0; i < 3; ++i) {
    if (std::abs(ans[i] + 1) < cSmall)
      ans[i] = -1;

    else if (std::abs(ans[i] - 1) < cSmall)
      ans[i] = 1;
  }

  return ans;
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

// Helper variables / function ends. The following are all member functions of Grid class

// --------------------------------------------------------------------------
// Return the index of the last element that has altitude smaller than or euqal to the input
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
// Set interpolation coefficients helper function for spherical grid
// Almost the copy of interp_sphere_linear_helper
// --------------------------------------------------------------------------

void Grid::set_interp_coef_sphere(const sphere_range &sr,
                                  const precision_t lon_in,
                                  const precision_t lat_in,
                                  const precision_t alt_in) {
  // WARNING: IF WE ARE DEALING WITH LESS THAN THE WHOLE EARTH, THEN ALL THE POINTS WITH
  // LONGITUDE = geo_grid_input.lon_max = settings["GeoGrid"]["MaxLon"]
  // OR LATITUDE = geo_grid_input.lat_max = settings["GeoGrid"]["MaxLat"]
  // ARE EXCLUDED.
  // TO FIX IT, EACH GRID SHOULD BE ABLE TO ACCESS THE MaxLon and MaxLat

  // The structure which will be put into the interp_coefs. Initialize in_grid to be false
  struct interp_coef_t coef;
  coef.in_grid = false;

  // Determine whether the point is inside this grid
  // Treat north pole specially because latitude is inclusive for both -cPI/2 and cPI/2
  if (lon_in < sr.lon_min || lon_in >= sr.lon_max || lat_in < sr.lat_min
      || lat_in > sr.lat_max || (lat_in == sr.lat_max && sr.lat_max != cPI / 2)
      || alt_in < sr.alt_min || alt_in > sr.alt_max) {
    interp_coefs.push_back(coef);
    return;
  }

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
  coef.iAlt = search_altitude(alt_in);
  coef.rAlt = (alt_in - geoAlt_scgc(0, 0, coef.iAlt))
              / (geoAlt_scgc(0, 0, coef.iAlt + 1) - geoAlt_scgc(0, 0, coef.iAlt));

  // Put the coefficient into the vector
  coef.in_grid = true;
  interp_coefs.push_back(coef);
}

// --------------------------------------------------------------------------
// Set interpolation coefficients helper function for cubesphere grid
// Almost the copy of interp_cubesphere_linear_helper
// --------------------------------------------------------------------------

void Grid::set_interp_coef_cubesphere(const cubesphere_range &cr,
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
    interp_coefs.push_back(coef);
    return;
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
                                            cr.col_max_exclusive)
      || alt_in < cr.alt_min || alt_in > cr.alt_max) {
    interp_coefs.push_back(coef);
    return;
  }

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
  // Use binary search to find the index for altitude
  coef.iAlt = search_altitude(alt_in);
  coef.rAlt = (alt_in - geoAlt_scgc(0, 0, coef.iAlt))
              / (geoAlt_scgc(0, 0, coef.iAlt + 1) - geoAlt_scgc(0, 0, coef.iAlt));

  // Put the coefficient into the vector
  coef.in_grid = true;
  interp_coefs.push_back(coef);
}

// --------------------------------------------------------------------------
// Set the interpolation coefficients
// --------------------------------------------------------------------------

bool Grid::set_interpolation_coefs(const std::vector<precision_t> &Lons,
                                   const std::vector<precision_t> &Lats,
                                   const std::vector<precision_t> &Alts) {
  // If this is not a geo grid, return false
  if (!IsGeoGrid)
    return false;

  // If the size of Lons, Lats and Alts are not the same, return false
  if (Lons.size() != Lats.size() || Lats.size() != Alts.size())
    return false;

  // Clear the previous interpolation coefficients
  interp_coefs.clear();

  // Handle according to whether it is cubesphere or not
  if (IsCubeSphereGrid) {
    // Calculate the range of the grid
    struct cubesphere_range cr;
    get_cubesphere_grid_range(cr);

    // Calculate the index and coefficients for each point
    for (size_t i = 0; i < Lons.size(); ++i)
      set_interp_coef_cubesphere(cr, Lons[i], Lats[i], Alts[i]);
  } else {
    // Calculate the range of the grid
    struct sphere_range sr;
    get_sphere_grid_range(sr);

    // Calculate the index and coefficients for each point
    for (size_t i = 0; i < Lons.size(); ++i)
      set_interp_coef_sphere(sr, Lons[i], Lats[i], Alts[i]);
  }

  return true;
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
