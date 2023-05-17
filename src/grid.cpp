// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>

#include "aether.h"

// --------------------------------------------------------------------------
// Initialize Grid class
// --------------------------------------------------------------------------

Grid::Grid(int nX_in, int nY_in, int nZ_in, int nGCs_in) {

  nX = nX_in + nGCs_in * 2;
  nLons = nX;
  nY = nY_in + nGCs_in * 2;
  nLats = nY;
  nZ = nZ_in + nGCs_in * 2;
  nAlts = nZ;
  nGCs = nGCs_in;

  geoLon_scgc.set_size(nX, nY, nZ);
  geoLat_scgc.set_size(nX, nY, nZ);
  geoAlt_scgc.set_size(nX, nY, nZ);
  geoLocalTime_scgc.set_size(nX, nY, nZ);


  geoLon_Left.set_size(nX + 1, nY, nZ);
  geoLat_Left.set_size(nX + 1, nY, nZ);

  geoLon_Down.set_size(nX, nY + 1, nZ);
  geoLat_Down.set_size(nX, nY + 1, nZ);

  geoLon_Corner.set_size(nX + 1, nY + 1, nZ + 1);
  geoLat_Corner.set_size(nX + 1, nY + 1, nZ + 1);
  geoAlt_Corner.set_size(nX + 1, nY + 1, nZ + 1);

  geoAlt_Below.set_size(nX, nY, nZ + 1);


  geoX_scgc.set_size(nX, nY, nZ);
  geoY_scgc.set_size(nX, nY, nZ);
  geoZ_scgc.set_size(nX, nY, nZ);

  magLon_scgc.set_size(nX, nY, nZ);
  magLat_scgc.set_size(nX, nY, nZ);
  magAlt_scgc.set_size(nX, nY, nZ);

  magX_scgc.set_size(nX, nY, nZ);
  magY_scgc.set_size(nX, nY, nZ);
  magZ_scgc.set_size(nX, nY, nZ);

  magLocalTime_scgc.set_size(nX, nY, nZ);

  radius_scgc.set_size(nX, nY, nZ);
  radius2_scgc.set_size(nX, nY, nZ);
  radius2i_scgc.set_size(nX, nY, nZ);

  dalt_center_scgc.set_size(nX, nY, nZ);
  dalt_lower_scgc.set_size(nX, nY, nZ);
  dalt_ratio_scgc.set_size(nX, nY, nZ);
  dalt_ratio_sq_scgc.set_size(nX, nY, nZ);

  dlat_center_scgc.set_size(nX, nY, nZ);
  dlat_center_dist_scgc.set_size(nX, nY, nZ);

  dlon_center_scgc.set_size(nX, nY, nZ);
  dlon_center_dist_scgc.set_size(nX, nY, nZ);

  sza_scgc.set_size(nX, nY, nZ);
  cos_sza_scgc.set_size(nX, nY, nZ);

  bfield_vcgc = make_cube_vector(nX, nY, nZ, 3);
  bfield_unit_vcgc = make_cube_vector(nX, nY, nZ, 3);
  bfield_mag_scgc.set_size(nX, nY, nZ);
  bfield_mag_scgc.zeros();

  GSE_XYZ_vcgc = make_cube_vector(nX, nY, nZ, 3);

  mag_pole_north_ll.set_size(2);
  mag_pole_south_ll.set_size(2);
  mag_pole_north_ll.zeros();
  mag_pole_south_ll.zeros();

  arma_cube tmp_col(1, 1, nZ);
  mag_pole_north_gse.push_back(tmp_col);
  mag_pole_north_gse.push_back(tmp_col);
  mag_pole_north_gse.push_back(tmp_col);

  mag_pole_south_gse.push_back(tmp_col);
  mag_pole_south_gse.push_back(tmp_col);
  mag_pole_south_gse.push_back(tmp_col);

  HasBField = 0;
}

// --------------------------------------------------------------------------
// write restart out files for the grid
// --------------------------------------------------------------------------

bool Grid::write_restart(std::string dir) {
  bool DidWork = true;

  // All Ensemble member grids should be the same, so only need to write
  // out the 0th member
  if (iMember == 0) {
    try {
      OutputContainer RestartContainer;
      RestartContainer.set_netcdf();
      RestartContainer.set_directory(dir);
      RestartContainer.set_version(0.1);
      RestartContainer.set_time(0.0);

      // Output Cell Centers
      RestartContainer.set_filename("grid_" + cGrid);
      RestartContainer.store_variable(longitude_name,
                                      longitude_unit,
                                      geoLon_scgc);
      RestartContainer.store_variable(latitude_name,
                                      latitude_unit,
                                      geoLat_scgc);
      RestartContainer.store_variable(altitude_name,
                                      altitude_unit,
                                      geoAlt_scgc);
      RestartContainer.write();
      RestartContainer.clear_variables();

      // Output Corners
      RestartContainer.set_filename("grid_corners_" + cGrid);
      RestartContainer.store_variable(longitude_name + " Corners",
                                      longitude_unit,
                                      geoLon_Corner);
      RestartContainer.store_variable(latitude_name + " Corners",
                                      latitude_unit,
                                      geoLat_Corner);
      RestartContainer.store_variable(altitude_name + " Corners",
                                      altitude_unit,
                                      geoAlt_Corner);
      RestartContainer.write();
      RestartContainer.clear_variables();

      // Output Left Sides
      RestartContainer.set_filename("grid_left_" + cGrid);
      RestartContainer.store_variable(longitude_name + " Left",
                                      longitude_unit,
                                      geoLon_Left);
      RestartContainer.store_variable(latitude_name + " Left",
                                      latitude_unit,
                                      geoLat_Left);
      RestartContainer.write();
      RestartContainer.clear_variables();

      // Output Down Sides
      RestartContainer.set_filename("grid_down_" + cGrid);
      RestartContainer.store_variable(longitude_name + " Down",
                                      longitude_unit,
                                      geoLon_Down);
      RestartContainer.store_variable(latitude_name + " Down",
                                      latitude_unit,
                                      geoLat_Down);
      RestartContainer.write();
      RestartContainer.clear_variables();

      // Output Below
      RestartContainer.set_filename("grid_below_" + cGrid);
      RestartContainer.store_variable(altitude_name + " Below",
                                      altitude_unit,
                                      geoAlt_Below);

      RestartContainer.write();
      RestartContainer.clear_variables();

    } catch (...) {
      std::cout << "Error writing grid restart file!\n";
      DidWork = false;
    }
  }

  DidWork = sync_across_all_procs(DidWork);
  return DidWork;
}

// --------------------------------------------------------------------------
// read restart out files for the grid
// - Returns true if everything worked ok
// --------------------------------------------------------------------------

bool Grid::read_restart(std::string dir) {

  bool DidWork = true;
  int64_t iVar;

  // While only the 0th ensemble member writes, all ensemble members have
  // to read the grid.

  try {
    OutputContainer RestartContainer;
    RestartContainer.set_netcdf();
    RestartContainer.set_directory(dir);
    RestartContainer.set_version(0.1);
    // Cell Centers:
    RestartContainer.set_filename("grid_" + cGrid);
    RestartContainer.read_container_netcdf();
    geoLon_scgc = RestartContainer.get_element_value(longitude_name);
    geoLat_scgc = RestartContainer.get_element_value(latitude_name);
    geoAlt_scgc = RestartContainer.get_element_value(altitude_name);
    // Down Edges:
    RestartContainer.set_filename("grid_below_" + cGrid);
    RestartContainer.read_container_netcdf();
    geoAlt_Below = RestartContainer.get_element_value(altitude_name +
                                                      " Below");
    // Cell Corners:
    RestartContainer.set_filename("grid_corners_" + cGrid);
    RestartContainer.read_container_netcdf();
    geoLon_Corner = RestartContainer.get_element_value(longitude_name +
                                                       " Corners");
    geoLat_Corner = RestartContainer.get_element_value(latitude_name +
                                                       " Corners");
    geoAlt_Corner = RestartContainer.get_element_value(altitude_name +
                                                       " Corners");
    // Left Edges:
    RestartContainer.set_filename("grid_left_" + cGrid);
    RestartContainer.read_container_netcdf();
    geoLon_Left = RestartContainer.get_element_value(longitude_name +
                                                     " Left");
    geoLat_Left = RestartContainer.get_element_value(latitude_name +
                                                     " Left");
    // Down Edges:
    RestartContainer.set_filename("grid_down_" + cGrid);
    RestartContainer.read_container_netcdf();
    geoLon_Down = RestartContainer.get_element_value(longitude_name +
                                                     " Down");
    geoLat_Down = RestartContainer.get_element_value(latitude_name +
                                                     " Down");

  } catch (...) {
    std::cout << "Error reading grid restart file!\n";
    DidWork = false;
  }

  DidWork = sync_across_all_procs(DidWork);
  return DidWork;
}

// --------------------------------------------------------------------------
// Report Grid boundaries:
// --------------------------------------------------------------------------

void Grid::report_grid_boundaries() {
  std::cout << "---------------------------------------------------\n";
  std::cout << "Grid Boundaries (min / max):\n";
  std::cout << "Lon : "
            << geoLon_scgc.min() << " / "
            << geoLon_scgc.max() << "\n";
  std::cout << "Lat : "
            << geoLat_scgc.min() << " / "
            << geoLat_scgc.max() << "\n";
  std::cout << "Alt : "
            << geoAlt_scgc.min() << " / "
            << geoAlt_scgc.max() << "\n";
  std::cout << "---------------------------------------------------\n";
}

// --------------------------------------------------------------------------
// Get whether the grid is a geographic grid (or magnetic - return 0)
// --------------------------------------------------------------------------

int Grid::get_IsGeoGrid() {
  return IsGeoGrid;
}

// --------------------------------------------------------------------------
// Get whether the grid is a geographic grid (or magnetic - return 0)
// --------------------------------------------------------------------------

bool Grid::get_HasBField() {
  return HasBField;
}

// --------------------------------------------------------------------------
// Set whether the grid is a geographic grid (or magnetic - set to 0)
// --------------------------------------------------------------------------

void Grid::set_IsGeoGrid(int value) {
  IsGeoGrid = value;
}

// --------------------------------------------------------------------------
// Get total number of grid points
// --------------------------------------------------------------------------

int64_t Grid::get_nPointsInGrid() {
  int64_t nPoints;
  nPoints = int64_t(nX) * int64_t(nY) * int64_t(nZ);
  return nPoints;
}

// --------------------------------------------------------------------------
// Get a bunch of sizes within the grid
// --------------------------------------------------------------------------

int64_t Grid::get_nX() {
  return nX;
}
int64_t Grid::get_nY() {
  return nY;
}
int64_t Grid::get_nZ() {
  return nZ;
}

int64_t Grid::get_nLons() {
  return nLons;
}
int64_t Grid::get_nLats() {
  return nLats;
}
int64_t Grid::get_nAlts() {
  return nAlts;
}

int64_t Grid::get_nGCs() {
  return nGCs;
}

// The size of a 2*2*2 arma cube
const arma::SizeCube Grid::unit_cube_size = arma::size(2, 2, 2);

// --------------------------------------------------------------------------
// Return the first index of two vectors on which they have different values
// --------------------------------------------------------------------------

int64_t Grid::first_diff_index(const arma_vec &a, const arma_vec &b) {
    int64_t i;
    for (i = 0; i < std::min(a.n_rows, b.n_rows); ++i) {
        if (std::abs(a[i] - b[i]) > cSmall) {
            return i;
        }
    }
    return i;
}

// --------------------------------------------------------------------------
// Return the index of the last element that has altitude smaller than or euqal to the input
// --------------------------------------------------------------------------

uint64_t Grid::search_altitude(const precision_t alt_in) {
    // Copy from std::upper_bound. Can't directly use it
    // mainly because geoAlt_scgc(0, 0, *) can't be formed as an iterator
    uint64_t first, last, len;
    first = nGCs;
    last = nAlts - nGCs;
    len = last - first;
    while (len > 0) {
        uint64_t half = len >> 1;
        uint64_t mid = first + half;
        if (geoAlt_scgc(0, 0, mid) > alt_in) {
            len = half;
        } else {
            first = mid + 1;
            len = len - half - 1;
        }
    }
    return first - 1;
}

// --------------------------------------------------------------------------
// Project a point described by lon and lat to a point on a surface of the 2-2-2 cube
// --------------------------------------------------------------------------

arma_vec Grid::sphere_to_cube(precision_t lon_in, precision_t lat_in) {
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
        if (std::abs(ans[i] + 1) < cSmall) {
            ans[i] = -1;
        } else if (std::abs(ans[i] - 1) < cSmall) {
            ans[i] = 1;
        }
    }

    return ans;
}

// --------------------------------------------------------------------------
// Assign any point on the surface of a cube a number within [0,5]
// --------------------------------------------------------------------------

int64_t Grid::get_cube_surface_number(precision_t x_in,
                                      precision_t y_in,
                                      precision_t z_in) {
    // The assigned number mainly follows from the iProc
    // i.e. 0 for left, 1 for front, 2 for right, 3 for back, 4 for below and 5 for top
    // The edge condition is a purely random choice
    // i.e. there are 8 corners and 6 surface, no perfect assignment
    if (z_in == 1) {
        return 5;
    } else if (y_in == -1 && x_in != 1) {
        return 0;
    } else if (x_in == 1 && y_in != 1) {
        return 1;
    } else if (y_in == 1 && x_in != -1) {
        return 2;
    } else if (x_in == -1 && y_in != -1) {
        return 3;
    } else if (z_in == -1) {
        return 4;
    } else {
        // The point is not on any of 6 surfaces of the a cube
        return -1;
    }
}

int64_t Grid::get_cube_surface_number(const arma_vec &point_in) {
    if (point_in.n_rows != 3) {
        // The input doesn't represent a point
        return -1;
    } else {
        return get_cube_surface_number(point_in[0],
                                       point_in[1],
                                       point_in[2]);
    }
}

// --------------------------------------------------------------------------
// Get the range of a spherical grid
// --------------------------------------------------------------------------

void Grid::get_sphere_grid_range(struct sphere_range &sr) {
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

void Grid::get_cubesphere_grid_range(struct cubesphere_range &cr) {
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
        if (cr.row_min == -1) {
            cr.row_min_exclusive = true;
        }
        if (cr.col_min == -1) {
            cr.col_min_exclusive = true;
        }
    } else if (cr.surface_number == 5) {
        // The top surface includes all of its 4 edges, so when the max
        // equals 1 for row or -1 for col, we need to turn exclusive to be false
        if (cr.row_min + cr.drow * (nLons - 2 * nGCs) == 1) {
            cr.row_max_exclusive = false;
        }
        if (cr.col_min + cr.dcol * (nLons - 2 * nGCs) == -1) {
            cr.col_max_exclusive = false;
        }
    }
}

// --------------------------------------------------------------------------
// Linear interpolation helper function for spherical grid
// --------------------------------------------------------------------------

precision_t Grid::interp_sphere_linear_helper(const arma_cube &data,
                                              const sphere_range &sr,
                                              const precision_t lon_in,
                                              const precision_t lat_in,
                                              const precision_t alt_in) {
    // WARNING: IF WE ARE DEALING WITH LESS THAN THE WHOLE EARTH, THEN ALL THE POINTS WITH
    // LONGITUDE = geo_grid_input.lon_max = settings["GeoGrid"]["MaxLon"]
    // OR LATITUDE = geo_grid_input.lat_max = settings["GeoGrid"]["MaxLat"]
    // ARE EXCLUDED.
    // TO FIX IT, EACH GRID SHOULD BE ABLE TO ACCESS THE MaxLon and MaxLat

    // Determine whether the point is inside this grid
    // Treat north pole specially because latitude is inclusive for both -cPI/2 and cPI/2
    if (lon_in < sr.lon_min || lon_in >= sr.lon_max || lat_in < sr.lat_min
        || lat_in > sr.lat_max || (lat_in == sr.lat_max && sr.lat_max != cPI/2)
        || alt_in < sr.alt_min || alt_in > sr.alt_max) {
        return cNinf;
    }

    // ASSUMPTION: LONGITUDE AND LATITUDE ARE LINEARLY SPACED, nGCs >= 1
    // For the cell containing it, directly calculate its x and y index
    // Find its z index using binary search
    uint64_t iLon, iLat, iAlt;
    precision_t rLon, rLat, rAlt;

    // The number of dLon between the innermost ghost cell and the given point
    rLon = (lon_in - sr.lon_min) / sr.dLon + 0.5;
    // Take the integer part
    iLon = static_cast<uint64_t>(rLon);
    // Calculate the fractional part, which is the ratio for Longitude
    rLon -= iLon;
    // The actual x-axis index of the bottom-left of the cube used for interpolation
    iLon += nGCs - 1;
    // Do the same for the Latitude
    rLat = (lat_in - sr.lat_min) / sr.dLat + 0.5;
    iLat = static_cast<uint64_t>(rLat);
    rLat -= iLat;
    iLat += nGCs - 1;

    // The altitude may not be linearly spaced, so use binary search to find
    // the first element smaller than or equal to the altitude of the give point
    // Implemented in search_altitude
    iAlt = search_altitude(alt_in);
    rAlt = (alt_in - geoAlt_scgc(0, 0, iAlt))
            / (geoAlt_scgc(0, 0, iAlt + 1) - geoAlt_scgc(0, 0, iAlt));

    std::cout << "iProc = " << iProc << " interpolates point ("
              << lon_in << ", " << lat_in << ", " << alt_in << ") successfully\n";

    // Return the estimated value
    return interpolate_unit_cube(data.subcube(iLon, iLat, iAlt, unit_cube_size),
                                 rLon,
                                 rLat,
                                 rAlt);
}

// --------------------------------------------------------------------------
// Linear interpolation helper function for spherical grid
// --------------------------------------------------------------------------

precision_t Grid::interp_cubesphere_linear_helper(const arma_cube &data,
                                                  const cubesphere_range &cr,
                                                  const precision_t lon_in,
                                                  const precision_t lat_in,
                                                  const precision_t alt_in) {
    // ASSUMPTION: THE SURFACES OF THE CUBE IS LINEARLY SPACED
    // I.E. init_geo_grid.cpp:106-137 WILL NEVER BE CHANGED

    // Find the projection point onto the cube and its surface number
    arma_vec point_in = sphere_to_cube(lon_in, lat_in);
    int64_t surface_in = get_cube_surface_number(point_in);

    // If the projection point is not on the surface of the grid, return cNinf
    if (surface_in != cr.surface_number) {
        return cNinf;
    }

    // Calculate the theoretical fractional row index and column index
    precision_t row_frac_index, col_frac_index, row_in, col_in;
    row_in = point_in(cr.row_direction);
    col_in = point_in(cr.col_direction);
    row_frac_index = (row_in - cr.row_min) / cr.drow;
    col_frac_index = (col_in - cr.col_min) / cr.dcol;

    // Return cNinf if it is out of range
    int64_t row_index_max, col_index_max;
    row_index_max = nLons - 2 * nGCs;
    col_index_max = nLats - 2 * nGCs;
    if (row_frac_index < 0 || (row_frac_index == 0 && cr.row_min_exclusive)
     || col_frac_index < 0 || (col_frac_index == 0 && cr.col_min_exclusive)
     || row_frac_index > row_index_max || (row_frac_index == row_index_max && cr.row_max_exclusive)
     || col_frac_index > col_index_max || (col_frac_index == col_index_max && cr.col_max_exclusive)
     || alt_in < cr.alt_min || alt_in > cr.alt_max) {
        return cNinf;
    }

    // Get the real integer index and the interpolation coefficient
    uint64_t row_index, col_index, alt_index;
    precision_t rRow, rCol, rAlt;
    // Add 0.5 because the data we have is at the center of the cell rather than corner of the cell
    row_frac_index += 0.5;
    // Take the integer part
    row_index = static_cast<uint64_t>(row_frac_index);
    // Calculate the fractional part, which is the coefficient
    rRow = row_frac_index - row_index;
    // The actual index considering the ghost cells
    row_index += nGCs - 1;
    // Do the same for the column
    col_frac_index += 0.5;
    col_index = static_cast<uint64_t>(col_frac_index);
    rCol = col_frac_index - col_index;
    col_index += nGCs - 1;
    // Use binary search to find the index for altitude
    alt_index = search_altitude(alt_in);
    rAlt = (alt_in - geoAlt_scgc(0, 0, alt_index))
            / (geoAlt_scgc(0, 0, alt_index + 1) - geoAlt_scgc(0, 0, alt_index));

    std::cout << "iProc = " << iProc << " interpolates point ("
              << lon_in << ", " << lat_in << ", " << alt_in << ") successfully\n";

    // Return the estimated value
    return interpolate_unit_cube(data.subcube(row_index, col_index, alt_index, unit_cube_size),
                                 rRow,
                                 rCol,
                                 rAlt);
}

// --------------------------------------------------------------------------
// Estimate the value of the point at (lon_in, lat_in, alt_in)
// --------------------------------------------------------------------------

precision_t Grid::interp_linear(const arma_cube &data,
                                const precision_t lon_in,
                                const precision_t lat_in,
                                const precision_t alt_in) {
    // Check that the size of the data is the same as the size of the grid
    if (data.n_rows != nLons || data.n_cols != nLats || data.n_slices != nAlts) {
        return cNinf;
    }

    // Check whether the grid is sphere or cubesphere
    if (IsCubeSphereGrid) {
        struct cubesphere_range cr;
        get_cubesphere_grid_range(cr);
        std::vector<precision_t> ans;
        return interp_cubesphere_linear_helper(data, cr, lon_in, lat_in, alt_in);
    } else {
        struct sphere_range sr;
        get_sphere_grid_range(sr);
        return interp_sphere_linear_helper(data, sr, lon_in, lat_in, alt_in);
    }
}

// --------------------------------------------------------------------------
// Estimate the value of points specified by three vectors of lon, lat and alt
// --------------------------------------------------------------------------

std::vector<precision_t> Grid::interp_linear(const arma_cube &data,
                                             const std::vector<precision_t> &Lons,
                                             const std::vector<precision_t> &Lats,
                                             const std::vector<precision_t> &Alts) {
    // If the size of Lons, Lats and Alts are not the same, return a vector
    // with only one element {cNinf}
    if (Lons.size() != Lats.size() || Lats.size() != Alts.size()) {
        return std::vector<precision_t>(1, cNinf);
    }
    // If the size of data is not the same as the size of grid, return a vector
    // with as many cNinf as the number of points
    if (data.n_rows != nLons || data.n_cols != nLats || data.n_slices != nAlts) {
        return std::vector<precision_t>(Alts.size(), cNinf);
    }

    // Check whether the grid is sphere or cubesphere
    std::vector<precision_t> ans;
    if (IsCubeSphereGrid) {
        struct cubesphere_range cr;
        get_cubesphere_grid_range(cr);
        for (int64_t i = 0; i < Lons.size(); ++i) {
            ans.push_back(interp_cubesphere_linear_helper(data, cr, Lons[i], Lats[i], Alts[i]));
        }
    } else {
        struct sphere_range sr;
        get_sphere_grid_range(sr);
        for (int64_t i = 0; i < Lons.size(); ++i) {
            ans.push_back(interp_sphere_linear_helper(data, sr, Lons[i], Lats[i], Alts[i]));
        }
    }

    return ans;
}
