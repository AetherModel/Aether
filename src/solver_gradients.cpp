// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "../include/aether.h"

// --------------------------------------------------------------------------
// Calculate the gradient in the longitudinal direction
// --------------------------------------------------------------------------

std::vector<arma_cube> calc_gradient_vector(arma_cube value_scgc, Grid grid) {

  std::vector<arma_cube> gradient_vcgc;

  gradient_vcgc.push_back(calc_gradient_lon(value_scgc, grid));
  gradient_vcgc.push_back(calc_gradient_lat(value_scgc, grid));
  gradient_vcgc.push_back(calc_gradient_alt(value_scgc, grid));
  return gradient_vcgc;
}

// --------------------------------------------------------------------------
// Calculate the gradient in the longitudinal direction
// --------------------------------------------------------------------------

arma_cube calc_gradient_lon(arma_cube value, Grid grid) {

  int64_t nLons = grid.get_nLons();
  int64_t nLats = grid.get_nLats();
  int64_t nAlts = grid.get_nAlts();
  int64_t iLon;

  arma_cube gradient(nLons, nLats, nAlts);
  gradient.zeros();

  // Interior:
  for (iLon = 1; iLon < nLons - 1; iLon++)
    gradient.row(iLon) =
      (value.row(iLon + 1) - value.row(iLon - 1)) /
      (2 * grid.dlon_center_dist_scgc.row(iLon));

  // Lower (one sided):
  iLon = 0;
  gradient.row(iLon) =
    (value.row(iLon + 1) - value.row(iLon)) /
    grid.dlon_center_dist_scgc.row(iLon);

  // Upper (one sided):
  iLon = nLons - 1;
  gradient.row(iLon) =
    (value.row(iLon) - value.row(iLon - 1)) /
    grid.dlon_center_dist_scgc.row(iLon);

  return gradient;
}

// --------------------------------------------------------------------------
// Calculate the gradient in the latitudinal direction
// --------------------------------------------------------------------------

arma_cube calc_gradient_lat(arma_cube value, Grid grid) {

  int64_t nLons = grid.get_nLons();
  int64_t nLats = grid.get_nLats();
  int64_t nAlts = grid.get_nAlts();
  int64_t iLat;

  arma_cube gradient(nLons, nLats, nAlts);
  gradient.zeros();

  // Interior:
  for (iLat = 1; iLat < nLats - 1; iLat++)
    gradient.col(iLat) =
      (value.col(iLat + 1) - value.col(iLat - 1)) /
      (2 * grid.dlat_center_dist_scgc.col(iLat));

  // Lower (one sided):
  iLat = 0;
  gradient.col(iLat) =
    (value.col(iLat + 1) - value.col(iLat)) /
    grid.dlat_center_dist_scgc.col(iLat);

  // Upper (one sided):
  iLat = nLats - 1;
  gradient.col(iLat) =
    (value.col(iLat) - value.col(iLat - 1)) /
    grid.dlat_center_dist_scgc.col(iLat);

  return gradient;
}

// --------------------------------------------------------------------------
// Calculate the gradient in the altitudinal direction
// --------------------------------------------------------------------------

arma_cube calc_gradient_alt(arma_cube value, Grid grid) {

  int64_t nLons = grid.get_nLons();
  int64_t nLats = grid.get_nLats();
  int64_t nAlts = grid.get_nAlts();
  int64_t nGCs = grid.get_nGCs();
  int64_t iAlt;

  arma_cube gradient(nLons, nLats, nAlts);
  gradient.zeros();

  arma_cube one_minus_r2 = 1.0 - grid.dalt_ratio_sq_scgc;

  // Central part
  for (iAlt = 1; iAlt < nAlts - 1; iAlt++)
    gradient.slice(iAlt) =
      (value.slice(iAlt + 1)
       - one_minus_r2.slice(iAlt) % value.slice(iAlt)
       - grid.dalt_ratio_sq_scgc.slice(iAlt) % value.slice(iAlt - 1)) /
      (grid.dalt_lower_scgc.slice(iAlt + 1) %
       (1.0 + grid.dalt_ratio_scgc.slice(iAlt)));

  // lower boundary
  iAlt = 0;
  gradient.slice(iAlt) =
    (value.slice(iAlt + 1) - value.slice(iAlt)) /
    grid.dalt_lower_scgc.slice(iAlt);

  // upper boundary
  iAlt = nAlts - 1;
  gradient.slice(iAlt) =
    (value.slice(iAlt) - value.slice(iAlt - 1)) /
    grid.dalt_lower_scgc.slice(iAlt);
  return gradient;
}

// --------------------------------------------------------------------------
// Calculate the 4th order gradient in the altitudinal direction
// --------------------------------------------------------------------------

arma_cube calc_gradient_alt_4th(arma_cube value, Grid grid) {

  int64_t nLons = grid.get_nLons();
  int64_t nLats = grid.get_nLats();
  int64_t nAlts = grid.get_nAlts();
  int64_t iAlt;

  arma_cube gradient(nLons, nLats, nAlts);
  gradient.zeros();

  for (iAlt = 2; iAlt < nAlts - 2; iAlt++) {
    gradient.slice(iAlt) =
      grid.MeshCoefm2.slice(iAlt) * value.slice(iAlt - 2) +
      grid.MeshCoefm1.slice(iAlt) * value.slice(iAlt - 1) +
      grid.MeshCoefp0.slice(iAlt) * value.slice(iAlt) +
      grid.MeshCoefp1.slice(iAlt) * value.slice(iAlt + 1) +
      grid.MeshCoefp2.slice(iAlt) * value.slice(iAlt + 2);
  }

  return gradient;
}

// --------------------------------------------------------------------------
// Calculate the gradient in cubesphere spatial discretization
// --------------------------------------------------------------------------
std::vector<arma_cube> calc_gradient_cubesphere(arma_cube value, Grid grid) {
  // Must be used for cubesphere (Probably need a boolean check)
  int64_t nXs = grid.get_nY();
  int64_t nYs = grid.get_nX();
  int64_t nGCs = grid.get_nGCs();
  int64_t nAlts = grid.get_nAlts();

  // Initialize two arma cubes for return
  arma_cube grad_lon(nXs, nYs, nAlts);
  arma_cube grad_lat(nXs, nYs, nAlts);

  for (int64_t iAlt = 0; iAlt < nAlts; iAlt++) {
    /** Extract Grid Features **/
    // Addition: Get a copy of dx dy
    arma_mat curr_refx = grid.refx_scgc.slice(iAlt);
    arma_mat curr_refy = grid.refy_scgc.slice(iAlt);

    // Get some dx dy metrics from the grid
    precision_t dx = grid.drefx(iAlt);
    precision_t dy = grid.drefy(iAlt);

    // Get values of current level
    arma_mat curr_value = value.slice(iAlt);

    /** Calculate Gradient with Central Difference Scheme **/
    // Since Reference grid is orthogonal, we only need gradient along reference xy direction
    // Then we convert gradient from xy direction ot lat-lon direction

    arma_mat grad_x_curr(nXs, nYs);
    arma_mat grad_y_curr(nXs, nYs);

    // Calc gradient on x and y direction (since reference grid)
    // Only update interior cells
    // May vectorize for future improvements

    if (nGCs >=
        2) { // if more than 1 nGCs, we do fourth order, some foolproofing in case we go into debug hell
      for (int j = nGCs; j < nYs - nGCs; j++) {
        for (int i = nGCs; i < nXs - nGCs; i++) {
          grad_x_curr(i, j) = (-curr_value(i + 2, j) + 8 * curr_value(i + 1,
                                                                      j) - 8 * curr_value(i - 1, j) + curr_value(i - 2, j)) * (1. / 12. / dx);
          grad_y_curr(i, j) = (-curr_value(i, j + 2) + 8 * curr_value(i,
                                                                      j + 1) - 8 * curr_value(i, j - 1) + curr_value(i, j - 2)) * (1. / 12. / dy);
        }
      }
    } else { // otherwise we do second order
      for (int j = nGCs; j < nYs - nGCs; j++) {
        for (int i = nGCs; i < nXs - nGCs; i++) {
          grad_x_curr(i, j) = (curr_value(i + 1, j) - curr_value(i - 1,
                                                                 j)) * (1. / 2. / dx);
          grad_y_curr(i, j) = (curr_value(i, j + 1) - curr_value(i,
                                                                 j - 1)) * (1. / 2. / dy);
        }
      }
    }

    // We then use A transformation matrices to convert grad_xy to grad_latlon
    // Ref -> Physical, we use A matrix
    grad_lon.slice(iAlt) = grad_x_curr % grid.A11_inv_scgc.slice(
                             iAlt) + grad_y_curr % grid.A21_inv_scgc.slice(iAlt);
    grad_lat.slice(iAlt) = grad_x_curr % grid.A12_inv_scgc.slice(
                             iAlt) + grad_y_curr % grid.A22_inv_scgc.slice(iAlt);
  }

  // Not initializing with array like procedure in case I get bugs
  std::vector<arma_cube> gradient;
  gradient.push_back(grad_lon);
  gradient.push_back(grad_lat);

  return gradient;
}