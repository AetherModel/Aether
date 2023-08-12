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
      (grid.dalt_lower_scgc.slice(iAlt + 1) % (1.0 + grid.dalt_ratio_scgc.slice(
                                                 iAlt)));

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


