// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md


#ifndef INCLUDE_INIT_MAG_GRID_H_
#define INCLUDE_INIT_MAG_GRID_H_

#include "aether.h"
#include "planets.h"
#include "grid.h"

bool init_dipole_grid(Grid &mGrid, Planets planet);

// Analytic solution to get from q,p dipole coords to r,theta
// q coordinate along b-field line
// p l-shell 
// return (r,theta)
std::pair<precision_t, precision_t> qp_to_r_theta(precision_t q, precision_t p);

// Take limits & specs of field line, fill it with points.
// std:: tuple <arma_mat, arma_mat, arma_vec> Grid::fill_field_lines (arma_vec lShells, 
//                                                 int64_t nAlts,
//                                                 int64_t nLats,
//                                                 precision_t Gamma);


// convert mag to geographic
std::vector <arma_cube> mag_to_geo(arma_cube magLon, 
                                   arma_cube magLat,
                                   arma_cube magAlt,
                                   Planets planet);


#endif  // INCLUDE_INIT_GEO_GRID_H_

