// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_SOLVERS_H_
#define INCLUDE_SOLVERS_H_

// The armadillo library is to allow the use of 3d cubes and other
// array types, with array math built in. This eliminates loops!
#include <armadillo>
using namespace arma;

arma_vec solver_conduction(arma_vec value,
			   arma_vec lambda,
			   arma_vec front,
			   precision_t dt,
			   arma_vec dx);


arma_cube solver_chemistry(arma_cube density,
			   arma_cube source,
			   arma_cube loss,
			   precision_t dt);

/// Set flag values that indicate whether the previous, next, closest,
/// or an interpolated value should be used.
const int iPrevious_ = 1;
const int iNext_ = 2;
const int iClosest_ = 3; 
const int iInterp_ = 4;

double interpolate_1d(double outX,
		      std::vector<double> inXs,
		      std::vector<double> inValues);

double interpolate_1d_get_index_doubles(double intime,
					std::vector<double> times);

// Overloading the interpolation function:
double interpolate_1d_w_index(std::vector<double> values,
			      double interpolation_index,
			      int interpolation_type);
double interpolate_1d_w_index(std::vector<float> values,
			      double interpolation_index,
			      int interpolation_type);
double interpolate_1d_w_index(std::vector<float> values,
			      float interpolation_index,
			      int interpolation_type);
double interpolate_1d_w_index(arma_vec values,
			      double interpolation_index,
			      int interpolation_type);
fmat interpolate_1d_w_index(std::vector<fmat> values,
			    double interpolation_index,
			    int interpolation_type);

arma_cube calc_gradient_lon(arma_cube value, Grid grid);
arma_cube calc_gradient_lat(arma_cube value, Grid grid);
arma_cube calc_gradient_alt(arma_cube value, Grid grid);
std::vector<arma_cube> calc_gradient_vector(arma_cube value_scgc, Grid grid);

// interpolation in 1D
precision_t linear_interpolation(const precision_t y0,
                                 const precision_t y1,
                                 const precision_t ratio);

// interpolation in 3D, data should be a cube of size 2-2-2
precision_t interpolate_unit_cube(const arma_cube &data,
                                  const precision_t xRatio,
                                  const precision_t yRatio,
                                  const precision_t zRatio);

#endif  // INCLUDE_SOLVERS_H_
