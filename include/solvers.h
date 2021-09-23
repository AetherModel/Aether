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
double interpolate_1d_w_index(fvec values,
			      double interpolation_index,
			      int interpolation_type);
fmat interpolate_1d_w_index(std::vector<fmat> values,
			    double interpolation_index,
			    int interpolation_type);

fcube calc_gradient_lon(fcube value, Grid grid);
fcube calc_gradient_lat(fcube value, Grid grid);
fcube calc_gradient_alt(fcube value, Grid grid);
std::vector<fcube> calc_gradient_vector(fcube value_scgc, Grid grid);
  
#endif  // INCLUDE_SOLVERS_H_
