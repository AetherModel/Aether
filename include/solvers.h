// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_SOLVERS_H_
#define INCLUDE_SOLVERS_H_

// The armadillo library is to allow the use of 3d cubes and other
// array types, with array math built in. This eliminates loops!
#include <armadillo>
using namespace arma;

struct projection_struct {

  arma_mat gradLR;
  arma_mat gradDU;
  arma_mat R;
  arma_mat L;
  arma_mat U;
  arma_mat D;
  arma_mat grad_edge_LR;
  arma_mat grad_edge_DU;
};

arma_vec limiter_mc(arma_vec &left, arma_vec &right, int64_t nPts, int64_t nGCs);
arma_vec calc_grad_1d(arma_vec &values,
		      						arma_vec &x,
		      						int64_t nPts,
		      						int64_t nGCs);
arma_mat calc_grad(arma_mat values, arma_mat x, int64_t nGCs, bool DoX);

void advect(Grid &grid,
            Times &time,
            Neutrals &neutrals);

arma_vec solver_conduction(
			   arma_vec value,
			   arma_vec lambda,
			   arma_vec front,
			   arma_vec source,
			   arma_vec dx,
			   precision_t dt,
			   int64_t nGCs,
			   bool return_diff);

arma_cube solver_chemistry(arma_cube density,
			   arma_cube source,
			   arma_cube loss,
			   precision_t dt);

std::vector<arma_cube> coriolis(std::vector<arma_cube> velocity,
                                precision_t rotation_rate,
                                arma_cube lat_scgc);

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
std::vector<arma_cube> calc_gradient_cubesphere(arma_cube value, Grid grid);
arma_cube calc_gradient_alt_4th(arma_cube value, Grid grid);

// interpolation in 1D
precision_t linear_interpolation(const precision_t y0,
                                 const precision_t y1,
                                 const precision_t ratio);

// interpolation in 3D, data should be a cube of size 2-2-2
precision_t interpolate_unit_cube(const arma_cube &data,
                                  const precision_t xRatio,
                                  const precision_t yRatio,
                                  const precision_t zRatio);

precision_t limiter_mc(precision_t dUp,
		       precision_t dDown,
		       precision_t beta);


  /**********************************************************************
     \brief Calculate dt (cell size / cMax) in each direction, and take min
     \param dt returns the neutral time-step
     \param grid The grid to define the neutrals on
   **/
  precision_t calc_dt(Grid grid, std::vector<arma_cube> cMax_vcgc);
  precision_t calc_dt_sphere(Grid grid, std::vector<arma_cube> cMax_vcgc);  
  precision_t calc_dt_cubesphere(Grid grid, std::vector<arma_cube> cMax_vcgc);  
  precision_t calc_dt_vertical(Grid grid, std::vector<arma_cube> cMax_vcgc);  

#endif  // INCLUDE_SOLVERS_H_
