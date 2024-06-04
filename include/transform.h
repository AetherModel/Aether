// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_TRANSFORM_H_
#define INCLUDE_TRANSFORM_H_

// The armadillo library is to allow the use of 3d cubes and other
// array types, with array math built in. This eliminates loops!
#include <armadillo>
#include "aether.h"
using namespace arma;

std::string mklower(std::string inString);
std::string mkupper(std::string inString);

void copy_cube_to_array(arma_cube cube_in,
                        float *array_out);
void copy_mat_to_array(arma_mat mat_in,
                        float *array_out,
                        bool isFortran);
void copy_array_to_mat(float *array_in,
                       arma_mat &mat_out,
                       bool isFortran);

void copy_vector_to_array(std::vector<float> vector_in,
			  int64_t nElements,
			  float *array_out);

// This is needed when sending strings to Fortran.
// We do this by copying the ascii numbers into an integer array,
// then in fortran copy them back into a character array.
int* copy_string_to_int(std::string inString);

arma_cube calc_magnitude(std::vector<arma_cube> xyz);
std::vector<arma_cube> transform_llr_to_xyz_3d(std::vector<arma_cube> llr);
std::vector<arma_cube> transform_xyz_to_llr_3d(std::vector<arma_cube> xyz);
std::vector<arma_cube> rotate_around_x_3d(std::vector<arma_cube> XYZ_in, precision_t angle);
std::vector<arma_cube> rotate_around_y_3d(std::vector<arma_cube> XYZ_in, precision_t angle);
std::vector<arma_cube> rotate_around_z_3d(std::vector<arma_cube> XYZ_in, precision_t angle);

void transform_llr_to_xyz(precision_t llr_in[3], precision_t xyz_out[3]);
void transform_rot_z(precision_t xyz_in[3], precision_t angle_in, precision_t xyz_out[3]);
void transform_rot_y(precision_t xyz_in[3], precision_t angle_in, precision_t xyz_out[3]);
void transform_float_vector_to_array(std::vector<float> input,
                                     precision_t output[3]);

void transform_vector_xyz_to_env(precision_t xyz_in[3],
                                 precision_t lon,
                                 precision_t lat,
                                 precision_t env_out[3]);

// These are not really transformations, but are placeholders

void vector_diff(precision_t vect_in_1[3],
                 precision_t vect_in_2[3],
                 precision_t vect_out[3]);

// -----------------------------------------------------------------------
// Simple 3-element vector addition
// -----------------------------------------------------------------------

void vector_add(precision_t vect_in_1[3],
                precision_t vect_in_2[3],
                precision_t vect_out[3]);

#endif  // INCLUDE_TRANSFORM_H_


