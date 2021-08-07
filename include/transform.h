// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_TRANSFORM_H_
#define INCLUDE_TRANSFORM_H_

// The armadillo library is to allow the use of 3d cubes and other
// array types, with array math built in. This eliminates loops!
#include <armadillo>
using namespace arma;

void copy_cube_to_array(fcube cube_in,
                        float *array_out);

void copy_vector_to_array(std::vector<float> vector_in,
			  int64_t nElements,
			  float *array_out);

std::vector<fcube> transform_llr_to_xyz_3d(std::vector<fcube> llr);
std::vector<fcube> rotate_around_x_3d(std::vector<fcube> XYZ_in, float angle);
std::vector<fcube> rotate_around_y_3d(std::vector<fcube> XYZ_in, float angle);
std::vector<fcube> rotate_around_z_3d(std::vector<fcube> XYZ_in, float angle);

void transform_llr_to_xyz(float llr_in[3], float xyz_out[3]);
void transform_rot_z(float xyz_in[3], float angle_in, float xyz_out[3]);
void transform_rot_y(float xyz_in[3], float angle_in, float xyz_out[3]);
void transform_float_vector_to_array(std::vector<float> input,
                                     float output[3]);

void transform_vector_xyz_to_env(float xyz_in[3],
                                 float lon,
                                 float lat,
                                 float env_out[3]);

// These are not really transformations, but are placeholders

void vector_diff(float vect_in_1[3],
                 float vect_in_2[3],
                 float vect_out[3]);

void vector_add(float vect_in_1[3],
                 float vect_in_2[3],
                 float vect_out[3]);


#endif  // INCLUDE_TRANSFORM_H_


