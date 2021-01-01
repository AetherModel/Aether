// (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef AETHER_INCLUDE_TRANSFORM_H_
#define AETHER_INCLUDE_TRANSFORM_H_

#include <armadillo>
using namespace arma;

void copy_cube_to_array(fcube cube_in,
			float *array_out);

void transform_llr_to_xyz_3d(fcube lat3d,
			     fcube lon3d,
			     fcube r3d,
			     fcube &x3d,
			     fcube &y3d,
			     fcube &z3d);

void transform_llr_to_xyz(float llr_in[3], float xyz_out[3]);
void transform_rot_z(float xyz_in[3], float angle_in, float xyz_out[3]);
void transform_rot_y(float xyz_in[3], float angle_in, float xyz_out[3]);
void transform_float_vector_to_array(std::vector<float> input,
				     float output[3]);

void transform_vector_xyz_to_env(float xyz_in[3],
				 float lon,
				 float lat,
				 float env_out[3]);

// These are not really transformations, but are 

void vector_diff(float vect_in_1[3],
		 float vect_in_2[3],
		 float vect_out[3]);


#endif // AETHER_INCLUDE_TRANSFORM_H_


