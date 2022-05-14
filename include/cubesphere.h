// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_CUBESPHERE_H_
#define INCLUDE_CUBESPHERE_H_

#include "aether.h"
#include <memory>
#include <armadillo>

/*************************************************
 * \brief A namespace with all quad sphere grid logic.
 *************************************************/
namespace CubeSphere {

/** Returns x, y and z locations of a cube from -1 to +1.
 **/
arma_mat map_cube2sphere(const arma_mat& cube);

/** Helper function the transforms cartesian coordinates to cube sphere grid
 * vertices in cartesian coordinates.
 */
arma_vec cube2sphere(const arma_vec& x, const arma_vec& y, const arma_vec& z);

/** Helper function the transforms cartesian coordinates to cube sphere grid
 * vertices in cartesian coordinates.
 */
precision_t cube2sphere(const precision_t& x, const precision_t& y, const precision_t& z);

/** Helper function that transforms cartesian coordinates of cube to cube sphere
 * using an x, y, z vector as input.
 */
arma_vec cube2sphere(const arma_vec& p);

/** Returns normalized x, y, z locations based on number of subdivisions per
 * face
 */
arma_cube sphere_face_locations(const uint64_t& face, const uint64_t& subdivisions);

/** Returns the normalized x, y, z coordinate locations of a cell in the
 * cubesphere grid.
 */
arma_vec cartesian_location(const uint64_t& face, const uint64_t& subdivisions,
                            const uint64_t& row, const uint64_t& column);

/// The normalized origins of each face of the cube (i.e. corner)
static const arma_mat ORIGINS = {
    {-1.0, -1.0, -1.0},
    { 1.0, -1.0, -1.0},
    { 1.0,  1.0, -1.0},
    {-1.0,  1.0, -1.0},
    {-1.0, -1.0, -1.0},
    { 1.0, -1.0,  1.0},
};

/// Normalized right steps in cube
static const arma_mat RIGHTS = {
		{ 2.0,  0.0, 0.0},
		{ 0.0,  2.0, 0.0},
		{-2.0,  0.0, 0.0},
		{ 0.0, -2.0, 0.0},
		{ 0.0,  2.0, 0.0},
		{ 0.0,  2.0, 0.0}
};

/// Normalized right steps in cube
static const arma_mat UPS = {
		{ 0.0, 0.0, 2.0},
		{ 0.0, 0.0, 2.0},
		{ 0.0, 0.0, 2.0},
		{ 0.0, 0.0, 2.0},
		{ 2.0, 0.0, 0.0},
		{-2.0, 0.0, 0.0}
};

} // CubeSphere::

#endif  // INCLUDE_CUBESPHERE_H_
