#ifndef INCLUDE_CUBESPHERE_H_
#define INCLUDE_CUBESPHERE_H_

#include "aether.h"
#include <armadillo>

/*************************************************
 * \brief A class to handle quad sphere grid logic.
 *************************************************/

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
precision_t cube2sphere_float(const precision_t x, const precision_t y, const precision_t z);

#endif  // INCLUDE_CUBESPHERE_H_
