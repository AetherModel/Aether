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
