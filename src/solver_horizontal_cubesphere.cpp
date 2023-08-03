// Copyright 2023, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

// Initial version: F. Cheng, July 2023

#include "../include/aether.h"

// --------------------------------------------------------------------------
// Convert spherical vector (velocities) to reference (contravariant) vector
// Units of the velocities and transformation laws must be the same
// u and v are spherical velocities
// u1 and u2 are contravariant velocities
// --------------------------------------------------------------------------
void sphvect2ref(arma_mat& u, arma_mat& v, arma_mat& u1, arma_mat& u2, mat_2x2 &A_inv_mat) {
    u1 = u % A_inv_mat.A11 + v % A_inv_mat.A12;
    u2 = u % A_inv_mat.A21 + v % A_inv_mat.A22;
}

// --------------------------------------------------------------------------
// Convert spherical vector (velocities) to reference (contravariant) vector
// Units of the velocities and transformation laws must be the same
// u and v are spherical velocities
// u1 and u2 are contravariant velocities
// --------------------------------------------------------------------------
void refvect2sph(arma_mat &u1, arma_mat &u2, arma_mat &u, arma_mat &v, mat_2x2 &A_mat) {
    u = u1 % A_mat.A11 + u2 % A_mat.A12;
    v = u1 % A_mat.A21 + u2 % A_mat.A22;
}