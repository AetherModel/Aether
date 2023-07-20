// Copyright 2023, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

// Initial version: F. Cheng, July 2023

#include "../include/aether.h"

// Structure for transformation matrix A
struct A_mat{
    arma_mat& A11; 
    arma_mat& A12; 
    arma_mat& A21; 
    arma_mat& A22;
};

// Structure for transformation matrix A_inv
struct A_inv_mat{
    arma_mat& A11_inv; 
    arma_mat& A12_inv; 
    arma_mat& A21_inv; 
    arma_mat& A22_inv;
};

// --------------------------------------------------------------------------
// Convert spherical vector (velocities) to reference (contravariant) vector
// Units of the velocities and transformation laws must be the same
// u and v are spherical velocities
// u1 and u2 are contravariant velocities
// --------------------------------------------------------------------------
void sphvect2ref(arma_mat& u, arma_mat& v, arma_mat& u1, arma_mat& u2, A_inv_mat& A_inv_mat) {
    u1 = u % A_inv_mat.A11_inv + v % A_inv_mat.A12_inv;
    u2 = u % A_inv_mat.A21_inv + v % A_inv_mat.A22_inv;
}

// --------------------------------------------------------------------------
// Convert spherical vector (velocities) to reference (contravariant) vector
// Units of the velocities and transformation laws must be the same
// u and v are spherical velocities
// u1 and u2 are contravariant velocities
// --------------------------------------------------------------------------
void refvect2sph(arma_mat& u1, arma_mat& u2, arma_mat& u, arma_mat& v, A_mat& A_mat) {
    u = u1 % A_mat.A11 + u2 % A_mat.A12;
    v = u1 % A_mat.A21 + u2 % A_mat.A22;
}