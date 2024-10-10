// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>

#include "aether.h"

// -----------------------------------------------------------------------
// This code solves the laplace equation in 1D:
//
// dV/dt = 1/front d/dx lambda dV/dx + source
//
// Some assumptions:
//  - Assume that lambda and front are scaled by radius squared!
//  - The spacing can be non-uniform.
//  - At this point, the bottom BC is fixed, while the top BC is zero gradient
//  - The dx variable is assumed to be distance between the CURRENT cell center
//    (i) and the cell center of the cell BELOW the current one (i-1).
// -----------------------------------------------------------------------

arma_vec solver_conduction(arma_vec value,
                           arma_vec lambda,
                           arma_vec front,
                           arma_vec source,
                           arma_vec dx,
                           precision_t dt,
                           int64_t nGCs,
                           bool return_diff) {

  int64_t nPts = value.n_elem;

  arma_vec di = lambda;
  arma_vec m = dt / front;

  // These are to allow for a stretched grid:
  // du is cell spacing in upper direction (so, lower, shifted by one):
  arma_vec du(nPts);
  du(span(0, nPts - 2)) = dx(span(1, nPts - 1));
  du(nPts - 1) = du(nPts - 2);
  // dl is lower cell spacing:
  arma_vec dl = dx;
  arma_vec r = du / dl;
  arma_vec du12 = du % du % (1 + r) % (1 + r);
  arma_vec du22 = 0.5 * (dl % du + du % du);
  arma_vec lou = di / du22;

  arma_vec conduction(nPts);
  conduction.zeros();

  int64_t i;

  for (i = nGCs; i < nPts - nGCs; i++)
    dl(i) = di(i + 1) - di(i - 1) * r(i) * r(i) - di(i) * (1.0 - r(i) * r(i));

  arma_vec a = di / du22 % r - dl / du12 % r % r;
  arma_vec c = di / du22 + dl / du12;
  arma_vec b = -1.0 / m - di / du22 % (1.0 + r) - dl / du12 % (1.0 - r % r);
  arma_vec d = -1.0 * (value / m + source % front * dt);

  // Lower BCs (fixed value):
  a(nGCs-1) = 0.0;
  b(nGCs-1) = -1.0;
  c(nGCs-1) = 0.0;
  d(nGCs-1) = -1.0 * value(nGCs-1);

  // Upper BCs:
  // This assumes a constant-gradient BC (need to change for ion and ele temps.

  i = nPts - nGCs;
  a(i) = 1.0 * (r(i) * (1.0 + r(i)) * di(i) * m(i) / du22(i));
  b(i) = -1.0 * (1.0 + r(i) * (1 + r(i)) * di(i) * m(i) / du22(i));
  c(i) = 0.0;
  d(i) = -1.0 * value(i);

  arma_vec cp(nPts, fill::zeros);
  arma_vec dp(nPts, fill::zeros);
  arma_vec result(nPts, fill::zeros);

  cp(nGCs - 1) = c(nGCs - 1) / b(nGCs - 1);

  for (i = nGCs; i <= nPts - nGCs; i++)
    cp(i) = c(i) / (b(i) - cp(i - 1) * a(i));

  dp(nGCs - 1) = d(nGCs - 1) / b(nGCs - 1);

  for (i = nGCs; i <= nPts - nGCs; i++)
    dp(i) = (d(i) - dp(i - 1) * a(i)) / (b(i) - cp(i - 1) * a(i));

  result(nPts - nGCs) = dp(nPts - nGCs);

  for (i = nPts - nGCs - 1; i > nGCs - 1; i--)
    result(i) = dp(i) - cp(i) * result(i + 1);

  if (return_diff) {
    conduction = result - value;
    for (i = 0; i < nGCs; i++) {
      conduction(i) = 0.0;
      conduction(nPts - i - 1) = 0.0;
    }
  } else {
    conduction = result;
    for (i = 0; i < nGCs; i++) {
      conduction(i) = value(i);
      conduction(nPts - nGCs + i) = conduction(nPts - nGCs - 1);
    }
  }

  return conduction;
}
