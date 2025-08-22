// Copyright 2024, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

// Initial version: F. Cheng, Feb 2024

#include "aether.h"

arma_vec Cubesphere_tools::limiter_mc(arma_vec &left,
                                      arma_vec &right,
                                      int64_t nPts,
                                      int64_t nGCs) {

  precision_t beta = 0.8;

  arma_vec s = left % right;
  arma_vec combined = (left + right) * 0.5;

  left = left * beta;
  right = right * beta;
  arma_vec limited = left;

  for (int64_t i = 1; i < nPts + 2 * nGCs - 1; i++) {
    if (s(i) < 0) {
      // Sign < 0 means opposite signed left and right:
      limited(i) = 0.0;
    } else {
      if (left(i) > 0 && right(i) > 0) {
        if (right(i) < limited(i))
          limited(i) = right(i);

        if (combined(i) < limited(i))
          limited(i) = combined(i);
      } else {
        if (right(i) > limited(i))
          limited(i) = right(i);

        if (combined(i) > limited(i))
          limited(i) = combined(i);
      }
    }
  }

  return limited;
}

void Cubesphere_tools::print(arma_vec values) {
  int64_t nP = values.n_elem;

  for (int64_t i = 0; i < nP; i++)
    std::cout << values(i) << " ";

  std::cout << "\n";
}

// ---------------------------------------------------------
// calc gradients at centers
//   - values and x defined at centers
// ---------------------------------------------------------

arma_vec Cubesphere_tools::calc_grad_1d(arma_vec &values,
                                        arma_vec &x,
                                        int64_t nPts,
                                        int64_t nGCs) {

  arma_vec gradients = values * 0.0;
  arma_vec gradL = values * 0.0;
  arma_vec gradR = values * 0.0;

  precision_t factor1 = 0.625;
  precision_t factor2 = 0.0416667;
  precision_t h;

  int64_t i;
  arma_vec hv = values * 0.0;

  i = nGCs - 1;
  h = 2.0 / (x(i + 1) - x(i));
  gradR(i) = h * (factor1 * (values(i + 1) - values(i)) -
                  factor2 * (values(i + 2) - values(i - 1)));
  gradL(i) = (values(i) - values(i - 1)) / (x(i) - x(i - 1));

  // This is attempting to vectorize the problem, but it seems to be slower?
  //  int64_t iS = nGCs;
  //  int64_t iE = nPts + nGCs - 1;
  //  hv.rows(iS, iE) = 2.0 / (x.rows(iS, iE) - x.rows(iS-1, iE-1));
  //  gradL.rows(iS, iE) = hv.rows(iS,iE) % (factor1 * (values.rows(iS, iE) -
  //                values.rows(iS-1, iE-1)) -
  //           factor2 * (values.rows(iS+1, iE+1) -
  //                values.rows(iS-2, iE-2)));
  //  hv.rows(iS, iE) = 2.0 / (x.rows(iS+1, iE+1) - x.rows(iS, iE));
  //  gradR.rows(iS, iE) = hv.rows(iS,iE) % (factor1 * (values.rows(iS+1, iE+1) -
  //                values.rows(iS, iE)) -
  //           factor2 * (values.rows(iS+2, iE+2) -
  //                values.rows(iS-1, iE-1)));

  for (i = nGCs; i < nPts + nGCs; i++) {
    h = 2.0 / (x(i) - x(i - 1));
    gradL(i) = h * (factor1 * (values(i) - values(i - 1)) -
                    factor2 * (values(i + 1) - values(i - 2)));
    h = 2.0 / (x(i + 1) - x(i));
    gradR(i) = h * (factor1 * (values(i + 1) - values(i)) -
                    factor2 * (values(i + 2) - values(i - 1)));
  }

  i = nPts + nGCs;
  h = 2.0 / (x(i) - x(i - 1));
  gradL(i) = h * (factor1 * (values(i) - values(i - 1)) -
                  factor2 * (values(i + 1) - values(i - 2)));
  gradR(i) = (values(i + 1) - values(i)) / (x(i + 1) - x(i));

  gradients = Cubesphere_tools::limiter_mc(gradL, gradR, nPts, nGCs);

  return gradients;
}

// ---------------------------------------------------------
// calc gradients at centers for 2d matrices
//   - values and x defined at centers
// ---------------------------------------------------------

arma_mat Cubesphere_tools::calc_grad(arma_mat values,
                                     arma_mat x,
                                     int64_t nGCs,
                                     bool DoX) {

  arma_mat v2d, x2d;

  if (DoX) {
    v2d = values;
    x2d = x;
  } else {
    v2d = values.t();
    x2d = x.t();
  }

  int64_t nX = v2d.n_rows;
  int64_t nY = v2d.n_cols;
  arma_mat grad2d = v2d * 0.0;

  int64_t nPts = nX - 2 * nGCs;
  arma_vec values1d(nX);
  arma_vec x1d(nX);

  for (int64_t j = 1; j < nY - 1; j++) {
    values1d = v2d.col(j);
    x1d = x2d.col(j);
    grad2d.col(j) = calc_grad_1d(values1d, x1d, nPts, nGCs);
  }

  arma_mat gradients;

  if (DoX)
    gradients = grad2d;
  else
    gradients = grad2d.t();

  return gradients;
}

// ---------------------------------------------------------
// Project gradients + values to the right face, from the left
//   returned values are on the i - 1/2 edges
//     (between i-1 and i cell center)
// ---------------------------------------------------------

arma_mat Cubesphere_tools::project_from_left(arma_mat values,
                                             arma_mat gradients,
                                             arma_mat x_centers,
                                             arma_mat x_edges,
                                             int64_t nGCs) {

  int64_t nX = values.n_rows;
  int64_t nY = values.n_cols;

  // Define at edges:
  arma_mat projected(nX + 1, nY);
  projected.zeros();

  // no gradient in the 0 or iEnd cells
  for (int64_t j = 0; j < nY; j++) {
    for (int64_t i = 1; i < nX - 1; i++) {
      projected(i + 1, j) = values(i, j) +
                            gradients(i, j) * (x_edges(i + 1, j) - x_centers(i, j));
    }

    projected(1, j) = projected(2, j);
    projected(0, j) = projected(1, j);
    projected(nX, j) = projected(nX - 1, j);
  }

  return projected;
}

// ---------------------------------------------------------
// Project gradients + values to the left face, from the right
//   returned values are on the i - 1 edges
//     (between i-1 and i cell center)
// ---------------------------------------------------------

arma_mat Cubesphere_tools::project_from_right(arma_mat values,
                                              arma_mat gradients,
                                              arma_mat x_centers,
                                              arma_mat x_edges,
                                              int64_t nGCs) {
  int64_t nX = values.n_rows;
  int64_t nY = values.n_cols;

  // Define at edges:
  arma_mat projected(nX + 1, nY);
  projected.zeros();

  // no gradient in the 0 or iEnd cells
  for (int64_t j = 0; j < nY; j++) {
    for (int64_t i = 1; i < nX - 1; i++) {
      projected(i, j) = values(i, j) +
                        gradients(i, j) * (x_edges(i, j) - x_centers(i, j));
    }

    projected(0, j) = projected(1, j);
    projected(nX - 1, j) = projected(nX - 2, j);
    projected(nX, j) = projected(nX - 1, j);
  }

  return projected;
}

// ---------------------------------------------------------
// Limiter on values
//   projected is assumed to be on the edge between the
//   i-1 and i cell (i-1/2)
//   limited is returned at edges
// ---------------------------------------------------------

arma_vec Cubesphere_tools::limiter_value(arma_vec projected,
                                         arma_vec values,
                                         int64_t nPts,
                                         int64_t nGCs) {

  int64_t iStart = 0;
  int64_t iEnd = nPts + 2 * nGCs;

  arma_vec limited = projected;

  precision_t mini, maxi;

  for (int64_t i = iStart + 1; i < iEnd - 1; i++) {

    mini = values(i - 1);

    if (values(i) < mini)
      mini = values(i);

    maxi = values(i - 1);

    if (values(i) > maxi)
      maxi = values(i);

    if (limited(i) < mini)
      limited(i) = mini;

    if (limited(i) > maxi)
      limited(i) = maxi;
  }

  return limited;
}

// // ---------------------------------------------------------
// // take gradients and project to all edges
// // ---------------------------------------------------------

// projection_struct Cubesphere_tools::project_to_edges(arma_mat &values,
//                                                      arma_mat &x_centers, arma_mat &x_edges,
//                                                      arma_mat &y_centers, arma_mat &y_edges,
//                                                      int64_t nGCs) {

//   int64_t nX = values.n_rows;
//   int64_t nY = values.n_cols;

//   projection_struct proj;

//   proj.gradLR = calc_grad(values, x_centers, nGCs, true);
//   proj.gradDU = calc_grad(values.t(), y_centers.t(), nGCs, true).t();

//   proj.R = project_from_left(values, proj.gradLR,
//                              x_centers, x_edges, nGCs);
//   // Left side of edge from left
//   proj.L = project_from_right(values, proj.gradLR,
//                               x_centers, x_edges, nGCs);
//   // Up side of edge from down (left)
//   proj.U = project_from_left(values.t(), proj.gradDU.t(),
//                              y_centers.t(), y_edges.t(), nGCs)
//            .t();
//   // Down side of edge from up (right)
//   proj.D = project_from_right(values.t(), proj.gradDU.t(),
//                               y_centers.t(), y_edges.t(), nGCs)
//            .t();

//   return proj;
// }
