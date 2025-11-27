// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "../include/aether.h"

// --------------------------------------------------------------------------
// Calculate the gradient in all directions, returning a vector
// It is assumed that this is lon, lat, rad
// --------------------------------------------------------------------------

std::vector<arma_cube> calc_gradient_vector(arma_cube value_scgc, Grid &grid) {

  std::vector<arma_cube> gradient_vcgc;

  if (report.test_verbose(4)) {
    std::cout << "grid shape (1 = sphere; 2 = cube; 3 = dipole) : " <<
              grid.iGridShape_ << "\n";
    display_vector("gradient, value : ", value_scgc.tube(9, 9));
  }

  if (grid.iGridShape_ == iCubesphere_)
    gradient_vcgc = calc_gradient_cubesphere(value_scgc, grid);
  else if (grid.iGridShape_ == iDipole_)
    gradient_vcgc = calc_gradient_dipole(value_scgc, grid);
  else {

    report.print(4, "Going into calc_gradient_lon");
    gradient_vcgc.push_back(calc_gradient2o_i(value_scgc, grid));

    if (report.test_verbose(4))
      display_vector("gradient[0] : ", gradient_vcgc[0].tube(9, 9));

    report.print(4, "Going into calc_gradient_lat");
    gradient_vcgc.push_back(calc_gradient2o_j(value_scgc, grid));

    if (report.test_verbose(4))
      display_vector("gradient[1] : ", gradient_vcgc[1].tube(9, 9));

    report.print(4, "Going into calc_gradient_alt");
    gradient_vcgc.push_back(calc_gradient_alt(value_scgc, grid));

    if (report.test_verbose(4))
      display_vector("gradient[2] : ", gradient_vcgc[2].tube(9, 9));
  }

  return gradient_vcgc;
}

// --------------------------------------------------------------------------
// Routines related to native i - direction
//  - 2nd order uniform grid
//  - 4th order uniform grid
//  - 2nd order stretched grid
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
// Calculate the 2nd order gradient in the native i direction
//   - these formulas assume that the grid is uniform.
// --------------------------------------------------------------------------

arma_cube calc_gradient2o_i(arma_cube value, Grid &grid) {

  std::string function = "calc_gradient2o_i";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int64_t nX = grid.get_nX();
  int64_t nY = grid.get_nY();
  int64_t nZ = grid.get_nZ();
  int64_t iX;

  arma_cube gradient(nX, nY, nZ);
  gradient.zeros();

  if (grid.get_HasXdim()) {
    // Interior:
    for (iX = 1; iX < nX - 1; iX++)
      gradient.row(iX) =
        (value.row(iX + 1) - value.row(iX - 1)) /
        (2.0 * grid.di_center_m_scgc.row(iX));

    // Lower (one sided):
    iX = 0;
    gradient.row(iX) =
      (value.row(iX + 1) - value.row(iX)) /
      grid.di_center_m_scgc.row(iX);

    // Upper (one sided):
    iX = nX - 1;
    gradient.row(iX) =
      (value.row(iX) - value.row(iX - 1)) /
      grid.di_center_m_scgc.row(iX);
  }
  report.exit(function);
  return gradient;
}

// --------------------------------------------------------------------------
// Calculate the 4th gradient in the native i direction
//   - these formulas assume that the grid is uniform.
// --------------------------------------------------------------------------

arma_cube calc_gradient4o_i(arma_cube value, Grid &grid) {

  int64_t nX = grid.get_nX();
  int64_t nY = grid.get_nY();
  int64_t nZ = grid.get_nZ();
  int64_t iX;

  arma_cube gradient(nX, nY, nZ);
  gradient.zeros();

  if (grid.get_HasXdim()) {
    // Interior, 5 point sencil:
    for (iX = 2; iX < nX - 2; iX++)
      gradient.row(iX) = (-value.row(iX + 2) +
                          8 * value.row(iX + 1) -
                          8 * value.row(iX - 1) +
                          value.row(iX - 2)) /
                         (12. * grid.di_center_m_scgc.row(iX));

    // Points just inside edges (2nd order):
    iX = 1;
    gradient.row(iX) =
      (value.row(iX + 1) - value.row(iX - 1)) /
      (2 * grid.di_center_m_scgc.row(iX));
    iX = nX - 2;
    gradient.row(iX) =
      (value.row(iX + 1) - value.row(iX - 1)) /
      (2 * grid.di_center_m_scgc.row(iX));

    // Points at the edges (1st order):
    iX = 0;
    gradient.row(iX) =
      (value.row(iX + 1) - value.row(iX)) /
      grid.di_center_m_scgc.row(iX);
    iX = nX - 1;
    gradient.row(iX) =
      (value.row(iX) - value.row(iX - 1)) /
      grid.di_center_m_scgc.row(iX);
  }

  return gradient;
}

// --------------------------------------------------------------------------
// Calculate the 2nd order gradient in the native i-direction,
//   assuming a stretched grid
// --------------------------------------------------------------------------

arma_cube calc_gradient_stretched_i(arma_cube value, Grid grid) {

  int64_t nX = grid.get_nX();
  int64_t nY = grid.get_nY();
  int64_t nZ = grid.get_nZ();
  int64_t iX;

  arma_cube gradient(nX, nY, nZ);
  gradient.zeros();

  if (grid.get_HasXdim()) {
    // Central part
    for (iX = 1; iX < nX - 1; iX++)
      gradient.row(iX) =
        (value.row(iX + 1)
         - grid.di_one_minus_r2.row(iX) % value.row(iX)
         - grid.di_ratio_sq.row(iX) % value.row(iX - 1)) /
        (grid.di_edge_m.row(iX + 1) %
         (1.0 + grid.di_ratio.row(iX)));

    // Points at the edges (1st order):
    iX = 0;
    gradient.row(iX) =
      (value.row(iX + 1) - value.row(iX)) /
      grid.di_center_m_scgc.row(iX);
    iX = nX - 1;
    gradient.row(iX) =
      (value.row(iX) - value.row(iX - 1)) /
      grid.di_center_m_scgc.row(iX);
  }

  return gradient;
}



// --------------------------------------------------------------------------
// Calculate the gradient in the longitudinal direction
// --------------------------------------------------------------------------

arma_cube calc_gradient_lon(arma_cube value, Grid &grid) {
  return calc_gradient2o_i(value, grid);
}

// --------------------------------------------------------------------------
// Calculate the 2nd order gradient in the native j direction
//   - these formulas assume that the grid is uniform.
// --------------------------------------------------------------------------

arma_cube calc_gradient2o_j(arma_cube value, Grid &grid) {

  std::string function = "calc_gradient2o_j";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int64_t nX = grid.get_nX();
  int64_t nY = grid.get_nY();
  int64_t nZ = grid.get_nZ();
  int64_t iY;

  arma_cube gradient(nX, nY, nZ);
  gradient.zeros();

  if (grid.get_HasYdim()) {
    // Interior:
    for (iY = 1; iY < nY - 1; iY++)
      gradient.col(iY) =
        (value.col(iY + 1) - value.col(iY - 1)) /
        (2 * grid.dj_center_m_scgc.col(iY));

    // Lower (one sided):
    iY = 0;
    gradient.col(iY) =
      (value.col(iY + 1) - value.col(iY)) /
      grid.dj_center_m_scgc.col(iY);

    // Upper (one sided):
    iY = nY - 1;
    gradient.col(iY) =
      (value.col(iY) - value.col(iY - 1)) /
      grid.dj_center_m_scgc.col(iY);
  }
  report.exit(function);
  return gradient;
}

// --------------------------------------------------------------------------
// Calculate the 4th order gradient in the native j direction
//   - these formulas assume that the grid is uniform.
// --------------------------------------------------------------------------

arma_cube calc_gradient4o_j(arma_cube value, Grid &grid) {

  int64_t nX = grid.get_nX();
  int64_t nY = grid.get_nY();
  int64_t nZ = grid.get_nZ();
  int64_t iY;

  arma_cube gradient(nX, nY, nZ);
  gradient.zeros();

  if (grid.get_HasYdim()) {
    // Interior:
    for (iY = 2; iY < nY - 2; iY++)
      gradient.col(iY) = (-value.col(iY + 2) +
                          8 * value.col(iY + 1) -
                          8 * value.col(iY - 1) +
                          value.col(iY - 2)) /
                         (12. * grid.dj_center_m_scgc.col(iY));

    // Points just inside edges (2nd order):
    iY = 1;
    gradient.col(iY) =
      (value.col(iY + 1) - value.col(iY - 1)) /
      (2 * grid.dj_center_m_scgc.col(iY));
    iY = nY - 2;
    gradient.col(iY) =
      (value.col(iY + 1) - value.col(iY - 1)) /
      (2 * grid.dj_center_m_scgc.col(iY));

    // Lower (one sided):
    iY = 0;
    gradient.col(iY) =
      (value.col(iY + 1) - value.col(iY)) /
      grid.dj_center_m_scgc.col(iY);

    // Upper (one sided):
    iY = nY - 1;
    gradient.col(iY) =
      (value.col(iY) - value.col(iY - 1)) /
      grid.dj_center_m_scgc.col(iY);
  }

  return gradient;
}

// --------------------------------------------------------------------------
// Calculate the 2nd order gradient in the native j-direction,
//   assuming a stretched grid
// --------------------------------------------------------------------------

arma_cube calc_gradient_stretched_j(arma_cube value, Grid grid) {

  int64_t nX = grid.get_nX();
  int64_t nY = grid.get_nY();
  int64_t nZ = grid.get_nZ();
  int64_t iY;

  arma_cube gradient(nX, nY, nZ);
  gradient.zeros();

  if (grid.get_HasYdim()) {
    // Central part
    for (iY = 1; iY < nY - 1; iY++)
      gradient.col(iY) =
        (value.col(iY + 1)
         - grid.dj_one_minus_r2.col(iY) % value.col(iY)
         - grid.dj_ratio_sq.col(iY) % value.col(iY - 1)) /
        (grid.dj_edge_m.col(iY + 1) %
         (1.0 + grid.dj_ratio.col(iY)));

    // Points at the edges (1st order):
    iY = 0;
    gradient.col(iY) =
      (value.col(iY + 1) - value.col(iY)) /
      grid.dj_center_m_scgc.col(iY);
    iY = nY - 1;
    gradient.col(iY) =
      (value.col(iY) - value.col(iY - 1)) /
      grid.dj_center_m_scgc.col(iY);
  }

  return gradient;
}

// --------------------------------------------------------------------------
// Calculate the gradient in the latitudinal direction
// --------------------------------------------------------------------------

arma_cube calc_gradient_lat(arma_cube value, Grid &grid) {
  return calc_gradient2o_j(value, grid);
}

// --------------------------------------------------------------------------
// Calculate the 2nd order gradient in the native k direction
//   - these formulas assume that the grid is uniform.
// --------------------------------------------------------------------------

arma_cube calc_gradient2o_k(arma_cube value, Grid &grid) {

  int64_t nX = grid.get_nX();
  int64_t nY = grid.get_nY();
  int64_t nZ = grid.get_nZ();
  int64_t iZ;

  arma_cube gradient(nX, nY, nZ);
  gradient.zeros();

  if (grid.get_HasZdim()) {
    // Interior:
    for (iZ = 1; iZ < nZ - 1; iZ++)
      gradient.slice(iZ) =
        (value.slice(iZ + 1) - value.slice(iZ - 1)) /
        (2 * grid.dk_center_m_scgc.slice(iZ));

    // Lower (one sided):
    iZ = 0;
    gradient.slice(iZ) =
      (value.slice(iZ + 1) - value.slice(iZ)) /
      grid.dk_center_m_scgc.slice(iZ);

    // Upper (one sided):
    iZ = nZ - 1;
    gradient.slice(iZ) =
      (value.slice(iZ) - value.slice(iZ - 1)) /
      grid.dk_center_m_scgc.slice(iZ);
  }

  return gradient;
}

// --------------------------------------------------------------------------
// Calculate the gradient in the altitudinal direction
// --------------------------------------------------------------------------

arma_cube calc_gradient_alt(arma_cube value, Grid &grid) {

  std::string function = "calc_gradient_alt";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int64_t nX = grid.get_nX();
  int64_t nY = grid.get_nY();
  int64_t nZ = grid.get_nZ();
  int64_t iK;

  arma_cube gradient(nX, nY, nZ);
  gradient.zeros();

  if (grid.get_HasZdim()) {
    // Central part
    for (iK = 1; iK < nZ - 1; iK++)
      gradient.slice(iK) =
        (value.slice(iK + 1)
         - grid.dk_one_minus_r2.slice(iK) % value.slice(iK)
         - grid.dk_ratio_sq.slice(iK) % value.slice(iK - 1)) /
        (grid.dk_edge_m.slice(iK + 1) %
         (1.0 + grid.dk_ratio.slice(iK)));

    // lower boundary
    iK = 0;
    gradient.slice(iK) =
      (value.slice(iK + 1) - value.slice(iK)) /
      grid.dk_edge_m.slice(iK);

    // upper boundary
    iK = nZ - 1;
    gradient.slice(iK) =
      (value.slice(iK) - value.slice(iK - 1)) /
      grid.dk_edge_m.slice(iK);
  }
  report.exit(function);

  return gradient;
}

// --------------------------------------------------------------------------
// Calculate the 4th order gradient in the altitudinal direction
// --------------------------------------------------------------------------

arma_cube calc_gradient_alt_4th(arma_cube value, Grid &grid) {

  int64_t nLons = grid.get_nLons();
  int64_t nLats = grid.get_nLats();
  int64_t nAlts = grid.get_nAlts();
  int64_t iAlt;

  arma_cube gradient(nLons, nLats, nAlts);
  gradient.zeros();

  for (iAlt = 2; iAlt < nAlts - 2; iAlt++) {
    gradient.slice(iAlt) =
      grid.MeshCoefm2.slice(iAlt) * value.slice(iAlt - 2) +
      grid.MeshCoefm1.slice(iAlt) * value.slice(iAlt - 1) +
      grid.MeshCoefp0.slice(iAlt) * value.slice(iAlt) +
      grid.MeshCoefp1.slice(iAlt) * value.slice(iAlt + 1) +
      grid.MeshCoefp2.slice(iAlt) * value.slice(iAlt + 2);
  }

  return gradient;
}

// --------------------------------------------------------------------------
// Calculate the 3rd order (on-sided) gradient in the altitudinal direction
//   - this is only defined for the bottom ghostcells!
// --------------------------------------------------------------------------

arma_mat project_onesided_alt_3rd(arma_cube value, Grid &grid, int64_t iAlt) {

  int64_t nLons = grid.get_nLons();
  int64_t nLats = grid.get_nLats();

  arma_mat gradient(nLons, nLats), valueOut(nLons, nLats);
  /*
  gradient =
      grid.MeshCoef1s3rdp1.slice(iAlt) % value.slice(iAlt + 1) +
      grid.MeshCoef1s3rdp2.slice(iAlt) % value.slice(iAlt + 2) +
      grid.MeshCoef1s3rdp3.slice(iAlt) % value.slice(iAlt + 3) +
      grid.MeshCoef1s3rdp4.slice(iAlt) % value.slice(iAlt + 4) +
      grid.MeshCoef1s3rdp5.slice(iAlt) % value.slice(iAlt + 5);
  */
  gradient = (value.slice(iAlt + 2) - value.slice(iAlt + 1)) /
             grid.dk_edge_m.slice(iAlt + 2);

  valueOut = value.slice(iAlt + 1) - gradient % grid.dk_edge_m.slice(iAlt + 1);
  return valueOut;
}

// --------------------------------------------------------------------------
// Calculate the gradient on the dipole grid
// - This is identical to the spherical grid, except the k/alt direction.
// --------------------------------------------------------------------------
std::vector<arma_cube> calc_gradient_dipole(arma_cube value_scgc, Grid grid) {

  std::vector<arma_cube> gradient_vcgc;

  report.print(3, "Calculating dipole griadient");

  report.print(4, "Going into calc_gradient_i (dipole)");
  gradient_vcgc.push_back(calc_gradient2o_i(value_scgc, grid));


  report.print(4, "Going into calc_gradient_j (dipole)");
  gradient_vcgc.push_back(calc_gradient2o_j(value_scgc, grid));


  report.print(4, "Going into calc_gradient_K (DIPOLE)");
  gradient_vcgc.push_back(calc_gradient2o_k(value_scgc, grid));

  return gradient_vcgc;
}

// --------------------------------------------------------------------------
// Calculate the gradient in cubesphere spatial discretization
// --------------------------------------------------------------------------
std::vector<arma_cube> calc_gradient_cubesphere(arma_cube value, Grid &grid) {
  // Must be used for cubesphere (Probably need a boolean check)
  int64_t nXs = grid.get_nY();
  int64_t nYs = grid.get_nX();
  int64_t nGCs = grid.get_nGCs();
  int64_t nAlts = grid.get_nAlts();

  // Initialize two arma cubes for return
  arma_cube grad_lon(nXs, nYs, nAlts);
  arma_cube grad_lat(nXs, nYs, nAlts);

  for (int64_t iAlt = 0; iAlt < nAlts; iAlt++) {
    /** Extract Grid Features **/
    // Addition: Get a copy of dx dy
    arma_mat curr_refx = grid.refx_scgc.slice(iAlt);
    arma_mat curr_refy = grid.refy_scgc.slice(iAlt);

    // Get some dx dy metrics from the grid
    precision_t dx = grid.drefx(iAlt);
    precision_t dy = grid.drefy(iAlt);

    // Get values of current level
    arma_mat curr_value = value.slice(iAlt);

    /** Calculate Gradient with Central Difference Scheme **/
    // Since Reference grid is orthogonal, we only need gradient along reference xy direction
    // Then we convert gradient from xy direction ot lat-lon direction

    arma_mat grad_x_curr(nXs, nYs);
    arma_mat grad_y_curr(nXs, nYs);

    // Calc gradient on x and y direction (since reference grid)
    // Only update interior cells
    // May vectorize for future improvements

    // if more than 1 nGCs, we do fourth order, some foolproofing in case we go into debug hell
    if (nGCs >= 2) {
      for (int j = nGCs; j < nYs - nGCs; j++) {
        for (int i = nGCs; i < nXs - nGCs; i++) {
          grad_x_curr(i, j) = (-curr_value(i + 2, j) +
                               8 * curr_value(i + 1, j) -
                               8 * curr_value(i - 1, j) +
                               curr_value(i - 2, j)) * (1. / 12. / dx);
          grad_y_curr(i, j) = (-curr_value(i, j + 2) +
                               8 * curr_value(i, j + 1) -
                               8 * curr_value(i, j - 1) +
                               curr_value(i, j - 2)) * (1. / 12. / dy);
        }
      }
    } else { // otherwise we do second order
      for (int j = nGCs; j < nYs - nGCs; j++) {
        for (int i = nGCs; i < nXs - nGCs; i++) {
          grad_x_curr(i, j) = (curr_value(i + 1, j) -
                               curr_value(i - 1, j)) * (1. / 2. / dx);
          grad_y_curr(i, j) = (curr_value(i, j + 1) -
                               curr_value(i, j - 1)) * (1. / 2. / dy);
        }
      }
    }

    // We then use A transformation matrices to convert grad_xy to grad_latlon
    // Ref -> Physical, we use A matrix
    grad_lon.slice(iAlt) = grad_x_curr % grid.A11_inv_scgc.slice(
                             iAlt) + grad_y_curr % grid.A21_inv_scgc.slice(iAlt);
    grad_lat.slice(iAlt) = grad_x_curr % grid.A12_inv_scgc.slice(
                             iAlt) + grad_y_curr % grid.A22_inv_scgc.slice(iAlt);
  }

  // Not initializing with array like procedure in case I get bugs
  std::vector<arma_cube> gradient;
  gradient.push_back(grad_lon);
  gradient.push_back(grad_lat);
  gradient.push_back(calc_gradient_alt(value, grid));

  return gradient;
}