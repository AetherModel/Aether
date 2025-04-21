
// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>

#include "aether.h"

bool test_gradient(Planets planet, Quadtree quadtree, json test_config,
                   Grid gGrid, Grid mGrid) {
  std::string function = "test_gradient";
  static int iFunction = -1;
  report.enter(function, iFunction);

  bool didWork;

  report.print(2, "Testing neutral grid");

  if (gGrid.IsCubeSphereGrid)
    didWork = test_gradient_cubesphere(planet, quadtree, gGrid);

  if (gGrid.IsDipole || gGrid.IsLatLonGrid)
    didWork = test_gradient_ijk(planet, gGrid);

  MPI_Barrier(aether_comm);

  if (!didWork && test_config["exit_on_fail"])
    throw std::string("Gradient test failed - neutral grid");

  report.print(2, "Testing ion grid");

  if (mGrid.IsCubeSphereGrid) // it's technically possible...
    didWork = test_gradient_cubesphere(planet, quadtree, mGrid);

  if (mGrid.IsDipole || mGrid.IsLatLonGrid)
    didWork = test_gradient_ijk(planet, mGrid);

  if (!didWork && test_config["exit_on_fail"])
    throw std::string("Gradient test failed - ion grid");



  report.exit(function);

  return didWork;
}

bool test_gradient_ijk(Planets planet, Grid grid) {

  std::string function = "test_gradient_dipole";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int64_t nIs = grid.get_nX();
  int64_t nJs = grid.get_nY();
  int64_t nKs = grid.get_nZ();
  int64_t nGCs = grid.get_nGCs();

  int64_t nX, nY, nZ;
  precision_t tol = 1e-3;
  bool didWork = true;

  arma_cube predicted_gradient, true_gradient;
  arma::uvec err_points;

  arma_cube gradient_error;
  gradient_error.set_size(nIs, nJs, nKs);
  gradient_error.zeros();

  int64_t nCellsTot = nX * nY * nZ;
  int64_t nCellsNGCs = (nX - 2 * nGCs) * (nY - 2 * nGCs) * (nZ - 2 * nGCs);

  report.print(2, "Beginning i-gradient");

  /////////////////////////////////////////////////////////////
  // Test the gradient in i-direction, d/dx(sin x) = cos(x)  //
  /////////////////////////////////////////////////////////////

  predicted_gradient = calc_gradient2o_i(sin(grid.i_center_scgc), grid);
  true_gradient = cos(grid.i_center_scgc) / (grid.di_center_m_scgc);

  gradient_error = abs(predicted_gradient - true_gradient) / abs(true_gradient);
  err_points = find(abs(gradient_error.subcube(nGCs, nGCs,
                                               nGCs, // don't look at ghost cells
                                               size(nIs - 2 * nGCs, nJs - 2 * nGCs, nKs - 2 * nGCs)))
                    > tol);

  didWork = all_finite(predicted_gradient, "Gradient_4o_i");

  std::cout << "(iproc " << iProc << ", gridtype: " << grid.get_gridtype()
            << ") => Points in i-gradient, above tol: " <<
            100.0 * err_points.n_elem / predicted_gradient.n_elem
            << "% (" << err_points.n_elem << ", " <<
            predicted_gradient.n_elem << ")\n";

  if (err_points.n_elem > true_gradient.n_elem * tol)
    didWork = false;

  report.print(2, "Beginning j-gradient");

  /////////////////////////////////////////////////////////////
  // Test the gradient in j-direction, d/dx(cos x) = -sin(x) //
  /////////////////////////////////////////////////////////////

  predicted_gradient = calc_gradient2o_j(cos(grid.j_center_scgc), grid);
  true_gradient = -1.0 * sin(grid.j_center_scgc) / grid.dj_center_m_scgc;

  gradient_error = predicted_gradient - true_gradient;
  err_points = find(abs(gradient_error.subcube(nGCs, nGCs, nGCs,
                                               size(nIs - 2 * nGCs, nJs - 2 * nGCs, nKs - 2 * nGCs)))
                    > tol);

  didWork = didWork && all_finite(predicted_gradient, "Gradient_2o_j");

  std::cout << "(iproc " << iProc << ", gridtype: " << grid.get_gridtype()
            << ") => Points in j-gradient, above tol: " <<
            100.0 * err_points.n_elem / predicted_gradient.n_elem
            << "% (" << err_points.n_elem << ", " <<
            predicted_gradient.n_elem << ")\n";

  if (err_points.n_elem > true_gradient.n_elem * tol)
    didWork = false;


  report.print(2, "Beginning k-gradient");

  //////////////////////////////////////////////////////
  // Test the gradient in k-direction, d/dx(x^2) = 2x //
  //////////////////////////////////////////////////////

  predicted_gradient = calc_gradient2o_k(grid.radius2_scgc, grid);
  true_gradient = 2.0 * grid.radius_scgc / grid.dk_center_m_scgc;

  gradient_error = (predicted_gradient - true_gradient);
  err_points = find(abs(gradient_error.subcube(nGCs, nGCs, nGCs,
                                               size(nIs - 2 * nGCs, nJs - 2 * nGCs, nKs - 2 * nGCs)))
                    > tol);

  didWork = didWork && all_finite(predicted_gradient, "Gradient_2o_k");

  std::cout << "(iproc " << iProc << ", gridtype: " << grid.get_gridtype()
            << ") => Points in k-gradient, above tol: " <<
            100.0 * err_points.n_elem / predicted_gradient.n_elem
            << "% (" << err_points.n_elem << ", " <<
            predicted_gradient.n_elem << ")\n";

  if (err_points.n_elem > true_gradient.n_elem * tol)
    didWork = false;

  report.exit(function);

  return didWork;
}


// This is non-functional. 
// Taken from src/main/main_test_gradient.cpp with enough edits to compile.
bool test_gradient_cubesphere(Planets planet, Quadtree quadtree, Grid grid) {


  std::string function = "test_gradient_cubesphere";
  static int iFunction = -1;
  report.enter(function, iFunction);

  // Set tolerance limit
  precision_t tol = 1e-5;

  // Print current side number
  std::string side_num = std::to_string(quadtree.iSide + 1);
  std::cout << "Initiating Test 1 for Side Number (1-based index): " << side_num
            << std::endl;

  /**
    * Extract some test data generated by Aether Model
    */

  // Cell center coordinates
  arma_mat aether_lon_cc = grid.geoLon_scgc.slice(0);
  arma_mat aether_lat_cc = grid.geoLat_scgc.slice(0);

  int64_t nXs = grid.get_nY();
  int64_t nYs = grid.get_nX();
  int64_t nGCs = grid.get_nGCs();
  int64_t nAlts = grid.get_nAlts();

  // Test scalar field and gradients
  arma_cube scgc(nXs, nYs, nAlts);
  arma_cube grad_lon_analytical(nXs, nYs, nAlts);
  arma_cube grad_lat_analytical(nXs, nYs, nAlts);

  // Radius Information
  precision_t planet_R = planet.get_radius(0);
  // radius of planet + altitude
  // just pick alt at (0,0) loction
  arma_vec R_Alts = grid.geoAlt_scgc.tube(0, 0) + planet_R;

  for (int iAlt = 0; iAlt < nAlts; iAlt++) {
    arma_mat curr_scalar(nXs, nYs, arma::fill::zeros); // setup zero mat
    arma_mat curr_grad_lon(nXs, nYs);
    arma_mat curr_grad_lat(nXs, nYs);
    precision_t A = 1;
    precision_t B = 1;

    for (int j = 0; j < nYs; j++) {
      for (int i = 0; i < nXs; i++) {
        precision_t curr_lat = aether_lat_cc(i, j);
        precision_t curr_lon = aether_lon_cc(i, j);

        curr_scalar(i, j) = std::sin(curr_lat);
        curr_grad_lon(i, j) = 0.;
        curr_grad_lat(i, j) = std::cos(
                                curr_lat); // Assume R=1, we will scale the numerical result
      }
    }

    scgc.slice(iAlt) = curr_scalar;
    grad_lon_analytical.slice(iAlt) = curr_grad_lon;
    grad_lat_analytical.slice(iAlt) = curr_grad_lat;
  }

  std::vector<arma_cube> test_res = calc_gradient_cubesphere(scgc, grid);

  // Perform Tests
  for (int iAlt = 0; iAlt < nAlts; iAlt++) {
    arma_mat curr_grad_lon = grad_lon_analytical.slice(iAlt);
    arma_mat curr_grad_lat = grad_lat_analytical.slice(iAlt);
    arma_mat curr_numgrad_lon = test_res[0].slice(iAlt);
    arma_mat curr_numgrad_lat = test_res[1].slice(iAlt);


    // Evaluate actual cells only
    for (int j = nGCs; j < nYs - nGCs; j++) {
      for (int i = nGCs; i < nXs - nGCs; i++) {
        if (std::abs(curr_grad_lat(i, j) - curr_numgrad_lat(i,
                                                            j) * R_Alts(iAlt)) > 1e-4) { // For float precision
          std::cout << "Found Incorrect latitudinal gradient for face " + side_num +
                    ", test f = sin(lat)" << std::endl;
          std::cout << std::abs(curr_grad_lat(i, j) - curr_numgrad_lat(i,
                                j)* R_Alts(iAlt)) << std::endl;
          std::cout << iAlt << std::endl;
          goto endloop1;
        }

        if (std::abs(curr_grad_lon(i, j) - curr_numgrad_lon(i,
                                                            j) * R_Alts(iAlt)) > 1e-4) { // For float precision
          std::cout << "Found Incorrect longitudinal gradient for face " + side_num +
                    ", test f = sin(lat)" << std::endl;
          goto endloop1;
        }
      }
    }
  }

endloop1:

  report.exit(function);
  report.times();

  return false;
}