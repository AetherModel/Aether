
// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>

#include "aether.h"


// Modularize the test function so it's easy to change
std::vector<arma_cube> test_func(Grid grid, Planets planet, bool debug) {

  // one element for each coord; compatibility w/ doing individual functions
  // out_vals has 6 elements:
  // - first 3 are the function & last 3 are expected gradient
  std::vector<arma_cube> out_vals;
  arma_cube one_elem, i_coords, j_coords, k_coords;

  if (grid.IsLatLonGrid) {
    i_coords = grid.geoLon_scgc;
    j_coords = grid.geoLat_scgc;
    k_coords = grid.radius_scgc;

    // use the func cos(i) * sin(j) * r^2
    one_elem = cos(i_coords) % sin(j_coords) % k_coords % k_coords;
    out_vals.push_back(one_elem);
    out_vals.push_back(one_elem);
    out_vals.push_back(one_elem);

    // The true gradient values:
    out_vals.push_back(-600.0 * sin(i_coords) % tan(j_coords) % k_coords);
    out_vals.push_back(cos(i_coords) % cos(j_coords) % k_coords);
    out_vals.push_back(2.0 * cos(i_coords) % sin(j_coords) % k_coords);
  }

  if (grid.IsDipole) {
    std::cout<<"I AMN DIPOLE\n";
    precision_t planetRadius = planet.get_radius(0.0);
    i_coords = grid.magLon_scgc;
    j_coords = grid.magP_scgc;
    k_coords = grid.magQ_scgc;

    // use the func cos(i) * sin(j) * r^2
    // one_elem = cos(grid.magLon_scgc) % sin(grid.magLat_scgc) % grid.radius_scgc % grid.radius_scgc;
    one_elem = i_coords % i_coords + j_coords % j_coords + k_coords % k_coords;
    out_vals.push_back(one_elem);
    out_vals.push_back(one_elem);
    out_vals.push_back(one_elem);

    arma_cube delT = pow(1 + 3.0 * cos(cPI/2. - grid.magLat_scgc) % cos(cPI/2. - grid.magLat_scgc),
                         0.5);

    // arma_cube r = grid.radius_scgc / planetRadius;
    // // The true gradient values:
    // out_vals.push_back(-k_coords % k_coords % sin(i_coords) % j_coords % j_coords
    //                    / (r % pow(cos(grid.magLat_scgc), 2.0))); // mayB sin
    // out_vals.push_back(2.0 * delT % j_coords % j_coords % cos(i_coords) % j_coords
    //                    / pow(cos(grid.magLat_scgc), 3.0)); // mayb sin???
    // out_vals.push_back(2.0 * delT % k_coords % cos(i_coords) % j_coords % j_coords
    //                    / pow(r, 3.0));

    out_vals.push_back(2.0 * i_coords
                       / (grid.radius_scgc % cos(grid.magLat_scgc))); // mayB sin
    out_vals.push_back(2.0 * j_coords
                       % delT / pow(cos(grid.magLat_scgc), 3.0)); // mayb sin???
    out_vals.push_back(2.0 * k_coords
                       % delT / pow(grid.radius_scgc, 3.0));

  }

  if (debug) {
    std::string numproc = tostr(iProc, 2);
    std::string gridshape = "gridshape-" + tostr(grid.iGridShape_, 2);
    i_coords.save(gridshape + "_proc-" + numproc + "_i_center.txt", arma_ascii);
    j_coords.save(gridshape + "_proc-" + numproc + "_j_center.txt", arma_ascii);
    k_coords.save(gridshape + "_proc-" + numproc + "_k_center.txt", arma_ascii);
  }

  return out_vals;
}

bool test_gradient(Planets planet, Quadtree quadtree, json test_config,
                   Grid gGrid, Grid mGrid) {
  std::string function = "test_gradient";
  static int iFunction = -1;
  report.enter(function, iFunction);

  bool didWork = true;
  bool debug = test_config["dump_debug_cubes"];

  report.print(2, "Testing neutral grid");

  if (gGrid.IsDipole || gGrid.IsLatLonGrid)
    didWork = didWork && test_gradient_ijk(planet, gGrid, debug);
  else
    report.error("Cubesphere gradient test not built yet sorry");

  if (!didWork && test_config["exit_on_fail"])
    throw std::string("Gradient test failed - neutral grid");

  report.print(2, "Testing ion grid");

  if (mGrid.IsCubeSphereGrid) // it's technically possible...
    didWork = didWork && test_gradient_cubesphere(planet, quadtree, mGrid);

  if (mGrid.IsDipole || mGrid.IsLatLonGrid)
    didWork = didWork && test_gradient_ijk(planet, mGrid, debug);

  // if (!didWork && test_config["exit_on_fail"])
  //   throw std::string("Gradient test failed - ion grid");

  report.exit(function);

  return didWork;
}

void send_message(std::string Message, int nGood, int nBad) {
  std::string newMessage;
  newMessage = "iProc: " + tostr(iProc, 2) + " " + Message;
  printf("%s has FAILED! (%i/%i); or (%f) perc\n", newMessage.data(), nBad, nGood,
         100.*nBad / nGood);
  return;
}


bool test_gradient_ijk(Planets planet, Grid grid, bool debug) {

  std::string function = "test_gradient_ijk";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int64_t nX = grid.get_nX();
  int64_t nY = grid.get_nY();
  int64_t nZ = grid.get_nZ();
  int64_t nGCs = grid.get_nGCs();

  // numbers of grid points without ghost cells:
  int64_t nI, nJ, nK;
  nI = nX - 2 * nGCs;
  nJ = nY - 2 * nGCs;
  nK = nZ - 2 * nGCs;

  bool didWork = true;

  int64_t nCellsTot = nX * nY * nZ;
  int64_t nCellsNGCs = nI * nJ * nK;

  report.print(2, "Beginning gradient test");

  std::vector<arma_cube> tmp, func_values, true_gradient, predicted_gradient;
  tmp = test_func(grid, planet, debug);

  func_values.push_back(tmp[0]);
  func_values.push_back(tmp[1]);
  func_values.push_back(tmp[2]);
  true_gradient.push_back(tmp[3]);
  true_gradient.push_back(tmp[4]);
  true_gradient.push_back(tmp[5]);

  if (grid.IsDipole) {
    predicted_gradient.push_back(calc_gradient2o_i(func_values[0], grid));
    predicted_gradient.push_back(calc_gradient2o_j(func_values[1], grid));
    predicted_gradient.push_back(calc_gradient2o_k(func_values[2], grid));
  } else {
    predicted_gradient.push_back(calc_gradient2o_i(func_values[0], grid));
    predicted_gradient.push_back(calc_gradient2o_j(func_values[1], grid));
    predicted_gradient.push_back(calc_gradient2o_k(func_values[2], grid));
  }

  arma::uvec bad_is, bad_js, bad_ks;

  // Look for values > 5% different from expected
  bad_is = find(abs(
                  (predicted_gradient[0].subcube(2, 2, 2, nI + 2, nJ + 2, nK + 2)
                   - true_gradient[0].subcube(2, 2, 2, nI + 2, nJ + 2, nK + 2))
                  / true_gradient[0].subcube(2, 2, 2, nI + 2, nJ + 2, nK + 2)) > 0.25);
  bad_js = find(abs(
                  (predicted_gradient[1].subcube(2, 2, 2, nI + 2, nJ + 2, nK + 2)
                   - true_gradient[1].subcube(2, 2, 2, nI + 2, nJ + 2, nK + 2))
                  / true_gradient[1].subcube(2, 2, 2, nI + 2, nJ + 2, nK + 2)) > 0.25);
  bad_ks = find(abs(
                  (predicted_gradient[2].subcube(2, 2, 2, nI + 2, nJ + 2, nK + 2)
                   - true_gradient[2].subcube(2, 2, 2, nI + 2, nJ + 2, nK + 2))
                  / true_gradient[2].subcube(2, 2, 2, nI + 2, nJ + 2, nK + 2)) > 0.25);

  // ghost cells are hard; if more than 1% of *real* cells are out of spec, the test fails.
  if (bad_is.n_elem > 0.1 * nCellsNGCs) {
    send_message("grad_i:", nCellsNGCs, bad_is.n_elem);
    didWork = false;
  }

  if (bad_js.n_elem > 0.1 * nCellsNGCs) {
    send_message("grad_j:", nCellsNGCs, bad_js.n_elem);
    didWork = false;
  }

  if (bad_ks.n_elem > 0.1 * nCellsNGCs) {
    send_message("grad_k:", nCellsNGCs, bad_ks.n_elem);
    didWork = false;
  }

  // Output if requested:
  std::string numproc = tostr(iProc, 2);

  if (debug) {
    std::string numproc = tostr(iProc, 2);
    std::string gridshape = "gridshape-" + tostr(grid.iGridShape_, 2) + "_iproc-";
    grid.di_center_m_scgc.save(gridshape + numproc + "_di_center_m.txt",
                               arma_ascii);
    grid.dj_center_m_scgc.save(gridshape + numproc + "_dj_center_m.txt",
                               arma_ascii);
    grid.dk_center_m_scgc.save(gridshape + numproc + "_dk_center_m.txt",
                               arma_ascii);

    func_values[0].save(gridshape + numproc + "_testfunc.txt", arma_ascii);
    true_gradient[0].save(gridshape + numproc + "_actual_grad_i.txt",
                          arma_ascii);
    true_gradient[1].save(gridshape + numproc + "_actual_grad_j.txt",
                          arma_ascii);
    true_gradient[2].save(gridshape + numproc + "_actual_grad_k.txt",
                          arma_ascii);

    predicted_gradient[0].save(gridshape + numproc + "_i-predicted-grad.txt",
                               arma_ascii);
    predicted_gradient[1].save(gridshape + numproc + "_j-predicted-grad.txt",
                               arma_ascii);
    predicted_gradient[2].save(gridshape + numproc + "_k-predicted-grad.txt",
                               arma_ascii);
  }

  // For completeness, check for non-finites
  didWork = didWork && all_finite(true_gradient, "TRUE GRADIENT");
  didWork = didWork && all_finite(func_values, "FUNCTION");
  didWork = didWork && all_finite(predicted_gradient, "AETHER'S GRADIENT");

  report.report_errors();

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