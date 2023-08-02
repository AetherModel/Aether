// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "../include/aether.h"

/**
 * Output function
 *
 * @param values Values
 * @param filename FileName
 * @param DoAppend
 */
void output(arma_mat &values,
            std::string filename,
            bool DoAppend)
{

  std::ofstream outfile;
  if (DoAppend)
    outfile.open(filename, std::ios_base::app);
  else
  {
    outfile.open(filename);
    int64_t nX = values.n_rows;
    int64_t nY = values.n_cols;
    outfile << nX << " " << nY << "\n";
  }
  outfile << values;
  outfile << "----";
  outfile.close();
}

/**
 * x and y behavior Checker 
 *
 * xs have same cols, ys have same rows
 *
 * @param x x coordinates
 * @param y y coordinates
 * @param tol tolerance
 * 
 */
bool xy_behavior_check(arma_mat &x, arma_mat &y, precision_t &tol)
{
  // Check x has same cols, y has same rows
  // Extract n_rows and n_cols
  int nx = x.n_rows;
  int ny = x.n_cols;

  // Test if y has same rows
  for (int i = 0; i < nx - 1; i++) {
    Row<precision_t> vec1 = y.row(i);
    Row<precision_t> vec2 = y.row(i+1);
    if (!is_approx_equal(vec1, vec2, tol)) {
      return false;
    }
  }
  
  // Test if x has same cols
  for (int i = 0; i < ny - 1; i++) {
    arma_vec vec1 = x.col(i);
    arma_vec vec2 = x.col(i+1);
    if (!is_approx_equal(vec1, vec2, tol)) {
      return false;
    }
  }
  return true;
}

/**
 * Equidistance Checker 
 *
 * @param x x coordinates
 * @param y y coordinates
 * @param tol tolerance
 * 
 */
bool equidistance_check(arma_mat &x, arma_mat &y, precision_t tol)
{
  // Only check first row/col of x and y since they are verified
  // to be the same
  arma_vec x_vec = x.col(0);
  arma_vec y_vec = conv_to<arma_vec>::from(y.row(0));

  int64_t x_size = x_vec.size();
  int64_t y_size = y_vec.size();

  // Check x
  // Trim first and last element to form two displaced vectors
  // for finding the distance
  arma_vec x_vec_left = x_vec.subvec(1,x_size-1);
  arma_vec x_vec_right = x_vec.subvec(0,x_size-2);

  arma_vec x_dist = x_vec_left - x_vec_right;

  if (!is_approx_constant(x_dist, tol)) {
    return false;
  }
  
  // Check y
  // Trim first and last element to form two displaced vectors
  // for finding the distance
  arma_vec y_vec_left = y_vec.subvec(1,y_size-1);
  arma_vec y_vec_right = y_vec.subvec(0,y_size-2);

  arma_vec y_dist = y_vec_left - y_vec_right;

   if (!is_approx_constant(y_dist, tol)) {
    return false;
  }

  return true;
}

/**
 * Range Checker 
 * 
 * MUST select edge coordinates as the function evals edge range
 *
 * @param x x coordinates
 * @param y y coordinates
 * @param tol tolerance
 * @param nGC number of ghost cells
 * 
 */
bool range_check(arma_mat &x, arma_mat &y, precision_t R, precision_t tol, int64_t nGCs) {
  // Check reference coordinate first
  arma_vec x_vec_nGCs = x.col(0)/R;
  arma_vec y_vec_nGCs = conv_to<arma_vec>::from(y.row(0))/R;

  int64_t x_size = x_vec_nGCs.size();
  int64_t y_size = y_vec_nGCs.size();

  // Trim nGCs and normalize
  arma_vec x_vec = x_vec_nGCs.subvec(nGCs,x_size-nGCs-1);
  arma_vec y_vec = y_vec_nGCs.subvec(nGCs,y_size-nGCs-1);

  // First should be -a, Last should be +a. 
  precision_t a = 1/sqrt(3);
  precision_t x_vec_head = x_vec(0);
  precision_t x_vec_tail = x_vec(x_vec.size()-1);
  precision_t y_vec_head = y_vec(0);
  precision_t y_vec_tail = y_vec(y_vec.size()-1);
  
  if (std::abs(x_vec_head+a) > tol) {
    return false;
  }
  if (std::abs(x_vec_tail-a) > tol) {
    return false;
  }

  if (std::abs(y_vec_head+a) > tol) {
    return false;
  }
  if (std::abs(y_vec_tail-a) > tol) {
    return false;
  }

  return true; 
}


int main() {

  int iErr = 0;
  std::string sError;
  bool DidWork = true;

  Times time;
  Report report;

  // Define the function and report:
  std::string function = "main";
  static int iFunction = -1;
  report.enter(function, iFunction);

  try {
  
    // Create inputs (reading the input file):
    Inputs input(time, report);
    if (!input.is_ok())
      throw std::string("input initialization failed!");
    
    Quadtree quadtree(input, report);
    if (!quadtree.is_ok())
      throw std::string("quadtree initialization failed!");
    
    // Initialize MPI and parallel aspects of the code:
    DidWork = init_parallel(input, quadtree, report);
    if (!DidWork)
      throw std::string("init_parallel failed!");
    
    // Initialize the planet:
    Planets planet(input, report);
    MPI_Barrier(aether_comm);
    if (!planet.is_ok())
      throw std::string("planet initialization failed!");

    // Initialize Geographic grid:
    Grid gGrid(input.get_nLonsGeo(),
	           input.get_nLatsGeo(),
	           input.get_nAltsGeo(),
	           nGeoGhosts);
    DidWork = gGrid.init_geo_grid(quadtree, planet, input, report);

    // Coordinate Generation Testing
    {
      // Set tolerance limit 
      precision_t tol = 1e-6;

      // Get number of ghost cells
      int nGCs_lcl = gGrid.get_nGCs();

      // Print current side number
      std::string side_num = std::to_string(quadtree.iSide + 1);
      std::cout << "Initiating Tests for Side Number (1-based index): " << side_num << std::endl;

      /**
        * Extract some test data generated by Aether Model
        * First slice of the cube should suffice as all coordinates are
        * generated the same way
        */
      
      // Grab radius of the altitude 0
      precision_t planet_R = planet.get_radius(0);
      arma_vec R_Alts = gGrid.geoAlt_scgc.tube(0, 0) + planet_R;

      precision_t R = R_Alts(0);

      // Cell center coordinates 
      arma_mat aether_lon_cc = gGrid.geoLon_scgc.slice(0);
      arma_mat aether_lat_cc = gGrid.geoLat_scgc.slice(0);
      arma_mat aether_x_cc = gGrid.refx_scgc.slice(0);
      arma_mat aether_y_cc = gGrid.refy_scgc.slice(0);
      arma_mat aether_A11_cc = gGrid.A11_scgc.slice(0);
      arma_mat aether_A12_cc = gGrid.A12_scgc.slice(0);
      arma_mat aether_A21_cc = gGrid.A21_scgc.slice(0);
      arma_mat aether_A22_cc = gGrid.A22_scgc.slice(0);
      arma_mat aether_A11_inv_cc = gGrid.A11_inv_scgc.slice(0);
      arma_mat aether_A12_inv_cc = gGrid.A12_inv_scgc.slice(0);
      arma_mat aether_A21_inv_cc = gGrid.A21_inv_scgc.slice(0);
      arma_mat aether_A22_inv_cc = gGrid.A22_inv_scgc.slice(0);
      arma_mat aether_g11_upper_cc = gGrid.g11_upper_scgc.slice(0);
      arma_mat aether_g12_upper_cc = gGrid.g12_upper_scgc.slice(0);
      arma_mat aether_g21_upper_cc = gGrid.g21_upper_scgc.slice(0);
      arma_mat aether_g22_upper_cc = gGrid.g22_upper_scgc.slice(0);

      // X Edge coordinates 
      arma_mat aether_lon_xedge = gGrid.geoLon_Left.slice(0);
      arma_mat aether_lat_xedge = gGrid.geoLat_Left.slice(0);
      arma_mat aether_x_xedge = gGrid.refx_Left.slice(0);
      arma_mat aether_y_xedge = gGrid.refy_Left.slice(0);
      arma_mat aether_A11_xedge = gGrid.A11_Left.slice(0);
      arma_mat aether_A12_xedge = gGrid.A12_Left.slice(0);
      arma_mat aether_A21_xedge = gGrid.A21_Left.slice(0);
      arma_mat aether_A22_xedge = gGrid.A22_Left.slice(0);
      arma_mat aether_A11_inv_xedge = gGrid.A11_inv_Left.slice(0);
      arma_mat aether_A12_inv_xedge = gGrid.A12_inv_Left.slice(0);
      arma_mat aether_A21_inv_xedge = gGrid.A21_inv_Left.slice(0);
      arma_mat aether_A22_inv_xedge = gGrid.A22_inv_Left.slice(0);
      arma_mat aether_g11_upper_xedge = gGrid.g11_upper_Left.slice(0);
      arma_mat aether_g12_upper_xedge = gGrid.g12_upper_Left.slice(0);
      arma_mat aether_g21_upper_xedge = gGrid.g21_upper_Left.slice(0);
      arma_mat aether_g22_upper_xedge = gGrid.g22_upper_Left.slice(0);

      // Y Edge coordinates 
      arma_mat aether_lon_yedge = gGrid.geoLon_Down.slice(0);
      arma_mat aether_lat_yedge = gGrid.geoLat_Down.slice(0);
      arma_mat aether_x_yedge = gGrid.refx_Down.slice(0);
      arma_mat aether_y_yedge = gGrid.refy_Down.slice(0);
      arma_mat aether_A11_yedge = gGrid.A11_Down.slice(0);
      arma_mat aether_A12_yedge = gGrid.A12_Down.slice(0);
      arma_mat aether_A21_yedge = gGrid.A21_Down.slice(0);
      arma_mat aether_A22_yedge = gGrid.A22_Down.slice(0);
      arma_mat aether_A11_inv_yedge = gGrid.A11_inv_Down.slice(0);
      arma_mat aether_A12_inv_yedge = gGrid.A12_inv_Down.slice(0);
      arma_mat aether_A21_inv_yedge = gGrid.A21_inv_Down.slice(0);
      arma_mat aether_A22_inv_yedge = gGrid.A22_inv_Down.slice(0);
      arma_mat aether_g11_upper_yedge = gGrid.g11_upper_Down.slice(0);
      arma_mat aether_g12_upper_yedge = gGrid.g12_upper_Down.slice(0);
      arma_mat aether_g21_upper_yedge = gGrid.g21_upper_Down.slice(0);
      arma_mat aether_g22_upper_yedge = gGrid.g22_upper_Down.slice(0);

      /**
      if (stoi(side_num) == 6) {
        output(aether_x_cc, "x_cc6.txt", false);
        output(aether_y_cc, "y_cc6.txt", false);
        output(aether_x_xedge, "x_xedge6.txt", false);
        output(aether_y_xedge, "y_xedge6.txt", false);
        output(aether_x_yedge, "x_yedge6.txt", false);
        output(aether_y_yedge, "y_yedge6.txt", false);
        output(aether_lat_cc, "lat_cc6.txt", false);
        output(aether_lon_cc, "lon_cc6.txt", false);
        output(aether_lat_xedge, "lat_xedge6.txt", false);
        output(aether_lon_xedge, "lon_xedge6.txt", false);
        output(aether_lat_yedge, "lat_yedge6.txt", false);
        output(aether_lon_yedge, "lon_yedge6.txt", false);
      }

      if (stoi(side_num) == 1) {
        output(aether_x_cc, "x_cc1.txt", false);
        output(aether_y_cc, "y_cc1.txt", false);
        output(aether_x_xedge, "x_xedge1.txt", false);
        output(aether_y_xedge, "y_xedge1.txt", false);
        output(aether_x_yedge, "x_yedge1.txt", false);
        output(aether_y_yedge, "y_yedge1.txt", false);
        output(aether_lat_cc, "lat_cc1.txt", false);
        output(aether_lon_cc, "lon_cc1.txt", false);
        output(aether_lat_xedge, "lat_xedge1.txt", false);
        output(aether_lon_xedge, "lon_xedge1.txt", false);
        output(aether_lat_yedge, "lat_yedge1.txt", false);
        output(aether_lon_yedge, "lon_yedge1.txt", false);
      }
      */
  
      // Test 0: Reference Coordinate behavior check
      // rows of y should be the same, cols of x should be the same
      if (!xy_behavior_check(aether_x_cc, aether_y_cc, tol)) {
        std::string err_msg = "XY Behavior Check Failed for Side " + side_num + " , Cell Centers";
        throw std::string(err_msg);
      }
      if (!xy_behavior_check(aether_x_xedge, aether_y_xedge, tol)) {
        std::string err_msg = "XY Behavior Check Failed for Side " + side_num + " , X edges";
        throw std::string(err_msg);
      }
      if (!xy_behavior_check(aether_x_yedge, aether_y_yedge, tol)) {
        std::string err_msg = "XY Behavior Check Failed for Side " + side_num + " , Y edges";
        throw std::string(err_msg);
      }

      // Test 1: Reference Coordinates are equidistant
      if (!equidistance_check(aether_x_cc, aether_y_cc, tol)) {
        std::string err_msg = "Equidistance Check Failed for Side " + side_num + " , Cell Centers";
        throw std::string(err_msg);
      }
      if (!equidistance_check(aether_x_xedge, aether_y_xedge, tol)) {
        std::string err_msg = "Equidistance Check Failed for Side " + side_num + " , X edges";
        throw std::string(err_msg);
      }
      if (!equidistance_check(aether_x_yedge, aether_y_yedge, tol)) {
        std::string err_msg = "Equidistance Check Failed for Side " + side_num + " , Y edges";
        throw std::string(err_msg);
      }

      // Test 2: Range Check
      if (!range_check(aether_x_xedge, aether_y_yedge, R, tol, nGCs_lcl)) {
        std::string err_msg = "Range Check Failed for Side " + side_num;
        throw std::string(err_msg);
      }

      std::cout << "All tests passed for side " + side_num << std::endl;
    }
    report.exit(function);
    report.times();

  } catch (std::string error) {
    //if (iProc == 0) {
      std::cout << error << "\n";
      std::cout << "---- Must Exit! ----\n";
    //}
  }
  MPI_Barrier(aether_comm);

    
  // End parallel tasks:
  iErr = MPI_Finalize();

  return iErr;
}