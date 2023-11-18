// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "../include/aether.h"

/**
 * Output function for debugging
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

int main() {

  int iErr = 0;
  std::string sError;
  bool DidWork = true;

  Times time;

  // Define the function and report:
  std::string function = "main";
  static int iFunction = -1;
  report.enter(function, iFunction);

  try {
  
    // Create inputs (reading the input file):
    input = Inputs(time);
    if (!input.is_ok())
      throw std::string("input initialization failed!");
    
    Quadtree quadtree;
    if (!quadtree.is_ok())
      throw std::string("quadtree initialization failed!");
    
    // Initialize MPI and parallel aspects of the code:
    DidWork = init_parallel(quadtree);
    if (!DidWork)
      throw std::string("init_parallel failed!");
    
    // Initialize the planet:
    Planets planet;
    MPI_Barrier(aether_comm);
    if (!planet.is_ok())
      throw std::string("planet initialization failed!");

    // Initialize Geographic grid:
    Grid gGrid(input.get_nLonsGeo(),
	           input.get_nLatsGeo(),
	           input.get_nAltsGeo(),
	           nGeoGhosts);
    DidWork = gGrid.init_geo_grid(quadtree, planet);

    // First check whether the initialization uses exactly 6 processes. 
    // The exactly 6 requirements is due to the checking of the range of reference coordinate system
    int world_size;
    MPI_Comm_size(aether_comm, &world_size);
    if (world_size != 6) {
      throw std::string("Comm size must be 6!!!");
    }

    // Gradient Test
    {
      // Set tolerance limit 
      precision_t tol = 1e-5;

      // Print current side number
      std::string side_num = std::to_string(quadtree.iSide + 1);
      std::cout << "Initiating Test 1 for Side Number (1-based index): " << side_num << std::endl;

      /**
        * Extract some test data generated by Aether Model
        */
      
      // Cell center coordinates 
      arma_mat aether_lon_cc = gGrid.geoLon_scgc.slice(0);
      arma_mat aether_lat_cc = gGrid.geoLat_scgc.slice(0);

      int64_t nXs = gGrid.get_nY();
      int64_t nYs = gGrid.get_nX();
      int64_t nGCs = gGrid.get_nGCs();
      int64_t nAlts = gGrid.get_nAlts();

      // Test scalar field and gradients
      arma_cube scgc(nXs, nYs, nAlts);
      arma_cube grad_lon_analytical(nXs, nYs, nAlts);
      arma_cube grad_lat_analytical(nXs, nYs, nAlts);

      // Radius Information
      precision_t planet_R = planet.get_radius(0);
      // radius of planet + altitude
      // just pick alt at (0,0) loction
      arma_vec R_Alts = gGrid.geoAlt_scgc.tube(0, 0) + planet_R;

      for (int iAlt = 0; iAlt < nAlts; iAlt++) {
        arma_mat curr_scalar(nXs, nYs, arma::fill::zeros); // setup zero mat
        arma_mat curr_grad_lon(nXs, nYs);
        arma_mat curr_grad_lat(nXs, nYs);
        precision_t A = 1;
        precision_t B = 1;

        for (int j = 0; j < nYs; j++)
        {
          for (int i = 0; i < nXs; i++)
          {
            precision_t curr_lat = aether_lat_cc(i, j);
            precision_t curr_lon = aether_lon_cc(i, j);

            curr_scalar(i,j) = std::sin(curr_lat);
            curr_grad_lon(i,j) = 0.;
            curr_grad_lat(i,j) = std::cos(curr_lat); // Assume R=1, we will scale the numerical result
          }
        }
        scgc.slice(iAlt) = curr_scalar;
        grad_lon_analytical.slice(iAlt) = curr_grad_lon;
        grad_lat_analytical.slice(iAlt) = curr_grad_lat;
      }

      std::vector<arma_cube> test_res = calc_gradient_cubesphere(scgc, gGrid);

      // Perform Tests
      for (int iAlt = 0; iAlt < nAlts; iAlt++) {
        arma_mat curr_grad_lon = grad_lon_analytical.slice(iAlt);
        arma_mat curr_grad_lat = grad_lat_analytical.slice(iAlt);
        arma_mat curr_numgrad_lon = test_res[0].slice(iAlt);
        arma_mat curr_numgrad_lat = test_res[1].slice(iAlt);
        

        // Evaluate actual cells only
        for (int j = nGCs; j < nYs - nGCs; j++)
        {
          for (int i = nGCs; i < nXs - nGCs; i++)
          {
            if (std::abs(curr_grad_lat(i,j) - curr_numgrad_lat(i,j) * R_Alts(iAlt)) > 1e-4) { // For float precision
              std::cout << "Found Incorrect latitudinal gradient for face " + side_num + ", test f = sin(lat)" << std::endl;
              std::cout << std::abs(curr_grad_lat(i,j) - curr_numgrad_lat(i,j)* R_Alts(iAlt)) << std::endl;
              std::cout << iAlt << std::endl;
              goto endloop1;
            }

            if (std::abs(curr_grad_lon(i,j) - curr_numgrad_lon(i,j) * R_Alts(iAlt)) > 1e-4) { // For float precision
              std::cout << "Found Incorrect longitudinal gradient for face " + side_num + ", test f = sin(lat)" << std::endl;
              goto endloop1;
            }
          }
        }
      }
    }
    endloop1:

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