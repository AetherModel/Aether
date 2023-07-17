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

    // This part extracts geogrid features and output
    {   
        output(gGrid.geoLon_scgc.slice(0), "lon.txt", false);
        output(gGrid.geoLat_scgc.slice(0), "lat.txt", false);
        output(gGrid.refx_scgc.slice(0), "x.txt", false);
        output(gGrid.refy_scgc.slice(0), "y.txt", false);
        std::cout << "Side" << std::endl;
        std::cout << quadtree.iSide << std::endl;
    }

    report.exit(function);
    report.times();

  } catch (std::string error) {
    if (iProc == 0) {
      std::cout << error << "\n";
      std::cout << "---- Must Exit! ----\n";
    }
  }

    
  // End parallel tasks:
  iErr = MPI_Finalize();

  return iErr;
}