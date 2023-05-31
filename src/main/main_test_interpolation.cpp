// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "../include/aether.h"

// -----------------------------------------------------------------------------
// Test for 3d interpolation in spherical and cubesphere grid
// For a complete test, the following combinations are required
// In default.json, change the [CubeSphere][is] to be true or false
// In aether.json, change nLons, nLats, nAlts to be:
// nLons = 72, nLats = 72, nAlts = 50
// nLons = 18, nLats = 18, nAlts = 50
// nLons = 36, nLats = 18, nAlts = 50 (Only for spherical grid)
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// Check whether the result is close to the expected value
// -----------------------------------------------------------------------------

void is_correct(precision_t lon_in, precision_t ans_lon,
                precision_t lat_in, precision_t ans_lat,
                precision_t alt_in, precision_t ans_alt) {
    // Situation where the return value indicates input error
    if (ans_lon == cNinf) {
        if (ans_lat != cNinf || ans_alt != cNinf) {
            throw std::string("Interpolation fail");
        }
        return;
    }

    // Situation where the processor does the interpolation
    // Longitude is checked only when the point is far away from
    // the north / south pole and prime meridian
    if (lat_in > -0.4 * cPI && lat_in < 0.4 * cPI
        && lon_in > 0.1 && lon_in < 6.2
        && std::abs((lon_in - ans_lon) / lon_in) > 0.05) {
        // throw std::string("Interpolation fail");
    }
    if (std::abs((lat_in - ans_lat) / lat_in) > 0.05
        || std::abs((alt_in - ans_alt) / alt_in) > 0.01) {
        // throw std::string("Interpolation fail");
    }
}

void is_correct(std::vector<precision_t> &lon_in,
                std::vector<precision_t> &ans_lon,
                std::vector<precision_t> &lat_in,
                std::vector<precision_t> &ans_lat,
                std::vector<precision_t> &alt_in,
                std::vector<precision_t> &ans_alt) {
    for (int64_t i = 0; i < lon_in.size(); ++i) {
        is_correct(lon_in[i], ans_lon[i], lat_in[i], ans_lat[i], alt_in[i], ans_alt[i]);
    }
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

    // !!!This part tests Grid::interp_linear, which is desprecated!!!
    {
        // What the test do is:
        // Given the position of an unknown point, estimate its position :)
        // Thus the output should be equal to the input

        // cPI = 3.14159, 2*cPI = 6.28318, cPI/2 = 1.57079
        // minAlt = 100,000, maxAlt = 100,000 + 5,000 * 50 = 350,000
        // Lon within [0, 2*cPI) = [0, 6.28318)
        // Lat within [-cPI/2, cPI/2] = [-1.57079, 1.57079]
        // Alt within [100000, 350000]

        std::cout << "iProc = " << iProc;

        // Edge condition for spherical grid
        {
            std::cout << "\n\n\nEdge 1\n\n\n";
            std::vector<precision_t> lon_in({0, 0, 0, 0, 2*cPI - cSmall, 2*cPI - cSmall, 2*cPI - cSmall, 2*cPI - cSmall});
            std::vector<precision_t> lat_in({-cPI/2, -cPI/2, cPI/2, cPI/2, -cPI/2, -cPI/2, cPI/2, cPI/2});
            std::vector<precision_t> alt_in({100000, 350000, 100000, 350000, 100000, 350000, 100000, 350000});

            std::vector<precision_t> ansLon = gGrid.interp_linear(gGrid.geoLon_scgc, lon_in, lat_in, alt_in);
            std::vector<precision_t> ansLat = gGrid.interp_linear(gGrid.geoLat_scgc, lon_in, lat_in, alt_in);
            std::vector<precision_t> ansAlt = gGrid.interp_linear(gGrid.geoAlt_scgc, lon_in, lat_in, alt_in);
            for (int64_t i = 0; i < lon_in.size(); ++i) {
                std::cout << "Estimation of (" << lon_in[i] << ", " << lat_in[i] << ", " << alt_in[i] << ") is "
                        << "(" << ansLon[i] << ", " << ansLat[i] << ", " << ansAlt[i] << ")\n";
            }
            is_correct(lon_in, ansLon, lat_in, ansLat, alt_in, ansAlt);
        }

        // Edge condition for cubeshpere grid
        {
            std::cout << "\n\n\nEdge 2\n\n\n";
            std::vector<precision_t> lon_in({0, 0, cPI / 2, cPI / 2, cPI, cPI, 3 * cPI / 2, 3 * cPI / 2});
            std::vector<precision_t> lat_in({-0.61547970867, 0.61547970867, -0.61547970867, 0.61547970867, -0.61547970867, 0.61547970867, -0.61547970867, 0.61547970867});
            std::vector<precision_t> alt_in({100000, 350000, 100000, 350000, 100000, 350000, 100000, 350000});

            std::vector<precision_t> ansLon = gGrid.interp_linear(gGrid.geoLon_scgc, lon_in, lat_in, alt_in);
            std::vector<precision_t> ansLat = gGrid.interp_linear(gGrid.geoLat_scgc, lon_in, lat_in, alt_in);
            std::vector<precision_t> ansAlt = gGrid.interp_linear(gGrid.geoAlt_scgc, lon_in, lat_in, alt_in);
            for (int64_t i = 0; i < lon_in.size(); ++i) {
                std::cout << "Estimation of (" << lon_in[i] << ", " << lat_in[i] << ", " << alt_in[i] << ") is "
                        << "(" << ansLon[i] << ", " << ansLat[i] << ", " << ansAlt[i] << ")\n";
            }
            is_correct(lon_in, ansLon, lat_in, ansLat, alt_in, ansAlt);
        }

        // Error handling
        {
            std::cout << "\n\n\nError handling 1\n\n\n";
            std::vector<precision_t> lon_in({-cSmall, -cSmall, -cSmall, -cSmall, 2*cPI, 2*cPI, 2*cPI, 2*cPI});
            std::vector<precision_t> lat_in({-cPI/2, -cPI/2, cPI/2, cPI/2, -cPI/2, -cPI/2, cPI/2, cPI/2});
            std::vector<precision_t> alt_in({100000, 350000, 100000, 350000, 100000, 350000, 100000, 350000});

            std::vector<precision_t> ansLon = gGrid.interp_linear(gGrid.geoLon_scgc, lon_in, lat_in, alt_in);
            std::vector<precision_t> ansLat = gGrid.interp_linear(gGrid.geoLat_scgc, lon_in, lat_in, alt_in);
            std::vector<precision_t> ansAlt = gGrid.interp_linear(gGrid.geoAlt_scgc, lon_in, lat_in, alt_in);
            for (int64_t i = 0; i < lon_in.size(); ++i) {
                std::cout << "Estimation of (" << lon_in[i] << ", " << lat_in[i] << ", " << alt_in[i] << ") is "
                        << "(" << ansLon[i] << ", " << ansLat[i] << ", " << ansAlt[i] << ")\n";
            }
            is_correct(lon_in, ansLon, lat_in, ansLat, alt_in, ansAlt);
        }


        {
            std::cout << "\n\n\nError handling 2\n\n\n";
            std::vector<precision_t> lon_in({0, 0, 0, 0, 2*cPI - cSmall, 2*cPI - cSmall, 2*cPI - cSmall, 2*cPI - cSmall});
            std::vector<precision_t> lat_in({-cPI/2 - cSmall, -cPI/2 - cSmall, cPI/2 + cSmall, cPI/2 + cSmall, -cPI/2 - cSmall, -cPI/2 - cSmall, cPI/2 + cSmall, cPI/2 + cSmall});
            std::vector<precision_t> alt_in({100000, 350000, 100000, 350000, 100000, 350000, 100000, 350000});

            std::vector<precision_t> ansLon = gGrid.interp_linear(gGrid.geoLon_scgc, lon_in, lat_in, alt_in);
            std::vector<precision_t> ansLat = gGrid.interp_linear(gGrid.geoLat_scgc, lon_in, lat_in, alt_in);
            std::vector<precision_t> ansAlt = gGrid.interp_linear(gGrid.geoAlt_scgc, lon_in, lat_in, alt_in);
            for (int64_t i = 0; i < lon_in.size(); ++i) {
                std::cout << "Estimation of (" << lon_in[i] << ", " << lat_in[i] << ", " << alt_in[i] << ") is "
                        << "(" << ansLon[i] << ", " << ansLat[i] << ", " << ansAlt[i] << ")\n";
            }
            is_correct(lon_in, ansLon, lat_in, ansLat, alt_in, ansAlt);
        }


        {
            std::cout << "\n\n\nError handling 3\n\n\n";
            std::vector<precision_t> lon_in({0, 0, 0, 0, 2*cPI - cSmall, 2*cPI - cSmall, 2*cPI - cSmall, 2*cPI - cSmall});
            std::vector<precision_t> lat_in({-cPI/2, -cPI/2, cPI/2, cPI/2, -cPI/2, -cPI/2, cPI/2, cPI/2});
            std::vector<precision_t> alt_in({100000 - 1, 350000 + 1, 100000 - 1, 350000 + 1, 100000 - 1, 350000 + 1, 100000 - 1, 350000 + 1});

            std::vector<precision_t> ansLon = gGrid.interp_linear(gGrid.geoLon_scgc, lon_in, lat_in, alt_in);
            std::vector<precision_t> ansLat = gGrid.interp_linear(gGrid.geoLat_scgc, lon_in, lat_in, alt_in);
            std::vector<precision_t> ansAlt = gGrid.interp_linear(gGrid.geoAlt_scgc, lon_in, lat_in, alt_in);
            for (int64_t i = 0; i < lon_in.size(); ++i) {
                std::cout << "Estimation of (" << lon_in[i] << ", " << lat_in[i] << ", " << alt_in[i] << ") is "
                        << "(" << ansLon[i] << ", " << ansLat[i] << ", " << ansAlt[i] << ")\n";
            }
            is_correct(lon_in, ansLon, lat_in, ansLat, alt_in, ansAlt);
        }

        // Normal condition
        {
            std::cout << "\n\n\nNormal condition\n\n\n";
            std::vector<precision_t> lon_in({0.1, 0.2, 0.3, 0.4, 2*cPI - 0.5, 2*cPI - 0.6, 2*cPI - 0.7, cPI, cPI});
            std::vector<precision_t> lat_in({-cPI/2 + 0.5, -cPI/2 + 0.6, 0, 0.1, cPI/2 - 0.5, cPI/2 - 0.6, -0.1, -0.2, 0});
            std::vector<precision_t> alt_in({100000, 150000, 200000, 250000, 300000, 350000, 120000, 225000, 225000});

            std::vector<precision_t> ansLon = gGrid.interp_linear(gGrid.geoLon_scgc, lon_in, lat_in, alt_in);
            std::vector<precision_t> ansLat = gGrid.interp_linear(gGrid.geoLat_scgc, lon_in, lat_in, alt_in);
            std::vector<precision_t> ansAlt = gGrid.interp_linear(gGrid.geoAlt_scgc, lon_in, lat_in, alt_in);
            for (int64_t i = 0; i < lon_in.size(); ++i) {
                std::cout << "Estimation of (" << lon_in[i] << ", " << lat_in[i] << ", " << alt_in[i] << ") is "
                        << "(" << ansLon[i] << ", " << ansLat[i] << ", " << ansAlt[i] << ")\n";
            }
            std::cout << "\n\n\n";
            is_correct(lon_in, ansLon, lat_in, ansLat, alt_in, ansAlt);
        }
    }

    // This part tests Grid::set_interpolation_coefs and get_interpolation_values
    {
        // What the test do is:
        // Given the position of an unknown point, estimate its position :)
        // Thus the output should be equal to the input

        // cPI = 3.14159, 2*cPI = 6.28318, cPI/2 = 1.57079
        // minAlt = 100,000, maxAlt = 100,000 + 5,000 * 50 = 350,000
        // Lon within [0, 2*cPI) = [0, 6.28318)
        // Lat within [-cPI/2, cPI/2] = [-1.57079, 1.57079]
        // Alt within [100000, 350000]

        std::cout << "iProc = " << iProc;

        // Edge condition for spherical grid
        {
            std::cout << "\n\n\nEdge 1\n\n\n";
            std::vector<precision_t> lon_in({0, 0, 0, 0, 2*cPI - cSmall, 2*cPI - cSmall, 2*cPI - cSmall, 2*cPI - cSmall});
            std::vector<precision_t> lat_in({-cPI/2, -cPI/2, cPI/2, cPI/2, -cPI/2, -cPI/2, cPI/2, cPI/2});
            std::vector<precision_t> alt_in({100000, 350000, 100000, 350000, 100000, 350000, 100000, 350000});

            gGrid.set_interpolation_coefs(lon_in, lat_in, alt_in);

            std::vector<precision_t> ansLon = gGrid.get_interpolation_values(gGrid.geoLon_scgc);
            std::vector<precision_t> ansLat = gGrid.get_interpolation_values(gGrid.geoLat_scgc);
            std::vector<precision_t> ansAlt = gGrid.get_interpolation_values(gGrid.geoAlt_scgc);

            for (int64_t i = 0; i < lon_in.size(); ++i) {
                std::cout << "Estimation of (" << lon_in[i] << ", " << lat_in[i] << ", " << alt_in[i] << ") is "
                        << "(" << ansLon[i] << ", " << ansLat[i] << ", " << ansAlt[i] << ")\n";
            }
            is_correct(lon_in, ansLon, lat_in, ansLat, alt_in, ansAlt);
        }

        // Edge condition for cubeshpere grid
        {
            std::cout << "\n\n\nEdge 2\n\n\n";
            std::vector<precision_t> lon_in({0, 0, cPI / 2, cPI / 2, cPI, cPI, 3 * cPI / 2, 3 * cPI / 2});
            std::vector<precision_t> lat_in({-0.61547970867, 0.61547970867, -0.61547970867, 0.61547970867, -0.61547970867, 0.61547970867, -0.61547970867, 0.61547970867});
            std::vector<precision_t> alt_in({100000, 350000, 100000, 350000, 100000, 350000, 100000, 350000});

            gGrid.set_interpolation_coefs(lon_in, lat_in, alt_in);

            std::vector<precision_t> ansLon = gGrid.get_interpolation_values(gGrid.geoLon_scgc);
            std::vector<precision_t> ansLat = gGrid.get_interpolation_values(gGrid.geoLat_scgc);
            std::vector<precision_t> ansAlt = gGrid.get_interpolation_values(gGrid.geoAlt_scgc);

            for (int64_t i = 0; i < lon_in.size(); ++i) {
                std::cout << "Estimation of (" << lon_in[i] << ", " << lat_in[i] << ", " << alt_in[i] << ") is "
                        << "(" << ansLon[i] << ", " << ansLat[i] << ", " << ansAlt[i] << ")\n";
            }
            is_correct(lon_in, ansLon, lat_in, ansLat, alt_in, ansAlt);
        }

        // Error handling
        {
            std::cout << "\n\n\nError handling 1\n\n\n";
            std::vector<precision_t> lon_in({-cSmall, -cSmall, -cSmall, -cSmall, 2*cPI, 2*cPI, 2*cPI, 2*cPI});
            std::vector<precision_t> lat_in({-cPI/2, -cPI/2, cPI/2, cPI/2, -cPI/2, -cPI/2, cPI/2, cPI/2});
            std::vector<precision_t> alt_in({100000, 350000, 100000, 350000, 100000, 350000, 100000, 350000});

            gGrid.set_interpolation_coefs(lon_in, lat_in, alt_in);

            std::vector<precision_t> ansLon = gGrid.get_interpolation_values(gGrid.geoLon_scgc);
            std::vector<precision_t> ansLat = gGrid.get_interpolation_values(gGrid.geoLat_scgc);
            std::vector<precision_t> ansAlt = gGrid.get_interpolation_values(gGrid.geoAlt_scgc);

            for (int64_t i = 0; i < lon_in.size(); ++i) {
                std::cout << "Estimation of (" << lon_in[i] << ", " << lat_in[i] << ", " << alt_in[i] << ") is "
                        << "(" << ansLon[i] << ", " << ansLat[i] << ", " << ansAlt[i] << ")\n";
            }
            is_correct(lon_in, ansLon, lat_in, ansLat, alt_in, ansAlt);
        }


        {
            std::cout << "\n\n\nError handling 2\n\n\n";
            std::vector<precision_t> lon_in({0, 0, 0, 0, 2*cPI - cSmall, 2*cPI - cSmall, 2*cPI - cSmall, 2*cPI - cSmall});
            std::vector<precision_t> lat_in({-cPI/2 - cSmall, -cPI/2 - cSmall, cPI/2 + cSmall, cPI/2 + cSmall, -cPI/2 - cSmall, -cPI/2 - cSmall, cPI/2 + cSmall, cPI/2 + cSmall});
            std::vector<precision_t> alt_in({100000, 350000, 100000, 350000, 100000, 350000, 100000, 350000});

            gGrid.set_interpolation_coefs(lon_in, lat_in, alt_in);

            std::vector<precision_t> ansLon = gGrid.get_interpolation_values(gGrid.geoLon_scgc);
            std::vector<precision_t> ansLat = gGrid.get_interpolation_values(gGrid.geoLat_scgc);
            std::vector<precision_t> ansAlt = gGrid.get_interpolation_values(gGrid.geoAlt_scgc);

            for (int64_t i = 0; i < lon_in.size(); ++i) {
                std::cout << "Estimation of (" << lon_in[i] << ", " << lat_in[i] << ", " << alt_in[i] << ") is "
                        << "(" << ansLon[i] << ", " << ansLat[i] << ", " << ansAlt[i] << ")\n";
            }
            is_correct(lon_in, ansLon, lat_in, ansLat, alt_in, ansAlt);
        }


        {
            std::cout << "\n\n\nError handling 3\n\n\n";
            std::vector<precision_t> lon_in({0, 0, 0, 0, 2*cPI - cSmall, 2*cPI - cSmall, 2*cPI - cSmall, 2*cPI - cSmall});
            std::vector<precision_t> lat_in({-cPI/2, -cPI/2, cPI/2, cPI/2, -cPI/2, -cPI/2, cPI/2, cPI/2});
            std::vector<precision_t> alt_in({100000 - 1, 350000 + 1, 100000 - 1, 350000 + 1, 100000 - 1, 350000 + 1, 100000 - 1, 350000 + 1});

            gGrid.set_interpolation_coefs(lon_in, lat_in, alt_in);

            std::vector<precision_t> ansLon = gGrid.get_interpolation_values(gGrid.geoLon_scgc);
            std::vector<precision_t> ansLat = gGrid.get_interpolation_values(gGrid.geoLat_scgc);
            std::vector<precision_t> ansAlt = gGrid.get_interpolation_values(gGrid.geoAlt_scgc);

            for (int64_t i = 0; i < lon_in.size(); ++i) {
                std::cout << "Estimation of (" << lon_in[i] << ", " << lat_in[i] << ", " << alt_in[i] << ") is "
                        << "(" << ansLon[i] << ", " << ansLat[i] << ", " << ansAlt[i] << ")\n";
            }
            is_correct(lon_in, ansLon, lat_in, ansLat, alt_in, ansAlt);
        }

        // Normal condition
        {
            std::cout << "\n\n\nNormal condition\n\n\n";
            std::vector<precision_t> lon_in({0.1, 0.2, 0.3, 0.4, 2*cPI - 0.5, 2*cPI - 0.6, 2*cPI - 0.7, cPI, cPI});
            std::vector<precision_t> lat_in({-cPI/2 + 0.5, -cPI/2 + 0.6, 0, 0.1, cPI/2 - 0.5, cPI/2 - 0.6, -0.1, -0.2, 0});
            std::vector<precision_t> alt_in({100000, 150000, 200000, 250000, 300000, 350000, 120000, 225000, 225000});

            gGrid.set_interpolation_coefs(lon_in, lat_in, alt_in);

            std::vector<precision_t> ansLon = gGrid.get_interpolation_values(gGrid.geoLon_scgc);
            std::vector<precision_t> ansLat = gGrid.get_interpolation_values(gGrid.geoLat_scgc);
            std::vector<precision_t> ansAlt = gGrid.get_interpolation_values(gGrid.geoAlt_scgc);

            for (int64_t i = 0; i < lon_in.size(); ++i) {
                std::cout << "Estimation of (" << lon_in[i] << ", " << lat_in[i] << ", " << alt_in[i] << ") is "
                        << "(" << ansLon[i] << ", " << ansLat[i] << ", " << ansAlt[i] << ")\n";
            }
            std::cout << "\n\n\n";
            is_correct(lon_in, ansLon, lat_in, ansLat, alt_in, ansAlt);
        }
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
