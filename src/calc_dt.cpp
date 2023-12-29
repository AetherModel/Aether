// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "../include/aether.h"

// --------------------------------------------------------------------------
//
// --------------------------------------------------------------------------


precision_t calc_dt(Grid grid, std::vector<arma_cube> cMax_vcgc) {

    std::string function = "calc_dt";
    static int iFunction = -1;
    report.enter(function, iFunction);

    precision_t dt;

    if (input.get_is_cubesphere())
        dt = calc_dt_cubesphere(grid, cMax_vcgc);
    else
        dt = calc_dt_sphere(grid, cMax_vcgc);

    report.exit(function);
    return dt;
}

// --------------------------------------------------------------------------
//
// --------------------------------------------------------------------------

precision_t calc_dt_sphere(Grid grid, std::vector<arma_cube> cMax_vcgc) {

    std::string function = "calc_dt_sphere";
    static int iFunction = -1;
    report.enter(function, iFunction);

    precision_t dt;

    arma_vec dta(4);
    arma_cube dtCube;

    // Longitudinal Direction:
    dtCube = grid.dlon_center_dist_scgc / cMax_vcgc[0];
    dta(0) = dtCube.min();

    // Latitudinal Direction:
    dtCube = grid.dlat_center_dist_scgc / cMax_vcgc[1];
    dta(1) = dtCube.min();

    // Vertical Direction:
    dta(2) = calc_dt_vertical(grid, cMax_vcgc);

    // Set a minimum dt:
    dta(3) = 10.0;

    dt = dta.min();

    if (report.test_verbose(3))
        std::cout << "dt (sphere) for neutrals : " << dt << "\n";

    if (report.test_verbose(4))
        std::cout << " derived from dt(x, y, z, extra) : " << dta << "\n";

    report.exit(function);
    return dt;
}

// --------------------------------------------------------------------------
//
// --------------------------------------------------------------------------

precision_t calc_dt_cubesphere(Grid grid, std::vector<arma_cube> cMax_vcgc) {

    std::string function = "calc_dt_sphere";
    static int iFunction = -1;
    report.enter(function, iFunction);

    precision_t dt;
    arma_vec dta(4);

    // Get some dimensions
    int64_t nAlts = grid.get_nAlts();
    int64_t nXs = grid.get_nLons();
    int64_t nYs = grid.get_nLats();

    // dtx dty for reference coordinate system
    arma_cube dtx(nXs, nYs, nAlts);
    arma_cube dty(nXs, nYs, nAlts);

    // A dummy constant one matrix
    arma_mat dummy_1(nXs, nYs, fill::ones);

    // Loop through altitudes
    for (int iAlt = 0; iAlt < nAlts; iAlt++) {
        // Conver cMax to contravariant velocity first
        arma_mat u1 = sqrt(
                          cMax_vcgc[0].slice(iAlt) % grid.A11_inv_scgc.slice(iAlt) %
                          cMax_vcgc[0].slice(iAlt) % grid.A11_inv_scgc.slice(iAlt) +
                          cMax_vcgc[1].slice(iAlt) % grid.A12_inv_scgc.slice(iAlt) %
                          cMax_vcgc[1].slice(iAlt) % grid.A12_inv_scgc.slice(iAlt));
        arma_mat u2 = sqrt(
                          cMax_vcgc[0].slice(iAlt) % grid.A21_inv_scgc.slice(iAlt) %
                          cMax_vcgc[0].slice(iAlt) % grid.A21_inv_scgc.slice(iAlt) +
                          cMax_vcgc[1].slice(iAlt) % grid.A22_inv_scgc.slice(iAlt) %
                          cMax_vcgc[1].slice(iAlt) % grid.A22_inv_scgc.slice(iAlt));
        dtx.slice(iAlt) = grid.drefx(iAlt) * dummy_1 / u1;
        dty.slice(iAlt) = grid.drefy(iAlt) * dummy_1 / u2;
    }
    // Take minimum dts in each direction:
    dta(0) = dtx.min();
    dta(1) = dty.min();
    // Vertical Direction:
    dta(2) = calc_dt_vertical(grid, cMax_vcgc);
    // Set a minimum dt:
    dta(3) = 10.0;
    // Take the minimum of all directions:
    dt = dta.min();

    if (report.test_verbose(3))
        std::cout << "dt (cubesphere) : " << dt << "\n";

    if (report.test_verbose(4))
        std::cout << " derived from dt(x, y, z, extra) : " << dta << "\n";

    report.exit(function);
    return dt;
}

// --------------------------------------------------------------------------
//
// --------------------------------------------------------------------------

precision_t calc_dt_vertical(Grid grid, std::vector<arma_cube> cMax_vcgc) {

    std::string function = "calc_dt_vertical";
    static int iFunction = -1;
    report.enter(function, iFunction);

    precision_t dt;
    if (input.get_nAltsGeo() > 1) {
        arma_cube dtz = grid.dalt_center_scgc / cMax_vcgc[2];
        dt = dtz.min();
    } else
        dt = 1e32;

    report.exit(function);
    return dt;
}