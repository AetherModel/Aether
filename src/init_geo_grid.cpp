// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>

#include "aether.h"
#include <math.h>

// ----------------------------------------------------------------------
// Create a spherical grid with lon/lat/alt coordinates
// ----------------------------------------------------------------------

void Grid::create_altitudes(Planets planet) {

  std::string function = "Grid::create_altitudes";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int64_t iLon, iLat, iAlt;

  arma_vec alt1d(nAlts);

  Inputs::grid_input_struct grid_input;

  grid_input = input.get_grid_inputs(gridType);

  if (grid_input.IsUniformAlt) {
    for (iAlt = 0; iAlt < nAlts; iAlt++)
      // Convert km to m:
      alt1d(iAlt) = (grid_input.alt_min + (iAlt - nGeoGhosts) * grid_input.daltKm) *
                    cKMtoM;
  } else {

    json neutrals = planet.get_neutrals();
    json temperatures = planet.get_temperatures();
    std::vector<double> input_alt;
    std::vector<double> input_temp;

    for (int i = 0; i < temperatures["alt"].size(); i++) {
      input_alt.push_back(double(temperatures["alt"][i]) * 1000.0);
      input_temp.push_back(temperatures["temp"][i]);
    }

    precision_t scale_height, temperature, gravity, radius, mass, density;
    int64_t nSp = neutrals["name"].size();
    arma_vec densities(nSp);
    arma_vec masses(nSp);
    arma_vec h(nSp);

    int64_t iSp;

    report.print(1, "Making non-uniform altitude grid!");

    if (grid_input.daltScale > 0.5) {
      if (report.test_verbose(0)) {
        std::cout << "-----------------------------------------------------\n";
        std::cout << "WARNING: daltScale is set to > 0.5, with non-uniform grid!\n";
        std::cout << "   daltScale = " << grid_input.daltScale << "\n";
        std::cout << "-----------------------------------------------------\n";
      }
    }

    // Convert to km
    double alt = grid_input.alt_min * cKMtoM;
    radius = planet.get_radius(0.0) + alt;
    precision_t mu = planet.get_mu();
    gravity = mu / (radius * radius);

    temperature = interpolate_1d(alt, input_alt, input_temp);

    mass = 0.0;
    density = 0.0;

    for (iSp = 0; iSp < nSp; iSp++) {
      masses(iSp) = double(neutrals["mass"][iSp]) * cAMU;
      densities[iSp] = neutrals["BC"][iSp];
      h(iSp) = cKB * temperature / (masses(iSp) * gravity);
      mass = mass + masses(iSp) * densities[iSp];
      density = density + densities[iSp];
    }

    // convert mass density into mass:
    mass = mass / density;
    scale_height = cKB * temperature / (mass * gravity);

    precision_t dalt = scale_height * grid_input.daltScale;
    precision_t dAltLimiter = dalt * 10.0;

    // Fills bottom ghost cells with constant dAlt
    // Fills bottom cell with actual desired bottom altitude
    for (iAlt = 0; iAlt <= nGeoGhosts; iAlt++) {
      alt1d(iAlt) = grid_input.alt_min * cKMtoM + (iAlt - nGeoGhosts) * dalt;

      if (report.test_verbose(1))
        std::cout << "iAlt : " << iAlt
                  << " Altitude : " << alt1d(iAlt) / 1000.0
                  << " (km)\n";
    }

    for (iAlt = nGeoGhosts + 1; iAlt < nAlts; iAlt++) {

      alt = alt1d(iAlt - 1);
      temperature = interpolate_1d(alt, input_alt, input_temp);
      radius = planet.get_radius(0.0) + alt;
      gravity = mu / (radius * radius);

      mass = 0.0;
      density = 0.0;

      for (iSp = 0; iSp < nSp; iSp++) {
        mass = mass + masses(iSp) * densities[iSp];
        density = density + densities[iSp];
      }

      // convert mass density into mass:
      mass = mass / density;
      scale_height = cKB * temperature / (mass * gravity);

      dalt = scale_height * grid_input.daltScale;

      if (dalt > dAltLimiter)
        dalt = dAltLimiter;

      alt1d(iAlt) = alt + dalt;

      h = cKB * temperature / (masses * gravity);
      densities = densities % exp(-dalt / h);

      if (report.test_verbose(1))
        std::cout << "iAlt : " << iAlt
                  << " Altitude : " << alt1d(iAlt) / 1000.0
                  << " (km)\n";
    }
  }

  // This takes cell centers and calculates the edges:
  arma_vec alt1d_below = calc_bin_edges(alt1d);

  for (iLon = 0; iLon < nLons; iLon++) {
    for (iLat = 0; iLat < nLats; iLat++) {
      geoAlt_scgc.tube(iLon, iLat) = alt1d;
      k_center_scgc.tube(iLon, iLat) = alt1d;
      geoAlt_Below.tube(iLon, iLat) = alt1d_below;
      k_edge_scgc.tube(iLon, iLat) = alt1d_below;
    }
  }

  for (iLon = 0; iLon < nLons + 1; iLon++) {
    for (iLat = 0; iLat < nLats + 1; iLat++) {
      geoAlt_Corner.tube(iLon, iLat) = alt1d_below;
      k_corner_scgc.tube(iLon, iLat) = alt1d_below;
    }
  }

  // All cells on the geographic grid *should* be ok
  isTooLowCell = find(geoAlt_scgc < grid_input.alt_min * cKMtoM);
  isPhysicalCell = find(geoAlt_scgc >= grid_input.alt_min * cKMtoM);
  // get the ghost cell indices on each lat/lon point.
  // may be redundant can fill lower with nGCs-1, but this is here for now
  arma::uvec theGCs;

  for (iLon = 0; iLon < nLons; iLon++) {
    for (iLat = 0; iLat < nLats; iLat++) {
      // find *last* cell below alt_min
      theGCs = find(geoAlt_scgc.tube(iLon, iLat) < grid_input.alt_min * cKMtoM);
      // Get the last element if the col-vec
      first_lower_gc(iLon, iLat) = theGCs(theGCs.n_elem - 1);
    }
  }

  first_upper_gc.fill(nAlts - nGCs * 2 - 1);

  report.exit(function);
  return;
}

// ----------------------------------------------------------------------
// Initialize the geographic grid.  At the moment, this is a simple
// Lon/Lat/Alt grid.  The grid structure is general enough that each
// of the lon, lat, and alt can be a function of the other variables.
// ----------------------------------------------------------------------

bool Grid::init_geo_grid(Quadtree quadtree,
                         Planets planet) {

  std::string function = "Grid::init_geo_grid";
  static int iFunction = -1;
  report.enter(function, iFunction);
  bool DidWork = true;

  IsGeoGrid = true;

  if (iGridShape_ == iCubesphere_) {
    report.print(0, "Creating Cubesphere Grid for : " + gridType);

    if (!Is0D & !Is1Dz)
      create_cubesphere_connection(quadtree);

    IsCubeSphereGrid = true;
  } else {
    report.print(0, "Creating Spherical Grid for : " + gridType);

    if (!Is0D & !Is1Dz)
      create_sphere_connection(quadtree);

    IsCubeSphereGrid = false;
  }

  //if (input.get_do_restart() & iGridShape_ != iCubesphere_) {
  //  report.print(1, "Restarting! Reading grid files!");
  //  DidWork = read_restart(input.get_restartin_dir());
  //} else {
  if (iGridShape_ == iCubesphere_)
    create_cubesphere_grid(quadtree);

  else
    create_sphere_grid(quadtree);

  //MPI_Barrier(aether_comm);
  create_altitudes(planet);

  // set the altitude of the lower boundary values:
  altitude_lower_bc = planet.get_altitude_of_bc();

  init_connection();

  //DidWork = write_restart(input.get_restartout_dir());
  //}

  // Calculate the radius (for spherical or non-spherical)
  fill_grid_radius(planet);

  // Correct the reference grid with correct length scale:
  // (with R = actual radius)
  if (iGridShape_ == iCubesphere_) {
    correct_xy_grid(planet);
    // New functions for equal-angular grid (center, left, down):
    report.print(2, "Scaling Cube by Radius");
    scale_cube_by_radius(cubeC);
    scale_cube_by_radius(cubeL);
    scale_cube_by_radius(cubeD);
    report.print(2, "Done Scaling Cube by Radius");
  }

  if (gridType == ionType_) {
    report.print(0, "--> Grid is Magnetic, so rotating");
    std::vector<arma_cube> llr, xyz, xyzRot1, xyzRot2;
    llr.push_back(geoLon_scgc);
    llr.push_back(geoLat_scgc);
    llr.push_back(radius_scgc);
    xyz = transform_llr_to_xyz_3d(llr);

    precision_t magnetic_pole_rotation = 265.0 * cDtoR;
    precision_t magnetic_pole_tilt = 10.0 * cDtoR;

    // Reverse our dipole rotations:
    xyzRot1 = rotate_around_y_3d(xyz, magnetic_pole_tilt);
    xyzRot2 = rotate_around_z_3d(xyzRot1, magnetic_pole_rotation);

    // transform back to lon, lat, radius:
    llr = transform_xyz_to_llr_3d(xyzRot2);

    geoLon_scgc = llr[0];
    geoLat_scgc = llr[1];
    geoAlt_scgc = llr[2] - planet.get_radius(0.0);

    IsGeoGrid = false;
  }

  // Calculate PFPC coordinates (i.e., XYZ from LLR)
  calc_xyz(planet);
  // Calculate grid spacing
  calc_grid_spacing(planet);
  //calculate radial unit vector (for spherical or oblate planet)
  calc_rad_unit(planet);
  // Calculate gravity (including J2 term, if desired)
  calc_gravity(planet);
  // Calculate magnetic field and magnetic coordinates:
  fill_grid_bfield(planet);

  write_restart(input.get_restartout_dir());

  // Throw a little message for students:
  report.student_checker_function_name(input.get_is_student(),
                                       input.get_student_name(),
                                       4, "");

  // The dipole grid has some variables that need to be set:
  IsClosed = false;
  setNorthAsDown = false;
  setSouthAsDown = false;

  report.exit(function);
  return DidWork;
}
