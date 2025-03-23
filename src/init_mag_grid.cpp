// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>

#include "aether.h"

// ----------------------------------------------------------------------
// Routine to convert p and q to r and theta. Can be solved iteratively,
// or with approach from (Swisdak, 2006), who solved it analytically:
//  https://arxiv.org/pdf/physics/0606044
//
// ----------------------------------------------------------------------

std::pair<precision_t, precision_t> qp_to_r_theta(precision_t q,
                                                  precision_t p) {

  // return quanties
  precision_t r, theta;
  // Intermediate quantities:
  precision_t term0, term1, term2, term3;

  term0 = 256.0 / 27.0 * pow(q, 2.0) * pow(p, 4.0);
  term1 = pow((1.0 + sqrt(1.0 + term0)), 2.0 / 3.0);
  term2 = pow(term0, 1.0 / 3.0);
  term3 = 0.5 * pow(((pow(term1, 2) + term1 * term2 + pow(term2, 2)) / term1),
                    3.0 / 2.0);

  r = p * (4.0 * term3) / ((1.0 + term3) * (1.0 + sqrt(2.0 * term3 - 1.0)));

  // now that r is determined we can solve for theta
  // theta = asin(sqrt(r/p));
  theta = acos(q * pow(r, 2.0));
  // Then make sure its the correct sign & direction (not colatitude)
  theta = cPI / 2 - theta;

  return {r, theta};
}

////////////////////////////////////////////
// convert cell coordinates to geographic //
////////////////////////////////////////////
std::vector <arma_cube> mag_to_geo(arma_cube magLon, arma_cube magLat,
                                   arma_cube magAlt,
                                   Planets planet) {
  std::string function = "Grid::mag_to_geo";
  static int iFunction = -1;
  report.enter(function, iFunction);

  std::vector<arma_cube> llr, xyz_mag, xyz_geo, xyzRot1, xyzRot2;
  llr.push_back(magLon);
  llr.push_back(magLat);
  llr.push_back(magAlt);
  xyz_mag = transform_llr_to_xyz_3d(llr);

  precision_t magnetic_pole_rotation = planet.get_dipole_rotation();
  precision_t magnetic_pole_tilt = planet.get_dipole_tilt();
  std::vector<precision_t> dipole_center = planet.get_dipole_center();

  // Reverse our dipole rotations:
  xyzRot1 = rotate_around_y_3d(xyz_mag, magnetic_pole_tilt);
  xyzRot2 = rotate_around_z_3d(xyzRot1, magnetic_pole_rotation);

  // offset dipole (not fully suported yet, so will be zero. 
  if ((dipole_center[0] != 0) || (dipole_center[1] != 0) ||
    (dipole_center[2] != 0)) {
  report.error("Dipole center != 0, but that is not supported yet. Setting to 0!");
  dipole_center = {0, 0, 0};
  }

  xyz_geo[0] = xyzRot2[0] + dipole_center[0];
  xyz_geo[1] = xyzRot2[1] + dipole_center[1];
  xyz_geo[2] = xyzRot2[2] + dipole_center[2];

  // transform back to lon, lat, radius:
  llr = transform_xyz_to_llr_3d(xyzRot2);

  report.exit(function);
  return llr;
}

// ----------------------------------------------------------------------
// Initialize the dipole grid.
// - inputs (min_apex, min_alt, LatStretch, FieldLineStretch, max_lat_dipole)
//   are read from input files. And the numbers of each coordinate.
// - nLats must be even!!
// ----------------------------------------------------------------------
bool Grid::init_dipole_grid(Quadtree quadtree_ion, Planets planet) {

  using namespace std;
  bool DidWork = true;

  string function = "Grid::init_dipole_grid";
  static int iFunction = -1;
  report.enter(function, iFunction);

  // turn the switch on!
  IsGeoGrid = false;
  IsMagGrid = true;
  IsCubeSphereGrid = false;
  IsDipole = true;

  // report.print(0, "Creating inter-node connections Grid");

  //if (!Is0D & !Is1Dz)
  //  create_sphere_connection(quadtree_ion);

  report.print(0, "Creating Dipole Grid");

  report.print(3, "Getting grid inputs for dipole grid");

  Inputs::grid_input_struct grid_input = input.get_grid_inputs("ionGrid");

  // Number of ghost cells:
  int64_t nGCs = get_nGCs();

  // Get inputs:

  precision_t min_lat = grid_input.lat_min;
  precision_t max_lat = grid_input.lat_max;

  precision_t min_alt = grid_input.alt_min * cKMtoM;
  precision_t max_alt = grid_input.alt_max * cKMtoM;

  // Normalize inputs to planet radius... (update when earth is oblate)
  precision_t planetRadius = planet.get_radius(0.0);
  // Altitude to begin modeling, normalized to planet radius
  precision_t min_alt_re = (min_alt + planetRadius) / planetRadius;
  precision_t max_alt_re = (max_alt + planetRadius) / planetRadius;

  if (nAlts % 2 != 0) {
    report.error("nAlts must be even!");
    DidWork = false;
  }

  // Get some coordinates and sizes in normalized coordinates:
  arma_vec lower_left_norm = quadtree_ion.get_vect("LL"); // origin
  arma_vec size_right_norm = quadtree_ion.get_vect("SR"); // lon_lims
  arma_vec size_up_norm = quadtree_ion.get_vect("SU");    // lat_extent
  report.print(3, "Got all settings. Initializing longitudes.");

  precision_t dlon = size_right_norm(0) * cPI / (nLons - 2 * nGCs);
  precision_t lon0 = lower_left_norm(0) * cPI;
  arma_vec lon1d(nLons);

  arma_vec lon1dLeft(nLons + 1);

  // if we are not doing anything in the lon direction, then set dlon to
  // something reasonable:
  if (!HasXdim)
    dlon = 1.0 * cDtoR;

  // Dimension iterators
  int64_t iLon, iLat, iAlt;

  /////////////////
  // Longitudes: //
  /////////////////

  // - Make a 1d vector
  // - copy it into the 3d cube
  for (iLon = 0; iLon < nLons; iLon++) {
    lon1d(iLon) = lon0 + (iLon - nGCs + 0.5) * dlon;
    lon1dLeft(iLon) = lon0 + (iLon - nGCs) * dlon; // corners
  }

  lon1dLeft(nLons) = lon0 + (nLons - nGCs) * dlon;

  for (iLat = 0; iLat < nLats; iLat++) {
    for (iAlt = 0; iAlt < nAlts; iAlt++) {
      // centers:
      magLon_scgc.subcube(0, iLat, iAlt, nLons - 1, iLat, iAlt) = lon1d;
      i_center_scgc.subcube(0, iLat, iAlt, nLons - 1, iLat, iAlt) = lon1d;
      // left edges
      magLon_Left.subcube(0, iLat, iAlt, nLons, iLat, iAlt) = lon1dLeft;
      i_edge_scgc.subcube(0, iLat, iAlt, nLons, iLat, iAlt) = lon1dLeft;
      // corners
      magLon_Corner.subcube(0, iLat, iAlt, nLons, iLat, iAlt) = lon1dLeft;
      i_corner_scgc.subcube(0, iLat, iAlt, nLons, iLat, iAlt) = lon1dLeft;
    }
  }

  if (magLon_scgc.has_nan())
    report.error("NAN IN MAGLON");

  report.print(3, "Done initializing longitudes, moving to latitude");

  ////////////////
  // Latitudes: //
  ////////////////

  // Invariant latitude is evenly spaced across each block.
  // Latitude limits are adjusted here, not in quadtree

  // - From the quadtree, we see the origin & extent of this block
  // - That is normalized, without any influence from settings
  // - Scale it with the latitude limits provided by the user
  // - Put invariant latitudes down, linearly, between this range.

  // This has to be done differently in the north & south hemisphere.
  // So note if we are in the southern hemisphere and invert it afterwards.

  bool isSouth = false;
  precision_t lat_origin = lower_left_norm(1);

  if (lat_origin < -0.01) { // handles some imprecision
    isSouth = true;
    lat_origin = -1.0 * lat_origin - size_up_norm(1);
  }

  precision_t lat0 = 2.0 * (max_lat - min_lat) * lat_origin;
  precision_t dlat = 2.0 * size_up_norm(1) * (max_lat -  min_lat) /
                     (nLats - nGCs);

  arma_vec lat1d(nLats);
  arma_vec lat1dDown(nLats + 1);

  for (iLat = 0; iLat < nLats; iLat++) {
    lat1d(iLat) = lat0 + (iLat - nGCs + 0.5) * dlat + min_lat; // centers
    lat1dDown(iLat) = lat0 + (iLat - nGCs) * dlat + min_lat; // corners & edges
  }

  lat1dDown(nLats) = lat0 + (nLats - nGCs) * dlat; // last corner

  // At the pole:
  // - put last ghost cell's corner at 89.9 degrees latitude
  // - put 2nd to last corner 1/2 way between 89.9 and the last real corner
  // - evenly space the ghost cells between these.

  // Check if we're touching the pole, need to look at original quadtree values
  if ((abs(lower_left_norm(1) + size_up_norm(1)) > 0.49) // north pole
      || (lower_left_norm(1) < -0.49)) { // south pole
    lat1dDown(nLats) = 89.9 * cDtoR;
    lat1dDown(nLats - 1) = (lat1dDown(nLats) + lat1dDown(nLats - 2)) / 2.0;
    lat1d(nLats - 1) = (lat1dDown(nLats) + lat1dDown(nLats - 1)) / 2.0;
    lat1d(nLats - 2) = (lat1dDown(nLats - 1) + lat1dDown(nLats - 2)) / 2.0;
  }

  // l-shells of centers
  arma_vec Pcenters = min_alt_re / pow(sin(cPI / 2 - lat1d), 2);

  // l-shells of corners
  arma_vec Pcorners = min_alt_re / pow(sin(cPI / 2 - lat1dDown), 2);

  report.print(3, "Done initializing invariant latitudes");

  ////////////////
  // Altitudes: //
  ////////////////

  // - Trace each field line from q_min to q_max. Identical for all field lines within this block.
  // - Obtain the minimum "altitude" (q) from the highest latitude
  // field line on each block
  // - Obtain the maximum "altitude" from the lowest latitude *open* field line.
  // (closed field lines are treated differently)
  //  - In other words, since we are tracing from q_min to q_max, use the highest field
  //    line to get q_min and the lowest for q_max. This forces all field lines to
  //    start & end within the bounds.
  // - Evenly space all points' "altitude" linear across these two values
  // - Altitude here refers to the dipole q-coordinate - cos(magLat)/r^2
  // - Blocks touching a pole or the equator are treated differently

  // Field lines close if:
  // - touching the (magnetic) equator
  // - minimum Lshell in this block is < max_alt (the q-value would be undefined)

  precision_t q_min;
  bool close_this_block = false;

  if (Pcorners.min() < max_alt_re) // invalid q's - Lshell < max_alt
    close_this_block = true;

  if (lat_origin < 0.01) // equator, with some imprecision
    close_this_block = true;

  if (close_this_block)
    q_min = 0; // q=0 at equator (for closed blocks)
  else
    // invLats are still all in North Hemisphere & increasing.
    // Use minimum p & alt to solve for q
    // q = sqrt((1-r/p)/r^4)
    q_min = pow(((1 - max_alt_re / Pcenters(nGCs)) / pow(max_alt_re, 4.0)), 0.5);

  // Trace each field line up to q_max, obtained from the lowest field line in the block
  precision_t q_max = pow(((1 - min_alt_re / Pcenters(nLats - nGCs)) / pow(
                             min_alt_re,
                             4.0)), 0.5);

  // Counter-intuitive, but the maximum value of q is actually where we start
  // (lowest altitude), since q=0 at equator.
  precision_t delQ = (q_max - q_min) / (nAlts - nGCs * 2.0);

  arma_vec magQ1d(nAlts);
  arma_vec magQ_corner_1d(nAlts + 1);

  for (iAlt = 0; iAlt < nAlts; iAlt ++) {
    magQ1d(iAlt) = q_min + (iAlt - nGCs + 0.5) * delQ;
    magQ_corner_1d(iAlt) = q_min + (iAlt - nGCs) * delQ;
  }

  magQ_corner_1d(nAlts) = q_min - nGCs * delQ;

  report.print(3,
               "Done generating points for magnetic grid. Plugging everything in");

  ////////////////////////////
  // That is the grid made. //
  ////////////////////////////
  // Now to store everything....
  // It's all done at the end to make things more simple earlier, but that makes this part messier.

  // temp holding names:
  std::pair<precision_t, precision_t> rtheta, rtheta_edge;
  precision_t radius, radius_edge, theta, theta_edge, invLat, invLat_edge,
              pcenter, pedge, qcenter, qedge;
  // we need to turn single floats into vectors/cubes:

  // We can solve for (r, theta) for each point on the (q,p) grid. Do that & store:
  // Currently the grid is symmetric in longitude.
  // Interte through centers & edges first, then do corners afterwards
  for (iLat = 0; iLat < nLats; iLat ++) {
    for (iAlt = 0; iAlt < nAlts; iAlt++) {
      // We have to reverse & negate things; want latitudes from south->north
      // and altitude low->high. Altitude is in the correct direction, so change how we
      // access values in the latitude dimension.

      if (isSouth) {
        qcenter = magQ1d(iAlt);
        pcenter = Pcenters(nLats - iLat - 1);

        qedge = magQ_corner_1d(iAlt);
        pedge = Pcorners(nLats - iLat - 1);

        invLat = lat1d(nLats - iLat - 1) * -1;
        invLat_edge = lat1dDown(nLats - iLat - 1) * -1;

        rtheta = qp_to_r_theta(qcenter, pcenter);
        rtheta_edge = qp_to_r_theta(qcenter, pcenter);

        // Flip hemisphere of latitude & q (cannot be done before qp_to_rtheta)
        radius = rtheta.first;
        theta = rtheta.second * -1.0;
        qcenter *= -1.0;

        radius_edge = rtheta_edge.first;
        theta_edge = rtheta.second * -1.0;
        qedge *= 1.0;
      } else {
        qcenter = magQ1d(iAlt);
        pcenter = Pcenters(iLat);

        qedge = magQ_corner_1d(iAlt);
        pedge = Pcorners(iLat);

        invLat = lat1d(iLat);
        invLat_edge = lat1dDown(iLat);

        rtheta = qp_to_r_theta(qcenter, pcenter);
        rtheta_edge = qp_to_r_theta(qedge, pedge);

        radius = rtheta.first;
        theta = rtheta.second;

        radius_edge = rtheta_edge.first;
        theta_edge = rtheta_edge.second;
      }

      for (iLon = 0; iLon < nLons; iLon ++) {
        magLat_scgc(iLon, iLat, iAlt) = theta;
        j_center_scgc(iLon, iLat, iAlt) = theta;

        magLat_Down(iLon, iLat, iAlt) = theta_edge;
        j_edge_scgc(iLon, iLat, iAlt) = theta_edge;

        magAlt_scgc(iLon, iLat, iAlt) = radius;
        k_center_scgc(iLon, iLat, iAlt) = radius;

        magAlt_Below(iLon, iLat, iAlt) = radius_edge;
        k_edge_scgc(iLon, iLat, iAlt) = radius_edge;

        // extra coordinates
        magP_scgc(iLon, iLat, iAlt) = pcenter;
        magQ_scgc(iLon, iLat, iAlt) = qcenter;
        magInvLat_scgc(iLon, iLat, iAlt) = invLat;
      }
    }
  }

  report.print(3, "Centers are in");

  precision_t radius_corner, theta_corner, invLat_corner,
              pcorner, qcorner;

  for (iLat = 0; iLat < nLats + 1; iLat ++) {
    for (iAlt = 0; iAlt < nAlts + 1; iAlt++) {

      // Same process as the centers & edges (above)
      if (isSouth) {
        qcorner = magQ_corner_1d(iAlt);
        pcorner = Pcorners(nLats - iLat);
        invLat_corner = lat1dDown(nLats - iLat) * -1;
        rtheta = qp_to_r_theta(qcorner, pcorner);

        radius_corner = rtheta.first;
        theta_corner = rtheta.second * -1;
        qcorner *= -1;
      } else {
        qcorner = magQ_corner_1d(iAlt);
        pcorner = Pcorners(iLat);
        invLat_corner = lat1dDown(iLat);
        rtheta = qp_to_r_theta(qcorner, pcorner);
        radius_corner = rtheta.first;
        theta_corner = rtheta.second;
      }

      for (iLon = 0; iLon < nLons + 1; iLon ++) {
        magLat_Corner(iLon, iLat, iAlt) = theta_corner;
        j_corner_scgc(iLon, iLat, iAlt) = theta_corner;

        magAlt_Corner(iLon, iLat, iAlt) = radius_corner;
        k_corner_scgc(iLon, iLat, iAlt) = radius_corner;

        magP_Corner(iLon, iLat, iAlt) = pcorner;
        magQ_Corner(iLon, iLat, iAlt) = qcorner;
        magInvLat_Corner(iLon, iLat, iAlt) = invLat_corner;
      }
    }
  }

  report.print(3, "Corners done too");

  // all distances, so far, are in units of planet radii, turn into meters.
  // Except for Q, leave that dimensionless.
  magAlt_scgc *= planetRadius;
  k_center_scgc *= planetRadius;
  magAlt_Below *= planetRadius;
  k_edge_scgc *= planetRadius;
  magP_scgc *= planetRadius;
  magAlt_Corner *= planetRadius;
  k_corner_scgc *= planetRadius;
  magP_Corner *= planetRadius;


  std::vector <arma_cube> llr = mag_to_geo(magLon_scgc, magLat_scgc, magAlt_scgc,
                                           planet);

  geoLon_scgc = llr[0];
  geoLat_scgc = llr[1];
  geoAlt_scgc = llr[2] - planet.get_radius(geoLat_scgc);
  report.print(4,
               "Done dipole -> geographic transformations for the dipole grid centers.");

  std::vector <arma_cube> llr_corner = mag_to_geo(magLon_Corner, magLat_Corner,
                                                  magAlt_Corner, planet);
  geoLon_Corner = llr_corner[0];
  geoLat_Corner = llr_corner[1];
  geoAlt_Corner = llr_corner[2] - planetRadius;
  report.print(4,
               "Done dipole -> geographic transformations for the dipole grid centers.");

  // Calculate the radius, of planet
  fill_grid_radius(planet);

  // Figure out what direction is radial:
  rad_unit_vcgc = make_cube_vector(nLons, nLats, nAlts, 3);
  gravity_vcgc = make_cube_vector(nLons, nLats, nAlts, 3);

  for (int iV = 0; iV < 3; iV++) {
    rad_unit_vcgc[iV].zeros();
    gravity_vcgc[iV].zeros();
  }

  arma_cube br = 2 * sin(abs(magLat_scgc));
  arma_cube bt = cos(magLat_scgc);
  arma_cube bm = sqrt(br % br + bt % bt);
  // Latitudinal direction of radial:
  arma_cube s = sign(magLat_scgc);
  s.elem(find(s == 0)).ones();

  rad_unit_vcgc[1] = bt / bm % s;
  rad_unit_vcgc[2] = -br / bm;

  precision_t mu = planet.get_mu();
  gravity_vcgc[1] = mu * rad_unit_vcgc[1] % radius2i_scgc;
  gravity_vcgc[2] = mu * rad_unit_vcgc[2] % radius2i_scgc;
  gravity_potential_scgc.set_size(nX, nY, nAlts);
  gravity_potential_scgc.zeros();
  gravity_mag_scgc = sqrt(
                       gravity_vcgc[0] % gravity_vcgc[0] +
                       gravity_vcgc[1] % gravity_vcgc[1] +
                       gravity_vcgc[2] % gravity_vcgc[2]);

  report.print(4, "Done gravity calculations for the dipole grid.");

  calc_dipole_grid_spacing(planet);

  report.print(4, "Done altitude spacing for the dipole grid.");

  // Calculate magnetic field and magnetic coordinates:
  fill_grid_bfield(planet);
  report.print(4, "Done filling dipole grid with b-field!");


  // put back into altitude. we've been carrying around radius:
  // magAlt_scgc = magAlt_scgc - planetRadius;
  // this breaks things more???

  report.exit(function);
  return DidWork;
}
