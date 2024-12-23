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
  // Then make sure its the correct sign & direction
  theta = cPI / 2 - theta;

  return {r, theta};
}

std::pair<arma_cube, arma_cube> qp_to_r_theta(arma_cube q, arma_cube p) {
  // return quanties
  arma_cube r, theta;
  // Intermediate quantities:
  arma_cube term0, term1, term2, term3;

  term0 = 256.0 / 27.0 * (q % q) % (p % p % p % p);
  term1 = pow((1.0 + sqrt(1.0 + term0)), 2.0 / 3.0);
  term2 = pow(term0, 1.0 / 3.0);
  term3 = 0.5 * pow(((term1 % term1 + term1 % term2 + term2 % term2) / term1),
                    3.0 / 2.0);

  r = p % (4.0 * term3) / ((1.0 + term3) % (1.0 + sqrt(2.0 * term3 - 1.0)));

  // now that r is determined we can solve for theta
  theta = asin(q % (r % r));

  return {r, theta};
}

arma_vec Grid::baselat_spacing(precision_t extent,
                               precision_t origin,
                               precision_t upper_lim,
                               precision_t lower_lim,
                               precision_t spacing_factor) {
  std::string function = "Grid::baselat_spacing";
  static int iFunction = -1;
  report.enter(function, iFunction);

  // intermediate latitude values
  precision_t lat_low, lat_high, lat_low0, lat_high0;
  // intermediate calculation values
  precision_t dlat, bb, aa, ang0, angq, nLats_here, extent_here;

  // Now we can allocate the return array,
  arma_vec Lats(nLats);

  // Noting the special case of 1 root node & 1 processor...
  bool DO_FLIPBACK = false;

  if (extent > 0.5) {
    DO_FLIPBACK = true;
    nLats_here = nLats / 2;
    extent_here = 0.5;
  } else
    nLats_here = nLats;

  // get the upper & lower latitude bounds for our division of the quadree
  if (origin < 0) {
    // negative origin: lat_high <=> lat_low
    lat_low, lat_high = -upper_lim, -lower_lim;
    lat_low0 = lat_low;
    lat_low = lat_low - (lat_high - lat_low) * (origin / 0.5);
    lat_high = lat_low0 - (lat_high - lat_low0) * (extent_here / .5 + origin / 0.5);
  } else {
    lat_low0 = lower_lim;
    lat_low = lat_low + (lat_high - lat_low) * (origin / 0.5);
    lat_high = lat_low0 + (lat_high - lat_low0) * (extent_here / .5 + origin / 0.5);
  }

  // normalized spacing in latitude
  // NOTE: spacing factor != 1 will not work yet. but framework is here...
  bb = (lat_high - lat_low) / (pow(lat_high, spacing_factor) - pow(lat_low,
                               spacing_factor));
  aa = lat_high - bb * pow(lat_high, spacing_factor);
  dlat = (lat_high - lat_low) / (nLats_here + 1);
  report.print(4, "baselats laydown!");

  // edge case for 1-D
  // In 1-D, the base latitudes will be 1/2 way between LatMax & minApex,
  // dlat is adjustable if it doesn't suit your needs.
  if (!HasYdim) {
    DO_FLIPBACK = false;
    dlat = 1.0 * cDtoR;
    nLats_here = nLats + 1;
  }

  for (int64_t j = 0; j < nLats_here; j++) {
    ang0 = lat_low + (j + 1) * dlat;
    angq = aa + bb * pow(ang0, spacing_factor);
    Lats[j] = angq;
  }

  report.print(5, "baselats flipback!");

  // In the flipback case (single processor, global sim), we want baselats
  // to be strictly increasing, same as geo grid!
  if (DO_FLIPBACK)
    for (int64_t j = 0; j < nLats_here; j++)
      Lats[j + nLats_here] = -1 * Lats[nLats_here - j - 1];

  report.print(4, "baselats flipback done!");

  report.exit(function);
  return Lats;
}

// // Gravity vectors in the dipole basis
// void calc_dipole_gravity(Planets planet){

// // rhat = -(2*cos/(del)) qhat + (sin/(del)) phat


// }


// === SPACING ALONG FIELD LINE === //
// Coordinates along the field line to begin modeling
// - Created in dipole (p,q) coordinates, stored as magnetic coords
// - North & south hemisphere base-latitudes, shouldn't be *too* hard to support offset
//   dipole and/or oblate Earth.
// isCorner is a bool, if false then the p's and q's are stored for later (p,q cell centers).
// Field line filling only needs to be redone for the "down" edges, left is the same p,q
// and then for "lower", we just shift the p,q after

void Grid::fill_field_lines(arma_vec baseLatsLoc,
                            precision_t min_altRe, precision_t Gamma,
                            Planets planet,
                            bool isCorner = false) {

  std::string function = "Grid::fill_field_lines";
  static int iFunction = -1;
  report.enter(function, iFunction);

  precision_t q_Start, delqp;

  // allocate & calculate some things outside of the main loop
  // - mostly just factors to make the code easier to read
  precision_t qp0, fb0, ft, delq, qp2, fa, fb, term0, term1, term2, term3;
  // exp_q_dist is the fraction of total q-distance to step for each pt along field line
  arma_vec exp_q_dist(nAlts);

  // corners/edges have one more lat dimension...
  int64_t nLatLoc = baseLatsLoc.n_elem;

  // temp holding of results from q,p -> r,theta conversion:
  std:: pair<precision_t, precision_t> r_theta;
  report.print(3, " calculating lshells!");

  // Find L-Shell for each baseLat
  // using L=R/sin2(theta), where theta is from north pole
  arma_vec Lshells(nLatLoc);

  for (int64_t iLat = 0; iLat < nLatLoc; iLat++)
    Lshells(iLat) = (min_altRe) / pow(sin(cPI / 2 - baseLatsLoc(iLat)), 2.0);

  report.print(3, "lshells calculated!");

  if (!isCorner) {
    for (int64_t iLon = 0; iLon < nLons; iLon ++) {
      for (int64_t iLat = 0; iLat < nLatLoc; iLat ++) {
        for (int64_t iAlt = 0; iAlt < nAlts; iAlt ++)
          magP_scgc(iLon, iLat, iAlt) = Lshells(iLat);
      }
    }
  } else {
    for (int64_t iLon = 0; iLon < nLons; iLon ++) {
      for (int64_t iLat = 0; iLat < nLatLoc; iLat ++) {
        for (int64_t iAlt = 0; iAlt < nAlts; iAlt ++)
          magP_Down(iLon, iLat, iAlt) = Lshells(iLat);
      }
    }
  }

  report.print(3, "dipole p-values stored for later.");

  for (int64_t iAlt = 0; iAlt < nAlts; iAlt++)
    exp_q_dist(iAlt) = Gamma + (1 - Gamma) * exp(-pow(((iAlt - nAlts) /
                                                       (nAlts / 5.0)), 2.0));

  report.print(3, "expQ");

  // mag alts and lats:
  arma_mat bAlts(nLatLoc, nAlts), bLats(nLatLoc, nAlts);

  for (int iLat = 0; iLat < nLatLoc; iLat++) {
    q_Start = -cos(cPI / 2 + baseLatsLoc(iLat)) / pow(min_altRe, 2.0);

    // calculate const stride in dipole coords, same as sami2/3 (huba & joyce 2000)
    // Note this is not the:
    // ==  >>   sinh(gamma*qi)/sinh(gamma*q_S)  <<  ==
    // but a different formula where the spacing is more easily controlled.
    // Doesn't have any lat/lon dependence so won't work for offset dipoles
    delqp = (-q_Start) / (nAlts + 1);
    delqp = min_altRe * delqp;

    for (int iAlt = 0; iAlt < nAlts; iAlt++) {
      qp0 = q_Start + iAlt * (delqp);
      fb0 = (1 - exp_q_dist(iAlt)) / exp(-q_Start / delqp - 1);
      ft = exp_q_dist(iAlt) - fb0 + fb0 * exp(-(qp0 - q_Start) / delqp);
      delq = qp0 - q_Start;

      // Q value at this point:
      qp2 = q_Start + ft * delq;

      if (isCorner) {
        // save the q for the "down" case:
        for (int64_t iLon = 0; iLon < nLons; iLon ++)
          magQ_Down(iLon, iLat, iAlt) = qp2;
      } else {
        for (int64_t iLon = 0; iLon < nLons; iLon ++)
          magQ_scgc(iLon, iLat, iAlt) = qp2;

        r_theta = qp_to_r_theta(qp2, Lshells(iLat));
        bAlts(iLat, iAlt) = r_theta.first;
        bLats(iLat, iAlt) = r_theta.second;
      }
    }
  }

  report.print(3, "QP-rtheta done!");

  if (isCorner) { // we don't need the rest, yet
    report.exit(function);
    return;
  }

  arma_vec rNorm1d(nAlts), lat1dAlong(nAlts);
  precision_t planetRadius;

  // rad_unit_vcgc = make_cube_vector(nLons, nLats, nAlts, 3);
  // This is wrong (same lat everywhere), but get_radius doesnt support oblate earth yet.
  planetRadius =  planet.get_radius(bLats.at(0));

  for (int64_t iLat = 0; iLat < nLatLoc; iLat++) {
    for (int64_t iLon = 0; iLon < nLons; iLon++) {
      // Not currently used. Dipole isn't offset. Leaving just in case.
      // Lon = magPhi_scgc(iLon, iLat, 1);

      for (int64_t iAlt = 0; iAlt < nAlts; iAlt++) {
        lat1dAlong(iAlt) = bLats(iLat, iAlt);
        rNorm1d(iAlt) = bAlts(iLat, iAlt);
      }

      // Lay things down in the same order as the geo grid.
      //centers only
      magAlt_scgc.tube(iLon,  iLat) = rNorm1d * planetRadius;
      magLat_scgc.tube(iLon, iLat) = lat1dAlong;
    }
  }

  report.exit(function);
  return;
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

  // Reverse our dipole rotations:
  xyzRot1 = rotate_around_y_3d(xyz_mag, magnetic_pole_tilt);
  xyzRot2 = rotate_around_z_3d(xyzRot1, magnetic_pole_rotation);

  // offset dipole (not yet implemented):
  // std::vector<precision_t> dipole_center = planet.get_dipole_center();
  // xyz_geo[0] = xyzRot2[0] + dipole_center[0];
  // xyz_geo[1] = xyzRot2[1] + dipole_center[1];
  // xyz_geo[2] = xyzRot2[2] + dipole_center[2];

  // transform back to lon, lat, radius:
  llr = transform_xyz_to_llr_3d(xyzRot2);

  report.exit(function);
  return llr;
}

// Use magP and magQ to make alt edges:
// This does the heavy lifting for the edges & corners of the dipole grid.
// These will be 1/2 way btwn each q point, which is pretty close to evenly spaced.
// They will not, however, line up from one field line to the next.
// It's not going to be *too* hard to get the corners to line up, but it messes with the
// orthogonality too much for me to figure out right now.
void Grid::dipole_alt_edges(Planets planet, precision_t min_altRe) {

  std::string function = "Grid::dipole_alt_edges";
  static int iFunction = -1;
  report.enter(function, iFunction);

  // P-coordinates will be the same along alt coord, we saved p-vals when we made them
  // in the fill field line function.
  precision_t pTmp;

  for (int64_t iLon = 0; iLon < nLons; iLon++) {
    for (int64_t iLat = 0; iLat < nLats + 1; iLat++) {
      pTmp = magP_Down(iLon, iLat, 0);

      for (int64_t iAlt = 0; iAlt < nAlts; iAlt ++)
        magP_Corner(iLon, iLat, iAlt) = pTmp;
    }
  }

  // Here are some shortcuts that exploit the symmetry.
  // This is done by each coord so cases like offset dipoles or oblate planets are easier later

  // first, use the fact that p is the same along each field line (alt)
  for (int64_t iLon = 0; iLon < nLons + 1; iLon++) {
    for (int64_t iLat = 0; iLat < nLats + 1; iLat++)
      magP_Corner(iLon, iLat, nAlts) = magP_Corner(iLon, iLat, nAlts - 1);
  }

  // And final step, use the longitude symmetry.
  // It's fine, until the dipole is offset. then the entire fill_field_lines needs to be redone.
  for (int64_t iAlt = 0; iAlt < nAlts + 1; iAlt++) {
    for (int64_t iLat = 0; iLat < nLats + 1; iLat++)
      magP_Corner(nLons, iLat, iAlt) = magP_Corner(nLons - 1, iLat, iAlt);
  }

  // For q-coord we'll avg q_down (from different baseLat) above and below the point...
  // May need to change the dipole spacing func's to get this working exactly though.
  // With how the field line pts are currently put in, this ends up being quite a hassle.
  // Not to mention, there would be a corner at q=0 (so r=A_LOT).
  // Top and bottom-most corners take the same q-step as the previous cell.
  precision_t qTmp;

  for (int64_t iLon = 0; iLon < nLons; iLon++) {
    for (int64_t iLat = 0; iLat < nLats + 1; iLat++) {
      for (int64_t iAlt = 1; iAlt < nAlts; iAlt ++)
        magQ_Corner(iLon, iLat, iAlt) = (magQ_Down(iLon, iLat,
                                                   iAlt - 1) + magQ_Down(iLon, iLat, iAlt)) / 2;

      magQ_Corner(iLon, iLat, 0) = (2 * magQ_Corner(iLon, iLat, 1) - magQ_Corner(iLon,
                                    iLat, 2));
    }
  }

  // for last (alt) corner, take the same step as the prev corner to the highest center.
  // this will force the highest corner to be above the last center
  for (int64_t iLon = 0; iLon < nLons; iLon++) {
    for (int64_t iLat = 0; iLat < nLats; iLat++) {
      qTmp = 2 * magQ_Corner(iLon, iLat, nAlts - 1) - magQ_Corner(iLon, iLat,
                                                                  nAlts - 2);
      magQ_Corner(iLon, iLat, nAlts) = qTmp;
    }
  }

  // last lon corner, copy previous. It's the same!
  for (int64_t iAlt = 0; iAlt < nAlts + 1; iAlt ++) {
    for (int64_t iLat = 0; iLat < nLats + 1; iLat++)
      magQ_Corner(nLons, iLat, iAlt) = magQ_Corner(nLons - 1, iLat, iAlt);
  }

  // Now we have (p,q) coords corners, convert to lon/lat/alt and we r off to the races
  std::pair <arma_cube, arma_cube> rtheta;
  precision_t planetRadius;
  rtheta = qp_to_r_theta(magQ_Corner, magP_Corner);
  magLat_Corner = rtheta.second;

  // Change if the dipole is offset and/or planet is oblate:
  planetRadius =  planet.get_radius(magLat_scgc.at(1));
  magAlt_Corner = rtheta.first * planetRadius;

  report.exit(function);
  return;
}


// -----------------------------------------------------------------------
// Convert XyzDipole to XyzGeo
//
// -----------------------------------------------------------------------

void Grid::convert_dipole_geo_xyz(Planets planet, precision_t XyzDipole[3],
                                  precision_t XyzGeo[3]) {

  std::string function = "Grid::convert_dipole_geo_xyz";
  static int iFunction = -1;
  report.enter(function, iFunction);

  precision_t XyzRemoveShift[3];
  precision_t XyzRemoveTilt[3];
  precision_t XyzRemoveRot[3];

  // get planetary parameters
  precision_t magnetic_pole_tilt = planet.get_dipole_tilt();
  precision_t magnetic_pole_rotation = planet.get_dipole_rotation();
  precision_t radius = planet.get_radius(0.0);


  // get the dipole shift, but normalize it to equatorial radius
  precision_t dipole_center[3];
  std::vector<precision_t> temp_dipole_center = planet.get_dipole_center();

  if ((temp_dipole_center[0] != 0) or (temp_dipole_center[1] != 0) or
      (temp_dipole_center[2] != 0)) {
    report.print(0,
                 "Dipole center != 0, but that is not supported yet. Setting to 0!");
    temp_dipole_center = {0, 0, 0};

  }

  transform_float_vector_to_array(temp_dipole_center, dipole_center);

  dipole_center[0] = dipole_center[0] / radius;
  dipole_center[1] = dipole_center[1] / radius;
  dipole_center[2] = dipole_center[2] / radius;

  // Remove Tilt
  transform_rot_y(XyzDipole, magnetic_pole_tilt, XyzRemoveTilt);

  // Remove Rot
  transform_rot_z(XyzRemoveTilt, magnetic_pole_rotation, XyzRemoveRot);

  // Remove Shift
  vector_add(XyzRemoveRot, dipole_center, XyzGeo);

  report.exit(function);
  return;

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

  report.print(0, "Creating inter-node connections Grid");

  if (!Is0D & !Is1Dz)
    create_sphere_connection(quadtree_ion);

  report.print(0, "Creating Dipole Grid");

  report.print(3, "Getting mgrid_inputs inputs in dipole grid");

  Inputs::grid_input_struct grid_input = input.get_grid_inputs("ionGrid");

  // Number of ghost cells:
  int64_t nGCs = get_nGCs();

  // Get inputs:
  precision_t min_alt = grid_input.alt_min * cKMtoM;
  precision_t LatStretch = grid_input.LatStretch;
  precision_t Gamma = grid_input.FieldLineStretch;
  precision_t min_apex = grid_input.min_apex * cKMtoM;
  precision_t max_lat = grid_input.max_blat;

  // Normalize inputs to planet radius... (update when earth is oblate)
  precision_t planetRadius = planet.get_radius(0.0);
  // Altitude to begin modeling, normalized to planet radius
  precision_t min_alt_re = (min_alt + planetRadius) / planetRadius;
  precision_t min_apex_re = (min_apex + planetRadius) / planetRadius;

  if (LatStretch != 1) {
    report.error("LatStretch values =/= 1 are not yet supported!");
    DidWork = false;
  }

  if (nAlts % 2 != 0) {
    report.error("nAlts must be even!");
    DidWork = false;
  }

  if (min_alt >= min_apex) {
    report.error("min_apex must be more than min_alt");
    DidWork = false;
  }

  // Get some coordinates and sizes in normalized coordinates:
  arma_vec lower_left_norm = quadtree_ion.get_vect("LL"); // origin
  arma_vec size_right_norm = quadtree_ion.get_vect("SR"); // lon_lims
  arma_vec size_up_norm = quadtree_ion.get_vect("SU");    //[1] = lat_lims
  report.print(3, "Initializing (dipole) longitudes");

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

  // Longitudes (symmetric, for now):
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
      // left edges
      magLon_Left.subcube(0, iLat, iAlt, nLons, iLat, iAlt) = lon1dLeft;
    }
  }


  for (iAlt = 0; iAlt < nAlts + 1; iAlt++) {
    for (iLat = 0; iLat < nLats + 1; iLat++) {
      // Corners
      magLon_Corner.subcube(0, iLat, iAlt, nLons, iLat, iAlt) = lon1dLeft;
    }
  }

  report.print(3, "Done initializing longitudes, moving to latitude");

  ////////////////
  // Latitudes: //
  ////////////////

  // min_lat calculated from min_apex
  precision_t min_lat = acos(sqrt(1 / min_apex_re));

  // latitude of field line base:
  // todo: needs support for variable stretching. it's like, halfway there.
  arma_vec baseLats = baselat_spacing(size_up_norm[1], lower_left_norm[1],
                                      max_lat, min_lat, 1.0);

  // downward sides (latitude shifted by 1/2 step):
  // TODO: This only works for linear latitude spacing, which is all that's supported right now.
  // When the exponential spacing (or something else) is fixed, this needs updating.
  precision_t dlat;
  dlat = baseLats(1) - baseLats(0);

  // put one cell halfway btwn each base latitude, leave 1st and last cell for now...
  for (int64_t iLat = 1; iLat < nLats; iLat ++)
    baseLats_down(iLat) = baseLats(iLat - 1) + (dlat * 0.5);

  // Put in 1st and last cell. Done this way so it's easier to put in supercell or something else
  baseLats_down(0) = baseLats(0) - dlat * 0.5;
  baseLats_down(nLats) = baseLats(nLats - 1) + dlat * 0.5;

  report.print(3, "baselats done!");

  // latitude & altitude of points on field lines (2D)
  // Cell centers
  fill_field_lines(baseLats, min_apex_re, Gamma, planet);
  // Corners (final bool argument) tells function to place stuff in the corner.
  // This is only down for the "down" edges, where the base latitudes are different.
  fill_field_lines(baseLats_down, min_apex_re, Gamma, planet, true);

  report.print(4, "Field-aligned Edges");
  dipole_alt_edges(planet, min_alt_re);

  report.print(3,
               "Done generating symmetric latitude & altitude spacing in dipole.");



  std::vector <arma_cube> llr = mag_to_geo(magLon_scgc, magLat_scgc, magAlt_scgc,
                                           planet);

  geoLon_scgc = llr[0];
  geoLat_scgc = llr[1];
  geoAlt_scgc = llr[2] - planetRadius;
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
  // s.elem(find(s == 0)).ones();

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
  magAlt_scgc = magAlt_scgc - planetRadius;

  report.exit(function);
  return DidWork;
}
