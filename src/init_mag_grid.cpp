// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>
#include <fstream>

#include "../include/aether.h"

// ----------------------------------------------------------------------
// Routine to convert p and q to r and theta. Can be solved iteratively,
// or with approach from (Swisdak, 2006), who solved it analytically:
//  https://arxiv.org/pdf/physics/0606044
//
// ----------------------------------------------------------------------
std::pair<precision_t, precision_t> Grid::qp_to_r_theta(precision_t q, precision_t p)
{
  // return quanties
  precision_t r, theta;
  // Intermediate quantities:
  precision_t term0, term1, term2, term3;

  term0 = 256.0 / 27.0 * pow(q, 2.0) * pow(p, 4.0);
  term1 = pow((1.0 + sqrt(1.0 + term0)), 2.0 / 3.0);
  term2 = pow(term0, 1.0 / 3.0);
  term3 = 0.5 * pow(((pow(term1, 2) + term1 * term2 + pow(term2, 2)) / term1), 3.0 / 2.0);

  r = p * (4.0 * term3) / (1.0 + term3) / (1.0 + sqrt(2.0 * term3 - 1.0));

  // now that r is determined we can solve for theta
  // theta = asin(sqrt(r/p));
  theta = asin(q * pow(r, 2.0));

  return {r, theta};
}

// ----------------------------------------------------------------------
// Initialize the dipole grid.
// - inputs (min_apex, min_alt, LatStretch, FieldLineStretch, max_lat_dipole)
//   are read from input files. And, of course, the numbers of each coordinate
// - nLats must be even!!
// ----------------------------------------------------------------------
void Grid::init_dipole_grid(Quadtree quadtree, Planets planet)
{

  std::string function = "Grid::init_dipole_grid";
  static int iFunction = -1;
  report.enter(function, iFunction);

  // turn the switch on!
  IsGeoGrid = false;
  IsMagGrid = true;

  // Dimension iterators
  int64_t iLon, iLat, iAlt;

  report.print(0, "Creating Dipole Grid");

  report.print(3, "Getting mgrid_inputs inputs in dipole grid");

  Inputs::grid_input_struct grid_input = input.get_grid_inputs("ionGrid");

  // Number of ghost cells:
  int64_t nGCs = get_nGCs();

  // Get inputs
  report.print(3, "Setting inputs in dipole grid");
  precision_t min_apex = grid_input.min_apex * cKMtoM;
  precision_t min_alt = grid_input.alt_min * cKMtoM;
  precision_t LatStretch = grid_input.LatStretch;
  precision_t Gamma = grid_input.FieldLineStretch;
  precision_t max_lat = grid_input.max_lat_dipole;

  // Normalize to planet radius...
  precision_t planetRadius = planet.get_radius(0.0);
  // L-Shell of minimum field line, normalized to planet radius
  precision_t min_lshell = (min_apex + planetRadius) / planetRadius;
  // Altitude to begin modeling, normalized to planet radius
  precision_t min_alt_re = (min_alt + planetRadius) / planetRadius;

  // Get some coordinates and sizes in normalized coordinates:
  arma_vec lower_left_norm = quadtree.get_vect("LL");
  arma_vec size_right_norm = quadtree.get_vect("SR");
  arma_vec size_up_norm = quadtree.get_vect("SU");

  // LONGITUDES:
  // - Make a 1d vector
  // - copy it into the 3d cube
  precision_t dlon = size_right_norm(0) * cPI / (nLons - 2 * nGCs);
  precision_t lon0 = lower_left_norm(0) * cPI;
  arma_vec lon1d(nLons);

  // Quick 1D check for longitudes:
  // Same as geo_grid here...
  if (!HasXdim) dlon = 1.0 * cDtoR;

  for (iLon = 0; iLon < nLons; iLon++)
    lon1d(iLon) = lon0 + (iLon - nGCs + 0.5) * dlon;

  for (iLat = 0; iLat < nLats; iLat++)
  {
    for (iAlt = 0; iAlt < nAlts; iAlt++)
      magLon_scgc.subcube(0, iLat, iAlt, nLons - 1, iLat, iAlt) = lon1d;
  }

  geoLon_scgc = magLon_scgc;

  // LATITUDES:

  // min_latitude calculated from min_lshell (min_apex input)
  precision_t min_lat = get_lat_from_r_and_lshell(1.0, min_lshell);

  // Lay down baseLat spacing according to an exponential factor.
  // some intermediates:
  precision_t del_lat, blat_min_, blat_max_, tmp_lat;
  // Integers for field-line loops:
  // - nF=nLats/2 (so nLats MUST be even)
  // - nZ=nAlts*2 - (we fill up HALF field lines)
  // - nZby2 = nAlts
  // -> Plus, experinemtal support for altitude ghost cells...
  int64_t nF = (nLats) / 2 , nZ = (nAlts) * 2, nZby2 = (nAlts);
  // lShells and baseLats are currently set for southern hemisphere then mirrored
  arma_vec Lshells(nF), baseLats(nF);

  blat_min_ = cos(pow(min_lat, 1.0 / LatStretch));
  blat_max_ = cos(pow(max_lat, 1.0 / LatStretch));
  del_lat = (blat_max_ - blat_min_) / (nF - nGCs * 2.0);

  // now make sure the user used 1 or an even number for nLats
  if (nLats % 2 != 0)
  {
    if (!HasYdim) 
    {
      del_lat = 1.0 * cDtoR;
      report.print(0, "Running in single latitude dimension. Experinental!!");
      nF = 1;
    }
    else
      report.error("Cannot use odd nLats with dipole grid!");
  }

  // loop over all cells - everything including the ghost cells
  // -> This means some points will go over the pole (baseLat > +- 90 degrees)
  //    They're taken care of in the conversion to geographic coordinates.
  for (int i = 0; i < nF; i++)
  {
    // first put down "linear" spacing
    tmp_lat = blat_min_ + del_lat * (i - nGCs + 0.5);
    // then scale it according to the exponent & acos
    tmp_lat = pow(acos(tmp_lat), LatStretch);
    // place values in array backwards, South pole -> equator.
    baseLats(nF - i - 1) = -tmp_lat;
  }
  report.print(3, "Done setting base latitudes for dipole grid.");

  // Find L-Shell for each baseLat
  // using L=R/sin2(theta), where theta is from north pole
  for (int i = 0; i < nF; i++)
    Lshells(i) = (min_alt_re) / pow(sin(cPI / 2 - baseLats(i)), 2.0);

  // SPACING ALONG FIELD LINE //
  // Coordinates along the field line to begin modeling
  // - In dipole (p,q) coordinates
  // - North & south hemisphere base
  precision_t q_S, q_N;
  // constant stride, scaled later
  precision_t delqp;

  // allocate & calculate some things outside of the main loop
  // - mistly just factors to make the code easier to read
  precision_t qp0, fb0, ft, delq, qp2, fa, fb, term0, term1, term2, term3, new_r;
  // exp_q_dist is the fraction of total q distance to step for each pt along field line
  arma_vec exp_q_dist(nZ), q_vals(nZ);
  // stored mag. coords temporarily in bAlts and bLats.
  arma_mat bAlts(nF, nZ), bLats(nF, nZ);

  // temp holding of results from q,p -> r,theta conversion:
  std::pair<precision_t, precision_t> r_theta;

  for (int i = 0; i < nAlts; i++)
    exp_q_dist(i) = Gamma + (1 - Gamma) * exp(-pow(((i - nZby2) / (nZ / 10.0)), 2.0));

  for (int i_nF = 0; i_nF < nF; i_nF++)
  {
    // min/max q
    q_S = -cos(cPI / 2 + baseLats(i_nF)) / pow(min_alt_re, 2.0);
    q_N = -q_S;

    // calculate const. stride similar to sami2/3 (huba & joyce 2000)
    // ==  >>   sinh(gamma*qi)/sinh(gamma*q_S)  <<  ==
    // inlo loop thru southern hemisphere, mirror in  north.
    for (int i_nZ = 0; i_nZ < nAlts; i_nZ++)
    {
      // This won't work for offset dipoles.
      // Doesn't have any lat/lon dependence.
      delqp = (q_N - q_S) / (nZ + 1);
      qp0 = q_S + i_nZ * (delqp);
      delqp = min_alt_re * delqp;
      fb0 = (1 - exp_q_dist(i_nZ)) / exp(-q_S / delqp - 1);
      ft = exp_q_dist(i_nZ) - fb0 + fb0 * exp(-(qp0 - q_S) / delqp);
      delq = qp0 - q_S;

      // Q value at this point:
      qp2 = q_S + ft * delq;

      r_theta = qp_to_r_theta(qp2, Lshells(i_nF));
      bAlts(i_nF, i_nZ) = r_theta.first;
      bLats(i_nF, i_nZ) = r_theta.second;
    }
  }
  report.print(3, "Done generating q-spacing for dipole grid.");

  arma_vec rNorm1d(nAlts), lat1dAlong(nAlts);
  arma_cube r3d(nLons, nLats, nAlts);

  // rad_unit_vcgc = make_cube_vector(nLons, nLats, nAlts, 3);

  for (iLat = 0; iLat < nLats / 2; iLat++)
  {
    for (iLon = 0; iLon < nLons; iLon++)
    {
      // Not currently used. Dipole isn't offset. Leaving just in case.
      // Lon = magPhi_scgc(iLon, iLat, 1);

      for (iAlt = 0; iAlt < nAlts; iAlt++)
      {
        lat1dAlong(iAlt) = bLats(iLat, iAlt);
        rNorm1d(iAlt) = bAlts(iLat, iAlt);
      }
      r3d.tube(iLon, iLat) = rNorm1d * planetRadius;
      r3d.tube(iLon, nLats - iLat - 1) = rNorm1d * planetRadius;
      magLat_scgc.tube(iLon, iLat) = lat1dAlong;
      magLat_scgc.tube(iLon, nLats - iLat - 1) = -lat1dAlong;
    }
  }
  report.print(3, "Done generating symmetric latitude & altitude spacing in dipole.");

  geoLat_scgc = magLat_scgc;
  magAlt_scgc = r3d - planetRadius;
  geoAlt_scgc = magAlt_scgc;

  report.print(4, "Beginning coordinate transformations of the dipole grid.");
  // Calculate the radius, of planet
  fill_grid_radius(planet);

  // Figure out what direction is radial:
  rad_unit_vcgc = make_cube_vector(nLons, nLats, nAlts, 3);
  gravity_vcgc = make_cube_vector(nLons, nLats, nAlts, 3);

  for (int iV = 0; iV < 3; iV++)
  {
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

  std::vector<arma_cube> llr, xyz, xyzRot1, xyzRot2;
  llr.push_back(magLon_scgc);
  llr.push_back(magLat_scgc);
  llr.push_back(r3d);
  xyz = transform_llr_to_xyz_3d(llr);

  precision_t magnetic_pole_rotation = planet.get_dipole_rotation();
  precision_t magnetic_pole_tilt = planet.get_dipole_tilt();

  // Reverse our dipole rotations:
  xyzRot1 = rotate_around_y_3d(xyz, magnetic_pole_tilt);
  xyzRot2 = rotate_around_z_3d(xyzRot1, magnetic_pole_rotation);

  // transform back to lon, lat, radius:
  llr = transform_xyz_to_llr_3d(xyzRot2);

  geoLon_scgc = llr[0];
  geoLat_scgc = llr[1];
  geoAlt_scgc = llr[2] - planetRadius;
  report.print(4, "Done dipole -> geographic transformations for the dipole grid.");

  calc_alt_grid_spacing();
  report.print(4, "Done altitude spacing for the dipole grid.");

  // Calculate magnetic field and magnetic coordinates:
  fill_grid_bfield(planet);
  report.print(4, "Done filling dipole grid with b-field!");

  report.exit(function);
  return;
}
