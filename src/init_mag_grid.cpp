// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>
#include <fstream>

#include "../include/aether.h"

// ----------------------------------------------------------------------
// Routine to find q_N and q_S for a given L
//
// ----------------------------------------------------------------------
std::pair<precision_t, precision_t> Grid::lshell_to_qn_qs(Planets planet,
                                                          precision_t Lshell,
                                                          precision_t Lon,
                                                          precision_t AltMin)
{
  std::string function = "Grid::lshell_to_qn_qs";
  static int iFunction = -1;
  report.enter(function, iFunction);

  precision_t qN, qS;

  precision_t XyzDipoleLeft[3], XyzDipoleMid[3], XyzDipoleRight[3];
  precision_t XyzGeoLeft[3], XyzGeoMid[3], XyzGeoRight[3];
  precision_t rGeoLeft, rGeoMid, rGeoRight;
  precision_t LlrDipoleLeft[3], LlrDipoleMid[3], LlrDipoleRight[3];
  precision_t ThetaTilt, PhiTilt;
  precision_t Lat, Radius, rMin;
  // Named dimension constants
  static int Lon_ = 0, Lat_ = 1, Radius_ = 2;

  // bound vars for bisection search
  precision_t ThetaRight, ThetaLeft, ThetaMid;
  precision_t rDipoleLeft, rDipoleMid, rDipoleRight;

  // Stopping condition for bisection search
  precision_t DeltaTheta;
  precision_t Tolerance = 1e-4;

  // status vars for bisection search
  int iStatusLeft, iStatusRight, iStatusMid;
  // note we normalize Lshell by equatorial radius
  precision_t RadiusEq = planet.get_radius(0.0);

  // loop for qN and qS
  for (int iQ = 0; iQ < 2; iQ++)
  {

    if (iQ == 0)
    {
      // set initial left, mid, right bounds for bisection search for qN
      ThetaRight = 0.5 * cPI;
      ThetaLeft = 1.0 * cDtoR;
      ThetaMid = 0.5 * (ThetaRight + ThetaLeft);
    }
    else
    {
      // set initial left, mid, right bounds for bisection search for qS
      ThetaLeft = 0.5 * cPI;
      ThetaRight = 179.0 * cDtoR;
      ThetaMid = 0.5 * (ThetaRight + ThetaLeft);
    }

    // Initial stopping condition stopping condition
    DeltaTheta = abs(ThetaLeft - ThetaRight);

    // start bisection search for qN
    while (DeltaTheta > Tolerance)
    {

      // find rDipole that cooresponds to these Left,Mid,Right
      // ThetaDipole values
      rDipoleLeft = Lshell * pow(sin(ThetaLeft), 2.0);
      rDipoleMid = Lshell * pow(sin(ThetaMid), 2.0);
      rDipoleRight = Lshell * pow(sin(ThetaRight), 2.0);

      // Compute XyzDipole for left, mid,right states
      LlrDipoleLeft[Lon_] = Lon;
      LlrDipoleLeft[Lat_] = 0.5 * cPI - ThetaLeft;
      LlrDipoleLeft[Radius_] = rDipoleLeft;
      transform_llr_to_xyz(LlrDipoleLeft, XyzDipoleLeft);

      LlrDipoleMid[Lon_] = Lon;
      LlrDipoleMid[Lat_] = 0.5 * cPI - ThetaMid;
      LlrDipoleMid[Radius_] = rDipoleMid;
      transform_llr_to_xyz(LlrDipoleMid, XyzDipoleMid);

      LlrDipoleRight[Lon_] = Lon;
      LlrDipoleRight[Lat_] = 0.5 * cPI - ThetaRight;
      LlrDipoleRight[Radius_] = rDipoleRight;
      transform_llr_to_xyz(LlrDipoleRight, XyzDipoleRight);

      // Transform to XyzGeo and unnormalize
      convert_dipole_geo_xyz(planet, XyzDipoleLeft, XyzGeoLeft);
      convert_dipole_geo_xyz(planet, XyzDipoleMid, XyzGeoMid);
      convert_dipole_geo_xyz(planet, XyzDipoleRight, XyzGeoRight);

      // cout << "XyzGeoLeft[0]" << XyzGeoLeft[0] << endl;
      // cout << "XyzGeoLeft[1]" << XyzGeoLeft[1] << endl;
      // cout << "XyzGeoLeft[2]" << XyzGeoLeft[2] << endl;

      XyzGeoLeft[0] = XyzGeoLeft[0] * RadiusEq;
      XyzGeoLeft[1] = XyzGeoLeft[1] * RadiusEq;
      XyzGeoLeft[2] = XyzGeoLeft[2] * RadiusEq;

      // abort;

      XyzGeoMid[0] = XyzGeoMid[0] * RadiusEq;
      XyzGeoMid[1] = XyzGeoMid[1] * RadiusEq;
      XyzGeoMid[2] = XyzGeoMid[2] * RadiusEq;

      XyzGeoRight[0] = XyzGeoRight[0] * RadiusEq;
      XyzGeoRight[1] = XyzGeoRight[1] * RadiusEq;
      XyzGeoRight[2] = XyzGeoRight[2] * RadiusEq;

      // Compute radius in geo coordinate for comparison to rmin
      rGeoLeft = sqrt(pow(XyzGeoLeft[0], 2) + pow(XyzGeoLeft[1], 2) + pow(XyzGeoLeft[2], 2));
      rGeoMid = sqrt(pow(XyzGeoMid[0], 2) + pow(XyzGeoMid[1], 2) + pow(XyzGeoMid[2], 2));
      rGeoRight = sqrt(pow(XyzGeoRight[0], 2) + pow(XyzGeoRight[1], 2) + pow(XyzGeoRight[2], 2));

      // get rmin for given latitude. Radius is lat dependent in general.
      // also find status in (0) or out (1) of rMin
      Lat = 0.5 * cPI - acos(XyzGeoLeft[2] / rGeoLeft);
      Radius = planet.get_radius(Lat);
      rMin = Radius + AltMin;
      if (rGeoLeft < rMin)
      {
        iStatusLeft = 0;
      }
      else
      {
        iStatusLeft = 1;
      }

      Lat = 0.5 * cPI - acos(XyzGeoMid[2] / rGeoMid);
      Radius = planet.get_radius(Lat);
      rMin = Radius + AltMin;
      if (rGeoMid < rMin)
      {
        iStatusMid = 0;
      }
      else
      {
        iStatusMid = 1;
      }

      Lat = 0.5 * cPI - acos(XyzGeoRight[2] / rGeoRight);
      Radius = planet.get_radius(Lat);
      rMin = Radius + AltMin;
      if (rGeoRight < rMin)
      {
        iStatusRight = 0;
      }
      else
      {
        iStatusRight = 1;
      }

      // Use status values to update left, right and mid values of theta
      if (iStatusMid == 0)
      {
        if (iStatusRight == 1)
        {
          // Mid becomes left and right stays right
          ThetaLeft = ThetaMid;
          ThetaMid = 0.5 * (ThetaLeft + ThetaRight);
        }
        else
        {
          // Mid becomes right and left stays left
          ThetaRight = ThetaMid;
          ThetaMid = 0.5 * (ThetaLeft + ThetaRight);
        }
      }
      else
      {
        if (iStatusRight == 0)
        {
          // Mid becomes left and right stays right
          ThetaLeft = ThetaMid;
          ThetaMid = 0.5 * (ThetaLeft + ThetaRight);
        }
        else
        {
          // Mid becomes right and left stays left
          ThetaRight = ThetaMid;
          ThetaMid = 0.5 * (ThetaLeft + ThetaRight);
        }
      }
      // Update stopping condition
      DeltaTheta = abs(ThetaLeft - ThetaRight);
    }

    // set the q value
    rDipoleMid = Lshell * pow(sin(ThetaMid), 2.0);
    if (iQ == 0)
    {
      qN = pow(rDipoleMid, -2.0) * cos(ThetaMid);
      // cout << "!!! For L = " << Lshell << endl;
      // cout << "!!! qN = " << qN << endl;
      // cout << "!!! ThetaMid = " << ThetaMid*cRtoD << endl;
    }
    else
    {
      qS = pow(rDipoleMid, -2.0) * cos(ThetaMid);
      // cout << "!!! qS = " << qS << endl;
    }
  }

  report.exit(function);
  return {qN, qS};
}

// -----------------------------------------------------------------------
// Convert XyzDipole to XyzGeo
//
// -----------------------------------------------------------------------

void Grid::convert_dipole_geo_xyz(Planets planet, precision_t XyzDipole[3], precision_t XyzGeo[3])
{
  precision_t XyzRemoveShift[3];
  precision_t XyzRemoveTilt[3];
  precision_t XyzRemoveRot[3];

  // get planetary parameters, use radius at equator for Lshell reference
  precision_t magnetic_pole_tilt = planet.get_dipole_tilt();
  precision_t magnetic_pole_rotation = planet.get_dipole_rotation();
  precision_t radius = planet.get_radius(0.0);

  // get the dipole shift, but normalize it to equatorial radius
  precision_t dipole_center[3];
  std::vector<precision_t> temp_dipole_center = planet.get_dipole_center();
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

  //  cout << "XyzDipole[0]" << XyzDipole[0] << endl;
  //  cout << "XyzDipole[1]" << XyzDipole[1] << endl;
  //  cout << "XyzDipole[2]" << XyzDipole[2] << endl;
  //
  //  cout << "XyzGeo[0]" << XyzGeo[0] << endl;
  //  cout << "XyzGeo[1]" << XyzGeo[1] << endl;
  //  cout << "XyzGeo[2]" << XyzGeo[2] << endl;
}

// ----------------------------------------------------------------------
// Routine to fill in the q values for a particular L and lon
// using equations 7-8 from Huba et al 2000
// ----------------------------------------------------------------------
std::pair<arma_vec, arma_vec> Grid::fill_dipole_q_line(precision_t qN_,
                                                       precision_t qS_,
                                                       precision_t Gamma_,
                                                       int64_t nZ_,
                                                       precision_t lShell_,
                                                       precision_t min_alt_)

{

  std::string function = "Grid::fill_dipole_q_line";
  static int iFunction = -1;
  report.enter(function, iFunction);

  precision_t r, theta;
  precision_t Dx;
  precision_t Llr[3], Xyz[3], q[nZ_];
  // precision_t Radius, rMin, alt_min, , delqp, qp0, fb0, ft;
  precision_t delqp, qp0, fb0, ft, delq;
  int64_t nZby2 = nZ_ / 2;
  arma_vec exp_q_dist(nZ_), r_vals(nZ_), lat_vals(nZ_);

  for (int i = 0; i < nZ_; i++)
  {
    exp_q_dist(i) = Gamma_ + (1 - Gamma_) * exp(-pow(((i - nZby2) / (nZ_ / 10.0)), 2.0));
  }
  for (int i_nZ = 0; i_nZ < nZ_; i_nZ++)
  {
    delqp = (qN_ - qS_) / nZ_;
    qp0 = qS_ + i_nZ * (delqp);
    delqp = min_alt_ * delqp;
    fb0 = (1 - exp_q_dist(i_nZ)) / exp(-qS_ / delqp - 1);
    ft = exp_q_dist(i_nZ) - fb0 + fb0 * exp(-(qp0 - qS_) / delqp);

    delq = qp0 - qS_;
    q[i_nZ] = qS_ + ft * delq;

    auto rtheta = qp_to_r_theta(-q[i_nZ], lShell_);
    r_vals[i_nZ] = rtheta.first;
    // r_vals[i_nZ + nZby2] = rtheta.first;
    lat_vals[i_nZ] = rtheta.second;
    // lat_vals[i_nZ + nZby2] = -rtheta.second;
  }

  report.exit(function);
  return {r_vals, lat_vals};
}

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
  for (iLon = 0; iLon < nLons; iLon++)
    lon1d(iLon) = lon0 + (iLon - nGCs + 0.5) * dlon;

  for (iLat = 0; iLat < nLats; iLat++)
  {
    for (iAlt = 0; iAlt < nAlts; iAlt++)
      magLon_scgc.subcube(0, iLat, iAlt, nLons - 1, iLat, iAlt) = lon1d;
  }

  geoLon_scgc = magLon_scgc;

  // LATITUDES

  // min_latitude calculated from min_lshell (min_apex input)
  precision_t min_lat = get_lat_from_r_and_lshell(1.0, min_lshell);

  // Lay down baseLat spacing according to an exponential factor.
  // some intermediates:
  precision_t del_lat, blat_min_, blat_max_, tmp_lat;
  // Integers for field-line loops:
  // - nF=nLats/2 (so nLarts MUST be even)
  // - nZ=nAlts
  // - nZby2 = nZ/2
  int64_t nF = nLats / 2, nZ = nAlts, nZby2 = nAlts / 2;
  // lShells and baseLats are currently set for southern hemisphere then mirrored
  arma_vec Lshells(nF), baseLats(nF);

  blat_min_ = cos(pow(min_lat, 1.0 / LatStretch));
  blat_max_ = cos(pow(max_lat, 1.0 / LatStretch));
  del_lat = (blat_max_ - blat_min_) / (nF - nGCs * 2.0);

  for (int i = 0; i < nF; i++)
  {
    // first put down "linear" spacing
    tmp_lat = blat_min_ + del_lat * (i - nGCs + 0.5);
    // then scale it according to the exponent & acos it
    tmp_lat = pow(acos(tmp_lat), LatStretch);
    // place values in array backwards, South pole -> equator.
    baseLats(nF - i - 1) = -tmp_lat;
  }

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

  for (int i = 0; i < nZ; i++)
    exp_q_dist(i) = Gamma + (1 - Gamma) * exp(-pow(((i - nZby2) / (nZ / 10.0)), 2.0));

  for (int i_nF = 0; i_nF < nF; i_nF++)
  {
    // min/max q
    q_S = -cos(cPI / 2 + baseLats(i_nF)) / pow(min_alt_re, 2.0);
    q_N = -q_S;

    // calculate const. stride similar to sami2/3 (huba & joyce 2000)
    // ==  >>   sinh(gamma*qi)/sinh(gamma*q_S)  <<  ==
    // inlo loop thru southern hemisphere, mirror in  north.
    for (int i_nZ = 0; i_nZ < nZ; i_nZ++)
    {
      // Should be using lshell_to_qn_qs, but it wasn't working right.
      // This does the same, but won't work for offset dipoles.
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

  geoLat_scgc = magLat_scgc;
  magAlt_scgc = r3d - planetRadius;
  geoAlt_scgc = magAlt_scgc;

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

  calc_alt_grid_spacing();

  // Calculate magnetic field and magnetic coordinates:
  fill_grid_bfield(planet);

  report.exit(function);
  return;
}
