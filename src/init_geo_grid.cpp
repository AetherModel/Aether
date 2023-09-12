// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>

#include "aether.h"
#include <math.h>

// ----------------------------------------------------------------------
// Create connectivity between the nodes for message passing for cubesphere
// ----------------------------------------------------------------------

void Grid::create_cubesphere_connection(Quadtree quadtree) {

  std::string function = "Grid::create_cubesphere_connection";
  static int iFunction = -1;
  report.enter(function, iFunction);

  IsLatLonGrid = false;

  arma_vec lower_left_norm = quadtree.get_vect("LL");
  arma_vec middle_norm = quadtree.get_vect("MID");
  arma_vec size_right_norm = quadtree.get_vect("SR");
  arma_vec size_up_norm = quadtree.get_vect("SU");

  // These points go off the edge to the next block in each direction:
  arma_vec down_norm = middle_norm - 0.51 * size_up_norm;
  arma_vec up_norm = middle_norm + 0.51 * size_up_norm;
  arma_vec left_norm = middle_norm - 0.51 * size_right_norm;
  arma_vec right_norm = middle_norm + 0.51 * size_right_norm;

  // Find those points in the quadtree to figure out which processor
  // they are on
  iProcYm = quadtree.find_point(down_norm) + iMember * nGrids;
  iProcYp = quadtree.find_point(up_norm) + iMember * nGrids;
  iProcXm = quadtree.find_point(left_norm) + iMember * nGrids;
  iProcXp = quadtree.find_point(right_norm) + iMember * nGrids;

  // Need to know which side the current block is on and which side each
  // of the blocks in the different directions is on.  Need this so we can
  // know how to unpack the variables after the message pass.
  iRoot = quadtree.find_root(middle_norm);
  iRootYm = quadtree.find_root(down_norm);
  iRootYp = quadtree.find_root(up_norm);
  iRootXm = quadtree.find_root(left_norm);
  iRootXp = quadtree.find_root(right_norm);

  // These should be the exact edge of the face.
  // The from and to processors should get these in the same place,
  // so they can be used to match which processor to send / receive info
  edge_Xp = middle_norm + size_right_norm / 2.0;
  edge_Xm = middle_norm - size_right_norm / 2.0;
  edge_Yp = middle_norm + size_up_norm / 2.0;
  edge_Ym = middle_norm - size_up_norm / 2.0;

  if (report.test_verbose(2))
    std::cout << "connectivity : "
              << "  iProc : " << iProc << "\n"
              << "  isnorth : " << DoesTouchNorthPole << "\n"
              << "  issouth : " << DoesTouchSouthPole << "\n"
              << "  iProcXp : " << iProcXp << "\n"
              << "  iProcYp : " << iProcYp << "\n"
              << "  iProcXm : " << iProcXm << "\n"
              << "  iProcYm : " << iProcYm << "\n"
              << "  iRoot   : " << iRoot << "\n"
              << "  iRootXp : " << iRootXp << "\n"
              << "  iRootYp : " << iRootYp << "\n"
              << "  iRootXm : " << iRootXm << "\n"
              << "  iRootYm : " << iRootYm << "\n";

  int64_t iProcSelf = quadtree.find_point(middle_norm);

  report.exit(function);
  return;
}

// ----------------------------------------------------------------------
// This function takes the normalized coordinates and makes latitude
// and longitude arrays from them.  It can do this for the corners or
// edges, depending on the offset.
// ----------------------------------------------------------------------

void fill_cubesphere_lat_lon_from_norms(Quadtree quadtree,
                                        arma_vec dr,
                                        arma_vec du,
                                        arma_vec ll,
                                        int64_t nGCs,
                                        precision_t left_off,
                                        precision_t down_off,
                                        arma_mat &lat2d,
                                        arma_mat &lon2d,
                                        arma_mat &refx,
                                        arma_mat &refy) {

  int64_t nX = lat2d.n_rows;
  int64_t nY = lat2d.n_cols;

  double xn, yn, zn, rn;
  double xp, yp, zp, rp, latp, lonp;

  double a = sqrt(3);

  arma_vec xyz, xyzn, xyz_wrapped;

  // Loop through each point and derive the coordinate
  for (int iDU = 0; iDU < nY; iDU++) {
    for (int iLR = 0; iLR < nX; iLR++) {

      // the offsets are so we can find cell centers, edges, and corners
      double iD = iDU - nGCs + down_off;
      double iL = iLR - nGCs + left_off;

      // This is the normalized coordinate:
      xyz = ll + dr * iL + du * iD;
      // Ghost cells could be off the edge, so wrap to other face:
      //xyz_wrapped = quadtree.wrap_point_cubesphere(xyz) * a;
      xyz_wrapped = xyz * a;
      // Normalize the coordinate to a unit vector:
      xyzn = normalise(xyz_wrapped);
      xp = xyzn(0);
      yp = xyzn(1);
      zp = xyzn(2);

      // Derive lat and lon from unit vector:
      latp = asin(zp);
      // offset for lon is to put the left edge of face 0 at 0 longitude:
      lonp = atan2(yp, xp) + 3 * cPI / 4;

      if (lonp > cTWOPI)
        lonp = lonp - cTWOPI;

      if (lonp < 0.0)
        lonp = lonp + cTWOPI;

      lat2d(iLR, iDU) = latp;
      lon2d(iLR, iDU) = lonp;

      // Identify sides, then apply correct transformation law
      // Face 1 to 4, equator faces with face 1 starting at the meridian
      // Face 5, North Pole (different from book def)
      // Face 6, South Pole (different from book def)
      // Note face number are subtracted by one to comply with
      // computer indexing
      // Lon are displaced by cPI/4 as coordinates are generated
      // with a right displacement of cPI/4
      if (quadtree.iSide == 1 - 1) {
        refx(iLR, iDU) = sqrt(3) / 3 * tan(lonp - cPI / 4.);
        refy(iLR, iDU) = sqrt(3) / 3 * tan(latp) / cos(lonp - cPI / 4.);
      } else if (quadtree.iSide == 2 - 1) {
        refx(iLR, iDU) = sqrt(3) / 3 * tan(lonp - cPI / 4. - cPI / 2.);
        refy(iLR, iDU) = sqrt(3) / 3 * tan(latp) / cos(lonp - cPI / 4. - cPI / 2.);
      } else if (quadtree.iSide == 3 - 1) {
        refx(iLR, iDU) = sqrt(3) / 3 * tan(lonp - cPI / 4. - cPI);
        refy(iLR, iDU) = sqrt(3) / 3 * tan(latp) / cos(lonp - cPI / 4. - cPI);
      } else if (quadtree.iSide == 4 - 1) {
        refx(iLR, iDU) = sqrt(3) / 3 * tan(lonp - cPI / 4 - 3 * cPI / 2.);
        refy(iLR, iDU) = sqrt(3) / 3 * tan(latp) / cos(lonp - cPI / 4. - 3 * cPI / 2.);
      } else if (quadtree.iSide == 5 - 1) {
        refx(iLR, iDU) = -sqrt(3) / 3 * sin(lonp - 3 * cPI / 4.) / tan(latp);
        refy(iLR, iDU) = -sqrt(3) / 3 * cos(lonp - 3 * cPI / 4.) / tan(latp);
      } else if (quadtree.iSide == 6 - 1) {
        refx(iLR, iDU) = sqrt(3) / 3 * sin(lonp - 3 * cPI / 4.) / tan(latp);
        refy(iLR, iDU) = -sqrt(3) / 3 * cos(lonp - 3 * cPI / 4.) / tan(latp);
      }
    }
  }

  return;
}

// ----------------------------------------------------------------------
// This function takes in lat-lon and reference xy coordinates to
// generate transformation and metric tensors
// ----------------------------------------------------------------------
void transformation_metrics(Quadtree quadtree,
                            arma_mat &lat2d,
                            arma_mat &lon2d,
                            arma_mat &refx,
                            arma_mat &refy,
                            arma_mat &A11,
                            arma_mat &A12,
                            arma_mat &A21,
                            arma_mat &A22,
                            arma_mat &A11_inv,
                            arma_mat &A12_inv,
                            arma_mat &A21_inv,
                            arma_mat &A22_inv,
                            arma_mat &g11_upper,
                            arma_mat &g12_upper,
                            arma_mat &g21_upper,
                            arma_mat &g22_upper,
                            arma_mat &sqrt_g) {
  int64_t nX = lat2d.n_rows;
  int64_t nY = lat2d.n_cols;
  // Assume R = 1 (since lat-lon/ xy generation assumes unit vect)
  double R = 1;
  double a = 1 / sqrt(3);
  double xref, yref, rref;
  double latp, lonp;
  double g;

  // Loop through each point and derive the coordinate
  for (int j = 0; j < nY; j++) {
    for (int i = 0; i < nX; i++) {
      xref = refx(i, j);
      yref = refy(i, j);
      rref = std::sqrt(xref * xref + yref * yref + a * a);

      latp = lat2d(i, j);
      lonp = lon2d(i, j);

      sqrt_g(i, j) = R * R * a / (rref * rref * rref);
      g = sqrt_g(i, j) * sqrt_g(i, j);

      // metric tensor with lower indices
      double front_factor = R * R / (rref * rref * rref * rref);
      double g11 = front_factor * (a * a + yref * yref);
      double g12 = -front_factor * xref * yref;
      double g21 = -front_factor * xref * yref;
      double g22 = front_factor * (a * a + xref * xref);

      // metric tensor with upper indices
      g11_upper(i, j) = g22 / g;
      g12_upper(i, j) = -g12 / g;
      g21_upper(i, j) = -g21 / g;
      g22_upper(i, j) = g11 / g;

      // Identify sides, then apply correct transformation law
      // Face 1 to 4, equator faces with face 1 starting at the meridian
      // Face 5, North Pole (different from book def)
      // Face 6, South Pole (different from book def)
      // Note face number are subtracted by one to comply with
      // computer indexing
      if (quadtree.iSide == 1 - 1) {
        double p1 = R * cos(latp) * cos(lonp - cPI / 4.) / a;
        A11(i, j) = p1 * cos(lonp - cPI / 4.);
        A12(i, j) = 0;
        A21(i, j) = -p1 * sin(latp) * sin(lonp - cPI / 4.);
        A22(i, j) = p1 * cos(latp);

        double p2 = a / cos(latp) / cos(lonp - cPI / 4.) / R;
        A11_inv(i, j) = p2 / cos(lonp - cPI / 4.);
        A12_inv(i, j) = 0;
        A21_inv(i, j) = p2 * tan(latp) * tan(lonp - cPI / 4.);
        A22_inv(i, j) = p2 / cos(latp);
      } else if (quadtree.iSide == 2 - 1) {
        double p1 = R * cos(latp) * cos(lonp - cPI / 4. - cPI / 2.) / a;
        A11(i, j) = p1 * cos(lonp - cPI / 4. - cPI / 2.);
        A12(i, j) = 0;
        A21(i, j) = -p1 * sin(latp) * sin(lonp - cPI / 4. - cPI / 2.);
        A22(i, j) = p1 * cos(latp);

        double p2 = a / cos(latp) / cos(lonp - cPI / 4. - cPI / 2.) / R;
        A11_inv(i, j) = p2 / cos(lonp - cPI / 4. - cPI / 2.);
        A12_inv(i, j) = 0;
        A21_inv(i, j) = p2 * tan(latp) * tan(lonp - cPI / 4. - cPI / 2.);
        A22_inv(i, j) = p2 / cos(latp);
      } else if (quadtree.iSide == 3 - 1) {
        double p1 = R * cos(latp) * cos(lonp - cPI / 4. - cPI) / a;
        A11(i, j) = p1 * cos(lonp - cPI / 4. - cPI);
        A12(i, j) = 0;
        A21(i, j) = -p1 * sin(latp) * sin(lonp - cPI / 4. - cPI);
        A22(i, j) = p1 * cos(latp);

        double p2 = a / cos(latp) / cos(lonp - cPI / 4. - cPI) / R;
        A11_inv(i, j) = p2 / cos(lonp - cPI / 4. - cPI);
        A12_inv(i, j) = 0;
        A21_inv(i, j) = p2 * tan(latp) * tan(lonp - cPI / 4. - cPI);
        A22_inv(i, j) = p2 / cos(latp);
      } else if (quadtree.iSide == 4 - 1) {
        double p1 = R * cos(latp) * cos(lonp - cPI / 4. - 3 * cPI / 2.) / a;
        A11(i, j) = p1 * cos(lonp - cPI / 4. - 3 * cPI / 2.);
        A12(i, j) = 0;
        A21(i, j) = -p1 * sin(latp) * sin(lonp - cPI / 4. - 3 * cPI / 2.);
        A22(i, j) = p1 * cos(latp);

        double p2 = a / cos(latp) / cos(lonp - cPI / 4. - 3 * cPI / 2.) / R;
        A11_inv(i, j) = p2 / cos(lonp - cPI / 4. - 3 * cPI / 2.);
        A12_inv(i, j) = 0;
        A21_inv(i, j) = p2 * tan(latp) * tan(lonp - cPI / 4. - 3 * cPI / 2.);
        A22_inv(i, j) = p2 / cos(latp);
      } else if (quadtree.iSide == 5 - 1) {
        double p1 = R * sin(latp) / a;
        A11(i, j) = p1 * cos(lonp - 3 * cPI / 4.);
        A12(i, j) = p1 * sin(lonp - 3 * cPI / 4.);
        A21(i, j) = -p1 * sin(latp) * sin(lonp - 3 * cPI / 4.);
        A22(i, j) = p1 * sin(latp) * cos(lonp - 3 * cPI / 4.);

        double p2 = a / R / sin(latp) / sin(latp);
        A11_inv(i, j) = p2 * sin(latp) * cos(lonp - 3 * cPI / 4.);
        A12_inv(i, j) = -p2 * sin(lonp - 3 * cPI / 4.);
        A21_inv(i, j) = p2 * sin(latp) * sin(lonp - 3 * cPI / 4.);
        A22_inv(i, j) = p2 * cos(lonp - 3 * cPI / 4.);
      } else if (quadtree.iSide == 6 - 1) {
        double p1 = R * sin(latp) / a;
        A11(i, j) = -p1 * cos(lonp - 3 * cPI / 4.);
        A12(i, j) = p1 * sin(lonp - 3 * cPI / 4.);
        A21(i, j) = p1 * sin(latp) * sin(lonp - 3 * cPI / 4.);
        A22(i, j) = p1 * sin(latp) * cos(lonp - 3 * cPI / 4.);

        double p2 = a / R / sin(latp) / sin(latp);
        A11_inv(i, j) = -p2 * sin(latp) * cos(lonp - 3 * cPI / 4.);
        A12_inv(i, j) = p2 * sin(lonp - 3 * cPI / 4.);
        A21_inv(i, j) = p2 * sin(latp) * sin(lonp - 3 * cPI / 4.);
        A22_inv(i, j) = p2 * cos(lonp - 3 * cPI / 4.);
      }
    }
  }
}

// ----------------------------------------------------------------------
// Create a geographic grid
//    - if restarting, read in the grid
//    - if not restarting, initialize the grid
// ----------------------------------------------------------------------

void Grid::create_cubesphere_grid(Quadtree quadtree) {

  std::string function = "Grid::create_cubesphere_grid";
  static int iFunction = -1;
  report.enter(function, iFunction);

  arma_vec dr(3), du(3), ll(3);
  double xn, yn, zn, rn;
  double xp, yp, zp, rp, latp, lonp;

  double a = sqrt(3);

  arma_vec lower_left_norm = quadtree.get_vect("LL");
  arma_vec size_right_norm = quadtree.get_vect("SR");
  arma_vec size_up_norm = quadtree.get_vect("SU");

  dr = size_right_norm / (nLons - 2 * nGCs);
  du = size_up_norm / (nLats - 2 * nGCs);
  ll = lower_left_norm;

  int64_t iAlt, iLon, iLat;

  // ---------------------------------------------
  // Cell Centers
  // ---------------------------------------------
  arma_mat lat2d(nLons, nLats);
  arma_mat lon2d(nLons, nLats);
  arma_mat refx(nLons, nLats);
  arma_mat refy(nLons, nLats);
  arma_mat A11(nLons, nLats);
  arma_mat A12(nLons, nLats);
  arma_mat A21(nLons, nLats);
  arma_mat A22(nLons, nLats);
  arma_mat A11_inv(nLons, nLats);
  arma_mat A12_inv(nLons, nLats);
  arma_mat A21_inv(nLons, nLats);
  arma_mat A22_inv(nLons, nLats);
  arma_mat g11_upper(nLons, nLats);
  arma_mat g12_upper(nLons, nLats);
  arma_mat g21_upper(nLons, nLats);
  arma_mat g22_upper(nLons, nLats);
  arma_mat sqrt_g(nLons, nLats);
  fill_cubesphere_lat_lon_from_norms(quadtree, dr, du, ll, nGCs, 0.5, 0.5,
                                     lat2d, lon2d, refx, refy);

  transformation_metrics(quadtree, lat2d, lon2d, refx, refy,
                         A11, A12, A21, A22, A11_inv, A12_inv,
                         A21_inv, A22_inv, g11_upper, g12_upper,
                         g21_upper, g22_upper, sqrt_g);

  for (iAlt = 0; iAlt < nAlts; iAlt++) {
    geoLon_scgc.slice(iAlt) = lon2d;
    geoLat_scgc.slice(iAlt) = lat2d;
    refx_scgc.slice(iAlt) = refx;
    refy_scgc.slice(iAlt) = refy;
    A11_scgc.slice(iAlt) = A11;
    A12_scgc.slice(iAlt) = A12;
    A21_scgc.slice(iAlt) = A21;
    A22_scgc.slice(iAlt) = A22;
    A11_inv_scgc.slice(iAlt) = A11_inv;
    A12_inv_scgc.slice(iAlt) = A12_inv;
    A21_inv_scgc.slice(iAlt) = A21_inv;
    A22_inv_scgc.slice(iAlt) = A22_inv;
    g11_upper_scgc.slice(iAlt) = g11_upper;
    g12_upper_scgc.slice(iAlt) = g12_upper;
    g21_upper_scgc.slice(iAlt) = g21_upper;
    g22_upper_scgc.slice(iAlt) = g22_upper;
    sqrt_g_scgc.slice(iAlt) = sqrt_g;
  }

  // ---------------------------------------------
  // Left Sides - edges on left side (no offset left)
  // ---------------------------------------------
  arma_mat lat2d_left(nLons + 1, nLats);
  arma_mat lon2d_left(nLons + 1, nLats);
  arma_mat refx_left(nLons + 1, nLats);
  arma_mat refy_left(nLons + 1, nLats);
  arma_mat A11_left(nLons + 1, nLats);
  arma_mat A12_left(nLons + 1, nLats);
  arma_mat A21_left(nLons + 1, nLats);
  arma_mat A22_left(nLons + 1, nLats);
  arma_mat A11_inv_left(nLons + 1, nLats);
  arma_mat A12_inv_left(nLons + 1, nLats);
  arma_mat A21_inv_left(nLons + 1, nLats);
  arma_mat A22_inv_left(nLons + 1, nLats);
  arma_mat g11_upper_left(nLons + 1, nLats);
  arma_mat g12_upper_left(nLons + 1, nLats);
  arma_mat g21_upper_left(nLons + 1, nLats);
  arma_mat g22_upper_left(nLons + 1, nLats);
  arma_mat sqrt_g_left(nLons + 1, nLats);
  fill_cubesphere_lat_lon_from_norms(quadtree, dr, du, ll, nGCs, 0.0, 0.5,
                                     lat2d_left, lon2d_left,
                                     refx_left, refy_left);

  transformation_metrics(quadtree,
                         lat2d_left, lon2d_left, refx_left, refy_left,
                         A11_left, A12_left, A21_left, A22_left,
                         A11_inv_left, A12_inv_left,
                         A21_inv_left, A22_inv_left,
                         g11_upper_left, g12_upper_left,
                         g21_upper_left, g22_upper_left,
                         sqrt_g_left);

  for (iAlt = 0; iAlt < nAlts; iAlt++) {
    geoLon_Left.slice(iAlt) = lon2d_left;
    geoLat_Left.slice(iAlt) = lat2d_left;
    refx_Left.slice(iAlt) = refx_left;
    refy_Left.slice(iAlt) = refy_left;
    A11_Left.slice(iAlt) = A11_left;
    A12_Left.slice(iAlt) = A12_left;
    A21_Left.slice(iAlt) = A21_left;
    A22_Left.slice(iAlt) = A22_left;
    A11_inv_Left.slice(iAlt) = A11_inv_left;
    A12_inv_Left.slice(iAlt) = A12_inv_left;
    A21_inv_Left.slice(iAlt) = A21_inv_left;
    A22_inv_Left.slice(iAlt) = A22_inv_left;
    g11_upper_Left.slice(iAlt) = g11_upper_left;
    g12_upper_Left.slice(iAlt) = g12_upper_left;
    g21_upper_Left.slice(iAlt) = g21_upper_left;
    g22_upper_Left.slice(iAlt) = g22_upper_left;
    sqrt_g_Left.slice(iAlt) = sqrt_g_left;
  }

  // ---------------------------------------------
  // Down Sides - edges on down side (no offset down)
  // ---------------------------------------------
  arma_mat lat2d_down(nLons, nLats + 1);
  arma_mat lon2d_down(nLons, nLats + 1);
  arma_mat refx_down(nLons, nLats + 1);
  arma_mat refy_down(nLons, nLats + 1);
  arma_mat A11_down(nLons, nLats + 1);
  arma_mat A12_down(nLons, nLats + 1);
  arma_mat A21_down(nLons, nLats + 1);
  arma_mat A22_down(nLons, nLats + 1);
  arma_mat A11_inv_down(nLons, nLats + 1);
  arma_mat A12_inv_down(nLons, nLats + 1);
  arma_mat A21_inv_down(nLons, nLats + 1);
  arma_mat A22_inv_down(nLons, nLats + 1);
  arma_mat g11_upper_down(nLons, nLats + 1);
  arma_mat g12_upper_down(nLons, nLats + 1);
  arma_mat g21_upper_down(nLons, nLats + 1);
  arma_mat g22_upper_down(nLons, nLats + 1);
  arma_mat sqrt_g_down(nLons, nLats + 1);

  fill_cubesphere_lat_lon_from_norms(quadtree, dr, du, ll, nGCs, 0.5, 0.0,
                                     lat2d_down, lon2d_down,
                                     refx_down, refy_down);

  transformation_metrics(quadtree,
                         lat2d_down, lon2d_down, refx_down, refy_down,
                         A11_down, A12_down, A21_down, A22_down,
                         A11_inv_down, A12_inv_down,
                         A21_inv_down, A22_inv_down,
                         g11_upper_down, g12_upper_down,
                         g21_upper_down, g22_upper_down,
                         sqrt_g_down);

  for (iAlt = 0; iAlt < nAlts; iAlt++) {
    geoLon_Down.slice(iAlt) = lon2d_down;
    geoLat_Down.slice(iAlt) = lat2d_down;
    refx_Down.slice(iAlt) = refx_down;
    refy_Down.slice(iAlt) = refy_down;
    A11_Down.slice(iAlt) = A11_down;
    A12_Down.slice(iAlt) = A12_down;
    A21_Down.slice(iAlt) = A21_down;
    A22_Down.slice(iAlt) = A22_down;
    A11_inv_Down.slice(iAlt) = A11_inv_down;
    A12_inv_Down.slice(iAlt) = A12_inv_down;
    A21_inv_Down.slice(iAlt) = A21_inv_down;
    A22_inv_Down.slice(iAlt) = A22_inv_down;
    g11_upper_Down.slice(iAlt) = g11_upper_down;
    g12_upper_Down.slice(iAlt) = g12_upper_down;
    g21_upper_Down.slice(iAlt) = g21_upper_down;
    g22_upper_Down.slice(iAlt) = g22_upper_down;
    sqrt_g_Down.slice(iAlt) = sqrt_g_down;
  }

  // ---------------------------------------------
  // Corners (lower left) - no offsets
  // ---------------------------------------------
  arma_mat lat2d_corner(nLons + 1, nLats + 1);
  arma_mat lon2d_corner(nLons + 1, nLats + 1);
  arma_mat refx_corner(nLons + 1, nLats + 1);
  arma_mat refy_corner(nLons + 1, nLats + 1);
  fill_cubesphere_lat_lon_from_norms(quadtree, dr, du, ll, nGCs, 0.0, 0.0,
                                     lat2d_corner, lon2d_corner,
                                     refx_corner, refy_corner);

  for (iAlt = 0; iAlt < nAlts + 1; iAlt++) {
    geoLon_Corner.slice(iAlt) = lon2d_corner;
    geoLat_Corner.slice(iAlt) = lat2d_corner;
    refx_Corner.slice(iAlt) = refx_corner;
    refy_Corner.slice(iAlt) = refy_corner;
  }

  report.exit(function);
  return;
}

// ----------------------------------------------------------------------
// Create connectivity between the nodes for message passing for sphere
// ----------------------------------------------------------------------

void Grid::create_sphere_connection(Quadtree quadtree) {

  std::string function = "Grid::create_sphere_connection";
  static int iFunction = -1;
  report.enter(function, iFunction);

  IsLatLonGrid = true;

  // Get some coordinates and sizes in normalized coordinates:
  arma_vec lower_left_norm = quadtree.get_vect("LL");
  arma_vec middle_norm = quadtree.get_vect("MID");
  arma_vec size_right_norm = quadtree.get_vect("SR");
  arma_vec size_up_norm = quadtree.get_vect("SU");

  // Move to the next block in 4 directions:
  arma_vec down_norm = middle_norm - 0.51 * size_up_norm;
  arma_vec up_norm = middle_norm + 0.51 * size_up_norm;
  arma_vec left_norm = middle_norm - 0.51 * size_right_norm;
  arma_vec right_norm = middle_norm + 0.51 * size_right_norm;

  // The first component could wrap around:
  right_norm(0) = fmod(right_norm(0), quadtree.limit_high(0));
  left_norm(0) = fmod((left_norm(0) + quadtree.limit_high(0)),
                      quadtree.limit_high(0));

  // These should be the exact edge of the face.
  // The from and to processors should get these in the same place,
  // so they can be used to match which processor to send / receive info
  edge_Xp = middle_norm + size_right_norm / 2.0;
  // wrap in longitude:
  edge_Xp(0) = fmod(edge_Xp(0), quadtree.limit_high(0));
  edge_Xm = middle_norm - size_right_norm / 2.0;
  edge_Yp = middle_norm + size_up_norm / 2.0;
  edge_Ym = middle_norm - size_up_norm / 2.0;

  iProcYm = quadtree.find_point(down_norm) + iMember * nGrids;
  iProcYp = quadtree.find_point(up_norm) + iMember * nGrids;
  iProcXm = quadtree.find_point(left_norm) + iMember * nGrids;
  iProcXp = quadtree.find_point(right_norm) + iMember * nGrids;

  iRoot = quadtree.find_root(middle_norm);
  iRootYm = quadtree.find_root(down_norm);
  iRootYp = quadtree.find_root(up_norm);
  iRootXm = quadtree.find_root(left_norm);
  iRootXp = quadtree.find_root(right_norm);

  // Check if touching South Pole:
  if (lower_left_norm(1) == quadtree.limit_low(1)) {
    DoesTouchSouthPole = true;

    // edges need to be adjusted to deal with longitudes, since the
    // pole will 180deg different for the from and to processors
    if (edge_Ym(0) < 1.0)
      edge_Ym(0) += 0.5;
    else
      edge_Ym(0) -= 0.5;
  }

  // Check if touching North Pole:
  if (lower_left_norm(1) + size_up_norm(1) == quadtree.limit_high(1)) {
    DoesTouchNorthPole = true;

    // edge need to be adjusted to deal with longitudes, since the
    // pole will 180deg different for the from and to processors
    if (edge_Yp(0) < 1.0)
      edge_Yp(0) += 0.5;
    else
      edge_Yp(0) -= 0.5;
  }

  if (report.test_verbose(2))
    std::cout << "connectivity : "
              << "  iProc : " << iProc << "\n"
              << "  isnorth : " << DoesTouchNorthPole << "\n"
              << "  issouth : " << DoesTouchSouthPole << "\n"
              << "  iProcYm : " << iProcYm << "\n"
              << "  iProcYp : " << iProcYp << "\n"
              << "  iProcXm : " << iProcXm << "\n"
              << "  iProcXp : " << iProcXp << "\n";

  report.exit(function);
  return;
}

// ----------------------------------------------------------------------
// Create a spherical grid with lon/lat/alt coordinates
// ----------------------------------------------------------------------

void Grid::create_sphere_grid(Quadtree quadtree) {

  std::string function = "Grid::create_simple_lat_lon_alt_grid";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int64_t iLon, iLat, iAlt;

  // Get some coordinates and sizes in normalized coordinates:
  arma_vec lower_left_norm = quadtree.get_vect("LL");
  arma_vec size_right_norm = quadtree.get_vect("SR");
  arma_vec size_up_norm = quadtree.get_vect("SU");

  precision_t dlon = size_right_norm(0) * cPI / (nLons - 2 * nGCs);
  precision_t lon0 = lower_left_norm(0) * cPI;
  arma_vec lon1d(nLons);

  // Longitudes:
  // - Make a 1d vector
  // - copy it into the 3d cube
  for (iLon = 0; iLon < nLons; iLon++)
    lon1d(iLon) = lon0 + (iLon - nGCs + 0.5) * dlon;

  for (iLat = 0; iLat < nLats; iLat++) {
    for (iAlt = 0; iAlt < nAlts; iAlt++)
      geoLon_scgc.subcube(0, iLat, iAlt, nLons - 1, iLat, iAlt) = lon1d;
  }

  precision_t dlat = size_up_norm(1) * cPI / (nLats - 2 * nGCs);
  precision_t lat0 = lower_left_norm(1) * cPI;
  arma_vec lat1d(nLats);

  // Latitudes:
  // - Make a 1d vector
  // - copy it into the 3d cube
  for (iLat = 0; iLat < nLats; iLat++)
    lat1d(iLat) = lat0 + (iLat - nGCs + 0.5) * dlat;

  for (iLon = 0; iLon < nLons; iLon++) {
    for (iAlt = 0; iAlt < nAlts; iAlt++)
      geoLat_scgc.subcube(iLon, 0, iAlt, iLon, nLats - 1, iAlt) = lat1d;
  }

  // ---------------------------------------------
  // Left Sides - edges on left side (no offset left)
  // ---------------------------------------------
  arma_mat lat2d_left(nLons + 1, nLats);
  arma_mat lon2d_left(nLons + 1, nLats);

  for (iLat = 0; iLat < nLats; iLat++) {
    for (iLon = 0; iLon < nLons + 1; iLon++) {
      lat2d_left(iLon, iLat) = lat0 + (iLat - nGCs + 0.5) * dlat;
      lon2d_left(iLon, iLat) = lon0 + (iLon - nGCs) * dlon;
    }
  }

  for (iAlt = 0; iAlt < nAlts; iAlt++) {
    geoLon_Left.slice(iAlt) = lon2d_left;
    geoLat_Left.slice(iAlt) = lat2d_left;
  }

  // ---------------------------------------------
  // Down Sides - edges on down side (no offset lat)
  // ---------------------------------------------
  arma_mat lat2d_down(nLons, nLats + 1);
  arma_mat lon2d_down(nLons, nLats + 1);

  for (iLat = 0; iLat < nLats + 1; iLat++) {
    for (iLon = 0; iLon < nLons; iLon++) {
      lat2d_down(iLon, iLat) = lat0 + (iLat - nGCs) * dlat;
      lon2d_down(iLon, iLat) = lon0 + (iLon - nGCs + 0.5) * dlon;
    }
  }

  for (iAlt = 0; iAlt < nAlts; iAlt++) {
    geoLon_Down.slice(iAlt) = lon2d_down;
    geoLat_Down.slice(iAlt) = lat2d_down;
  }

  // ---------------------------------------------
  // Corner Sides - corner (no offset lat or lon)
  // ---------------------------------------------
  arma_mat lat2d_corner(nLons + 1, nLats + 1);
  arma_mat lon2d_corner(nLons + 1, nLats + 1);

  for (iLat = 0; iLat < nLats + 1; iLat++) {
    for (iLon = 0; iLon < nLons + 1; iLon++) {
      lat2d_corner(iLon, iLat) = lat0 + (iLat - nGCs) * dlat;
      lon2d_corner(iLon, iLat) = lon0 + (iLon - nGCs) * dlon;
    }
  }

  for (iAlt = 0; iAlt < nAlts + 1; iAlt++) {
    geoLon_Corner.slice(iAlt) = lon2d_corner;
    geoLat_Corner.slice(iAlt) = lat2d_corner;
  }

  report.exit(function);
  return;
}

// ----------------------------------------------------------------------
// Create a spherical grid with lon/lat/alt coordinates
// ----------------------------------------------------------------------

void Grid::create_altitudes(Planets planet) {

  std::string function = "Grid::create_altitudes";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int64_t iLon, iLat, iAlt;

  arma_vec alt1d(nAlts);

  Inputs::grid_input_struct grid_input = input.get_grid_inputs();

  if (grid_input.IsUniformAlt) {
    for (iAlt = 0; iAlt < nAlts; iAlt++)
      alt1d(iAlt) = grid_input.alt_min + (iAlt - nGeoGhosts) * grid_input.dalt;
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

    if (grid_input.dalt > 0.5) {
      if (report.test_verbose(0)) {
        std::cout << "-----------------------------------------------------\n";
        std::cout << "WARNING: dAlt is set to > 0.5, with non-uniform grid!\n";
        std::cout << "   dAlt = " << grid_input.dalt << "\n";
        std::cout << "-----------------------------------------------------\n";
      }
    }

    double alt = grid_input.alt_min;
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

    precision_t dalt = scale_height * grid_input.dalt;
    precision_t dAltLimiter = dalt * 10.0;

    // Fills bottom ghost cells with constant dAlt
    // Fills bottom cell with actual desired bottom altitude
    for (iAlt = 0; iAlt <= nGeoGhosts; iAlt++) {
      alt1d(iAlt) = grid_input.alt_min + (iAlt - nGeoGhosts) * dalt;

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

      dalt = scale_height * grid_input.dalt;

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
      geoAlt_Below.tube(iLon, iLat) = alt1d_below;
    }
  }

  for (iLon = 0; iLon < nLons + 1; iLon++) {
    for (iLat = 0; iLat < nLats + 1; iLat++)
      geoAlt_Corner.tube(iLon, iLat) = alt1d_below;
  }

  report.exit(function);
  return;
}

// ----------------------------------------------------------------------
// Corrects xy grid by scaling the R used in xy coordinate generation
// and transformation laws, as in previous generation R = 1.
// This function should only be used when cubesphere is used.
// Assumes radius of planet and altitude are constant
// ----------------------------------------------------------------------

void Grid::correct_xy_grid(Planets planet) {
  std::string function = "Grid::correct_xy_grid";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int64_t iAlt;

  // initialize grid drefx drefy
  drefx = arma_vec(nAlts); 
  drefy = arma_vec(nAlts);

  // Planet.get_radius() takes in latitude
  // but at current stage is unimplemented
  // Anyway, we use equator radius as assumption for CubeSphere
  // CubeSphere must be a perfect sphere!!
  precision_t planet_R = planet.get_radius(0);

  // radius of planet + altitude
  // just pick alt at (0,0) loction
  arma_vec R_Alts = geoAlt_scgc.tube(0, 0) + planet_R;

  for (iAlt = 0; iAlt < nAlts; iAlt++) {
    precision_t R = R_Alts(iAlt);
    refx_scgc.slice(iAlt) *= R;
    refy_scgc.slice(iAlt) *= R;
    A11_scgc.slice(iAlt) *= R;
    A12_scgc.slice(iAlt) *= R;
    A21_scgc.slice(iAlt) *= R;
    A22_scgc.slice(iAlt) *= R;
    A11_inv_scgc.slice(iAlt) /= R;
    A12_inv_scgc.slice(iAlt) /= R;
    A21_inv_scgc.slice(iAlt) /= R;
    A22_inv_scgc.slice(iAlt) /= R;
    g11_upper_scgc.slice(iAlt) /= R * R;
    g12_upper_scgc.slice(iAlt) /= R * R;
    g21_upper_scgc.slice(iAlt) /= R * R;
    g22_upper_scgc.slice(iAlt) /= R * R;
    sqrt_g_scgc.slice(iAlt) *= R * R;

    // Addition: Get a copy of dx dy
    arma_mat curr_refx = refx_scgc.slice(iAlt);
    arma_mat curr_refy = refy_scgc.slice(iAlt);

    drefx(iAlt) = curr_refx(1, 0) - curr_refx(0, 0);
    drefy(iAlt) = curr_refy(0, 1) - curr_refy(0, 0);

    refx_Left.slice(iAlt) *= R;
    refy_Left.slice(iAlt) *= R;
    A11_Left.slice(iAlt) *= R;
    A12_Left.slice(iAlt) *= R;
    A21_Left.slice(iAlt) *= R;
    A22_Left.slice(iAlt) *= R;
    A11_inv_Left.slice(iAlt) /= R;
    A12_inv_Left.slice(iAlt) /= R;
    A21_inv_Left.slice(iAlt) /= R;
    A22_inv_Left.slice(iAlt) /= R;
    g11_upper_Left.slice(iAlt) /= R * R;
    g12_upper_Left.slice(iAlt) /= R * R;
    g21_upper_Left.slice(iAlt) /= R * R;
    g22_upper_Left.slice(iAlt) /= R * R;
    sqrt_g_Left.slice(iAlt) *= R * R;

    refx_Down.slice(iAlt) *= R;
    refy_Down.slice(iAlt) *= R;
    A11_Down.slice(iAlt) *= R;
    A12_Down.slice(iAlt) *= R;
    A21_Down.slice(iAlt) *= R;
    A22_Down.slice(iAlt) *= R;
    A11_inv_Down.slice(iAlt) /= R;
    A12_inv_Down.slice(iAlt) /= R;
    A21_inv_Down.slice(iAlt) /= R;
    A22_inv_Down.slice(iAlt) /= R;
    g11_upper_Down.slice(iAlt) /= R * R;
    g12_upper_Down.slice(iAlt) /= R * R;
    g21_upper_Down.slice(iAlt) /= R * R;
    g22_upper_Down.slice(iAlt) /= R * R;
    sqrt_g_Down.slice(iAlt) *= R * R;

    refx_Corner.slice(iAlt) *= R;
    refy_Corner.slice(iAlt) *= R;
  }

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

  IsGeoGrid = 1;

  IsCubeSphereGrid = input.get_is_cubesphere();

  if (input.get_is_cubesphere())
    create_cubesphere_connection(quadtree);
  else
    create_sphere_connection(quadtree);

  if (input.get_do_restart()) {
    report.print(1, "Restarting! Reading grid files!");
    DidWork = read_restart(input.get_restartin_dir());
  } else {
    if (input.get_is_cubesphere())
      create_cubesphere_grid(quadtree);
    else
      create_sphere_grid(quadtree);

    MPI_Barrier(aether_comm);
    create_altitudes(planet);

    init_connection();

    DidWork = write_restart(input.get_restartout_dir());
  }

  // Calculate the radius (for spherical or non-spherical)
  fill_grid_radius(planet);
  // Calculate grid spacing
  calc_grid_spacing(planet);
  //calculate radial unit vector (for spherical or oblate planet)
  calc_rad_unit(planet);
  // Calculate gravity (including J2 term, if desired)
  calc_gravity(planet);

  // Calculate magnetic field and magnetic coordinates:
  fill_grid_bfield(planet);

  // Correct the reference grid with correct length scale:
  // (with R = actual radius)
  if (input.get_is_cubesphere())
    correct_xy_grid(planet);

  // Throw a little message for students:
  report.student_checker_function_name(input.get_is_student(),
                                       input.get_student_name(),
                                       4, "");

  report.exit(function);
  return DidWork;
}
