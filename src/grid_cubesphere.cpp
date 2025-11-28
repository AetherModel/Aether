// Copyright 2025, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md


#include "aether.h"

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

void Grid::init_cubesphere_grid(Quadtree quadtree,
                                arma_vec dr,
                                arma_vec du,
                                arma_vec ll,
                                precision_t left_off,
                                precision_t down_off,
                                cubesphere_chars &cubeX) {

  std::string function = "Grid::init_cubesphere_grid";
  static int iFunction = -1;
  report.enter(function, iFunction);

  precision_t dnu, dxi, nu0, xi0;

  // The du, dr, and ll were meant to be used on the cube
  // and not really on the equal-angle grid.  So, we probably
  // want to rethink these...
  if (quadtree.iSide == 0) {
    dnu = du[2];
    dxi = dr[0];
    nu0 = ll[0];
    xi0 = ll[2];
  }

  if (quadtree.iSide == 1) {
    dnu = du[2];
    dxi = dr[1];
    nu0 = ll[1];
    xi0 = ll[2];
  }

  if (quadtree.iSide == 2) {
    dnu = du[2];
    dxi = -dr[0];
    nu0 = -ll[0];
    xi0 = ll[2];
  }

  if (quadtree.iSide == 3) {
    dnu = du[2];
    dxi = -dr[1];
    nu0 = -ll[1];
    xi0 = ll[2];
  }

  if (quadtree.iSide == 4) {
    dnu = du[0];
    dxi = dr[1];
    nu0 = ll[1];
    xi0 = ll[0];
  }

  if (quadtree.iSide == 5) {
    dnu = -du[0];
    dxi = dr[1];
    nu0 = ll[1];
    xi0 = -ll[0];
  }

  // Normalized from -1 to 1 -> -pi/4 to pi/4
  dnu = dnu * cPI / 4.0;
  dxi = dxi * cPI / 4.0;
  nu0 = nu0 * cPI / 4.0;
  xi0 = xi0 * cPI / 4.0;

  cubeX.dnu = dnu;
  cubeX.dxi = dxi;

  int64_t iDU, iLR;
  precision_t iD, iL;
  int64_t nXp = nX, nYp = nY;

  // If we are shifting the grid over and doing edges, we
  // need to increase the number of points by 1 in that
  // direction:
  if (left_off < cSmall)
    nXp++;

  if (down_off < cSmall)
    nYp++;

  // These are convenient for the solver:
  cubeX.nXt = nXp;
  cubeX.nYt = nYp;
  cubeX.nGCs = nGCs;
  cubeX.iXfirst_ = nGCs;
  cubeX.iXlast_ = nXp - nGCs;
  cubeX.iYfirst_ = nGCs;
  cubeX.iYlast_ = nYp - nGCs;

  // these are coordinates:
  cubeX.lat.resize(nXp, nYp);
  cubeX.lon.resize(nXp, nYp);
  cubeX.nu.resize(nXp, nYp);
  cubeX.xi.resize(nXp, nYp);

  cubeX.X.resize(nXp, nYp);
  cubeX.Y.resize(nXp, nYp);
  cubeX.Z.resize(nXp, nYp);
  cubeX.C.resize(nXp, nYp);
  cubeX.D.resize(nXp, nYp);
  cubeX.d.resize(nXp, nYp);

  // These are dependent on radius,
  // but that is not included at this time:
  cubeX.dlx.resize(nXp, nYp, nZ);
  cubeX.dln.resize(nXp, nYp, nZ);
  cubeX.dS.resize(nXp, nYp, nZ);
  cubeX.R.resize(nZ);

  // These are matricies for rotating vectors:
  cubeX.Apn.resize(nXp, nYp);
  cubeX.Apx.resize(nXp, nYp);
  cubeX.Atn.resize(nXp, nYp);
  cubeX.Atx.resize(nXp, nYp);
  cubeX.Axt.resize(nXp, nYp);
  cubeX.Axp.resize(nXp, nYp);
  cubeX.Ant.resize(nXp, nYp);
  cubeX.Anp.resize(nXp, nYp);

  // These are for computing normals to the cell edges (horizontal)
  cubeX.nXiLon.resize(nXp, nYp);
  cubeX.nXiLat.resize(nXp, nYp);
  cubeX.nNuLon.resize(nXp, nYp);
  cubeX.nNuLat.resize(nXp, nYp);

  precision_t det, dmo, latp, lonp;

  // Loop through each point and derive the coordinate
  for (iDU = 0; iDU < nY; iDU++) {
    for (iLR = 0; iLR < nX; iLR++) {

      // the offsets are so we can find cell centers, edges, and corners
      iD = iDU - nGCs + down_off;
      iL = iLR - nGCs + left_off;

      // Define local coordinates:
      // Xi is LR (x), Nu is UD (y)
      cubeX.nu(iLR, iDU) = (nu0 + dnu * iD);
      cubeX.xi(iLR, iDU) = (xi0 + dxi * iL);

      cubeX.X(iLR, iDU) = tan(cubeX.xi(iLR, iDU));
      cubeX.Y(iLR, iDU) = tan(cubeX.nu(iLR, iDU));

      // Transformation from 3D Cartesian to LatLong
      // lonp = std::atan2(y_cart, x_cart) + cPI/2.0;
      if (quadtree.iSide == 0) {
        lonp = std::atan(cubeX.X(iLR, iDU));
        // Theta in Ronchi is from the north pole, so lat is 90 - theta
        latp = std::atan(1.0 / cubeX.Y(iLR, iDU) / std::cos(lonp));
      }

      if (quadtree.iSide == 1) {
        lonp = std::atan(-1.0 / cubeX.X(iLR, iDU));

        if (lonp < 0)
          lonp = cPI + lonp;

        // Theta in Ronchi is from the north pole, so lat is 90 - theta
        latp = std::atan(1.0 / cubeX.Y(iLR, iDU) / std::sin(lonp));
      }

      if (quadtree.iSide == 2) {
        lonp = std::atan(cubeX.X(iLR, iDU)) + cPI;
        // Theta in Ronchi is from the north pole, so lat is 90 - theta
        latp = std::atan(-1.0 / cubeX.Y(iLR, iDU) / std::cos(lonp));
      }

      if (quadtree.iSide == 3) {
        lonp = std::atan(-1.0 / cubeX.X(iLR, iDU));

        if (lonp > 0)
          lonp = lonp + cPI;
        else
          lonp = 2 * cPI + lonp;

        // Theta in Ronchi is from the north pole, so lat is 90 - theta
        latp = std::atan(-1.0 / cubeX.Y(iLR, iDU) / std::sin(lonp));
      }

      if (quadtree.iSide == 4) {
        lonp = std::atan2(cubeX.X(iLR, iDU), cubeX.Y(iLR, iDU));
        latp = std::atan2(-cubeX.Y(iLR, iDU), cos(lonp) );
      }

      if (quadtree.iSide == 5) {
        lonp = std::atan2(-cubeX.X(iLR, iDU), cubeX.Y(iLR, iDU));
        latp = -std::atan2(-cubeX.Y(iLR, iDU), cos(lonp) );
      }

      if (latp > 0)
        latp = cPI / 2 - latp;
      else
        latp = -(cPI / 2 + latp);

      if (lonp > cTWOPI)
        lonp = lonp - cTWOPI;

      if (lonp < 0.0)
        lonp = lonp + cTWOPI;

      // Fill Computed coords
      cubeX.lat(iLR, iDU) = latp;
      cubeX.lon(iLR, iDU) = lonp;

      cubeX.d(iLR, iDU) =
        1 +
        cubeX.X(iLR, iDU) * cubeX.X(iLR, iDU) +
        cubeX.Y(iLR, iDU) * cubeX.Y(iLR, iDU);

      cubeX.C(iLR, iDU) =
        sqrt(1 + cubeX.X(iLR, iDU) * cubeX.X(iLR, iDU));
      cubeX.D(iLR, iDU) =
        sqrt(1 + cubeX.Y(iLR, iDU) * cubeX.Y(iLR, iDU));

      if (quadtree.iSide < 4) {
        cubeX.Axt(iLR, iDU) = 0.0;
        cubeX.Axp(iLR, iDU) =
          cubeX.C(iLR, iDU) * cubeX.D(iLR, iDU) /
          sqrt(cubeX.d(iLR, iDU));
        cubeX.Ant(iLR, iDU) = -1.0;
        cubeX.Anp(iLR, iDU) =
          cubeX.X(iLR, iDU) * cubeX.Y(iLR, iDU) /
          sqrt(cubeX.d(iLR, iDU));
      } else {
        if (cubeX.d(iLR, iDU) < 1.0001)
          cubeX.d(iLR, iDU) = 1.0001;

        dmo = 1.0 / std::sqrt(cubeX.d(iLR, iDU) - 1);

        if (quadtree.iSide == 4) {
          cubeX.Axt(iLR, iDU) =
            - dmo * cubeX.D(iLR, iDU) * cubeX.X(iLR, iDU);
          cubeX.Axp(iLR, iDU) =
            dmo * cubeX.D(iLR, iDU) * cubeX.Y(iLR, iDU) /
            sqrt(cubeX.d(iLR, iDU));
          cubeX.Ant(iLR, iDU) =
            - dmo * cubeX.C(iLR, iDU) * cubeX.Y(iLR, iDU);
          cubeX.Anp(iLR, iDU) =
            - dmo * cubeX.C(iLR, iDU) * cubeX.X(iLR, iDU) /
            sqrt(cubeX.d(iLR, iDU));

        } else {
          // iFace == 5
          cubeX.Axt(iLR, iDU) =
            dmo * cubeX.D(iLR, iDU) * cubeX.X(iLR, iDU);
          cubeX.Axp(iLR, iDU) =
            - dmo * cubeX.D(iLR, iDU) *
            cubeX.Y(iLR, iDU) /
            sqrt(cubeX.d(iLR, iDU));
          cubeX.Ant(iLR, iDU) =
            dmo * cubeX.C(iLR, iDU) * cubeX.Y(iLR, iDU);
          cubeX.Anp(iLR, iDU) =
            dmo * cubeX.C(iLR, iDU) *
            cubeX.X(iLR, iDU) /
            sqrt(cubeX.d(iLR, iDU));
        }
      }

      // Calculate inverse of matrix for calculating Ax and An from At and Ap:
      det = 1.0 / (cubeX.Axt(iLR, iDU) * cubeX.Anp(iLR, iDU) -
                   cubeX.Axp(iLR, iDU) * cubeX.Ant(iLR, iDU));

      cubeX.Atx(iLR, iDU) = det * cubeX.Anp(iLR, iDU);
      cubeX.Atn(iLR, iDU) = - det * cubeX.Axp(iLR, iDU);
      cubeX.Apx(iLR, iDU) = - det * cubeX.Ant(iLR, iDU);
      cubeX.Apn(iLR, iDU) = det * cubeX.Axt(iLR, iDU);

      // These (dlx and dln) need to be multiplied by radius
      cubeX.dlx(iLR, iDU, 0) =
        cubeX.D(iLR, iDU) * dxi /
        cubeX.d(iLR, iDU) /
        (cos(cubeX.xi(iLR, iDU)) * cos(cubeX.xi(iLR, iDU)));
      cubeX.dln(iLR, iDU, 0) =
        cubeX.C(iLR, iDU) * dnu /
        cubeX.d(iLR, iDU) /
        (cos(cubeX.nu(iLR, iDU)) * cos(cubeX.nu(iLR, iDU)));

      // Need to multiply dS * radius ^ 2
      cubeX.dS(iLR, iDU, 0) =
        dxi * dnu /
        (sqrt(cubeX.d(iLR, iDU) * cubeX.d(iLR, iDU) * cubeX.d(iLR, iDU)) *
         cos(cubeX.xi(iLR, iDU)) * cos(cubeX.xi(iLR, iDU)) *
         cos(cubeX.nu(iLR, iDU)) * cos(cubeX.nu(iLR, iDU)));

    }
  }

  // Calculate norms given the values above:
  arma_mat e1Lat, e1Lon, e2Lat, e2Lon, m, one, zero;
  m.resize(nXp, nYp);
  one.resize(nXp, nYp);
  one.fill(1.0);
  zero.resize(nXp, nYp);
  zero.fill(0.0);

  // define e1 as the LR (xi) direction:
  e1Lat.resize(nXp, nYp);
  e1Lon.resize(nXp, nYp);
  convert_vector_xn_to_ll(one, zero, e1Lon, e1Lat, cubeX);
  m = sqrt(e1Lon % e1Lon + e1Lat % e1Lat);

  // Rotate by 90 deg (CCW) to get the norm:
  cubeX.nNuLon = -e1Lat / m;
  cubeX.nNuLat = e1Lon / m;

  // define e2 as the DU (nu) direction:
  e2Lat.resize(nXp, nYp);
  e2Lon.resize(nXp, nYp);
  convert_vector_xn_to_ll(zero, one, e2Lon, e2Lat, cubeX);
  m = sqrt(e2Lon % e2Lon + e2Lat % e2Lat);
  // Rotate by 90 deg (CW) to get the norm:
  cubeX.nXiLon = e2Lat / m;
  cubeX.nXiLat = -e2Lon / m;

  report.exit(function);
  return;
}

// ---------------------------------------------------------
// Convert vector from Alat, Alon to Axi, Anu
//   -> Using equation (7) of Ronchi et al:
// ---------------------------------------------------------

void Grid::convert_vector_xn_to_ll(arma_mat aXi,
                                   arma_mat aNu,
                                   arma_mat &aLon,
                                   arma_mat &aLat,
                                   cubesphere_chars grid) {

  // Ronchi defines aPhi = aLon, aTheta = -aLat
  aLat = -(grid.Atx % aXi + grid.Atn % aNu);
  aLon = grid.Apx % aXi + grid.Apn % aNu;

  return;
}

// ---------------------------------------------------------
// Convert vector from Alat, Alon to Axi, Anu
//   -> Using equation (7) of Ronchi et al:
// ---------------------------------------------------------

void Grid::convert_vector_ll_to_xn(arma_mat aLon,
                                   arma_mat aLat,
                                   arma_mat &aXi,
                                   arma_mat &aNu,
                                   cubesphere_chars grid) {

  // Ronchi defines aPhi = aLon, aTheta = -aLat
  aXi = -grid.Axt % aLat + grid.Axp % aLon;
  aNu = -grid.Ant % aLat + grid.Anp % aLon;
  return;
}



// ----------------------------------------------------------------------
// This function scales the deltas in the grid by the radius
//   - This assumes that radius is not dependent on lat / lon!!!
// ----------------------------------------------------------------------

void Grid::scale_cube_by_radius(cubesphere_chars &cubeX) {

  int64_t iZ;

  for (iZ = 1; iZ < nZ; iZ++) {
    cubeX.R(iZ) = radius_scgc(nGCs, nGCs, iZ);
    // These are distances:
    cubeX.dlx.slice(iZ) =
      cubeX.dlx.slice(0) * cubeX.R(iZ);
    cubeX.dln.slice(iZ) =
      cubeX.dln.slice(0) * cubeX.R(iZ);
    // This is an area:
    cubeX.dS.slice(iZ) =
      cubeX.dS.slice(0) * cubeX.R(iZ) * cubeX.R(iZ);
  }

  // Lastly, scale the 0th slice
  iZ = 0;
  cubeX.R(iZ) = radius_scgc(nGCs, nGCs, iZ);
  cubeX.dlx.slice(iZ) =
    cubeX.dlx.slice(0) * cubeX.R(iZ);
  cubeX.dln.slice(iZ) =
    cubeX.dln.slice(0) * cubeX.R(iZ);
  // This is an area:
  cubeX.dS.slice(iZ) =
    cubeX.dS.slice(0) * cubeX.R(iZ) * cubeX.R(iZ);

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
                                        arma_mat & lat2d,
                                        arma_mat & lon2d,
                                        arma_mat & refx,
                                        arma_mat & refy) {

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

      // refx and refy are the X, Y coordinates on the side of the CUBE
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
                            arma_mat & lat2d,
                            arma_mat & lon2d,
                            arma_mat & refx,
                            arma_mat & refy,
                            arma_mat & A11,
                            arma_mat & A12,
                            arma_mat & A21,
                            arma_mat & A22,
                            arma_mat & A11_inv,
                            arma_mat & A12_inv,
                            arma_mat & A21_inv,
                            arma_mat & A22_inv,
                            arma_mat & g11_upper,
                            arma_mat & g12_upper,
                            arma_mat & g21_upper,
                            arma_mat & g22_upper,
                            arma_mat & sqrt_g,
                            arma_mat & refx_angle,
                            arma_mat & refy_angle) {

  int64_t nX = lat2d.n_rows;
  int64_t nY = lat2d.n_cols;
  // Assume R = 1 (since lat-lon/ xy generation assumes unit vect)
  double R = 1;
  double a = R / sqrt(3);
  double xref, yref, rref, xy;
  double latp, lonp;
  double g;

  // Loop through each point and derive the coordinate
  for (int j = 0; j < nY; j++) {
    for (int i = 0; i < nX; i++) {
      xref = refx(i, j);
      yref = refy(i, j);
      xy = std::sqrt(xref * xref + yref * yref);
      rref = std::sqrt(xref * xref + yref * yref + a * a);

      // Want to calculate angles based on x, y, z
      refx_angle(i, j) = asin(xref / rref);
      refy_angle(i, j) = asin(yref / rref);

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
      } else if (quadtree.iSide == 6 -
                 1) { // Face 5 and 6 are flipped than Nair's formulation
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
      } else if (quadtree.iSide == 5 -
                 1) { // Face 5 and 6 are flipped than Nair's formulation
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

  // This function builds the equal-angle grid, but doesn't
  // scale them with altitude, since that has not been created, yet:
  init_cubesphere_grid(quadtree, dr, du, ll, 0.5, 0.5, cubeC);
  init_cubesphere_grid(quadtree, dr, du, ll, 0.0, 0.5, cubeL);
  init_cubesphere_grid(quadtree, dr, du, ll, 0.5, 0.0, cubeD);

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
  arma_mat refx_angle_temp(nLons, nLats);
  arma_mat refy_angle_temp(nLons, nLats);
  fill_cubesphere_lat_lon_from_norms(quadtree, dr, du, ll, nGCs, 0.5, 0.5,
                                     lat2d, lon2d, refx, refy);

  transformation_metrics(quadtree, lat2d, lon2d, refx, refy,
                         A11, A12, A21, A22, A11_inv, A12_inv,
                         A21_inv, A22_inv, g11_upper, g12_upper,
                         g21_upper, g22_upper, sqrt_g,
                         refx_angle_temp, refy_angle_temp);

  for (iAlt = 0; iAlt < nAlts; iAlt++) {
    refx_angle.slice(iAlt) = refx_angle_temp;
    refy_angle.slice(iAlt) = refy_angle_temp;

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
  arma_mat refx_angle_left_temp(nLons + 1, nLats);
  arma_mat refy_angle_left_temp(nLons + 1, nLats);
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
                         sqrt_g_left,
                         refx_angle_left_temp, refy_angle_left_temp);

  refx_angle_Left = refx_angle_left_temp;
  refy_angle_Left = refy_angle_left_temp;

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
  arma_mat refx_angle_down_temp(nLons, nLats + 1);
  arma_mat refy_angle_down_temp(nLons, nLats + 1);

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
                         sqrt_g_down,
                         refx_angle_down_temp, refy_angle_down_temp);
  refx_angle_Down = refx_angle_down_temp;
  refy_angle_Down = refy_angle_down_temp;

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

    // Addition: Get a copy of dx dy
    arma_mat curr_refx = refx_scgc.slice(iAlt);
    arma_mat curr_refy = refy_scgc.slice(iAlt);

    drefx(iAlt) = curr_refx(1, 0) - curr_refx(0, 0);
    drefy(iAlt) = curr_refy(0, 1) - curr_refy(0, 0);

    refx_Left.slice(iAlt) *= R;
    refy_Left.slice(iAlt) *= R;
    refx_Down.slice(iAlt) *= R;
    refy_Down.slice(iAlt) *= R;
    refx_Corner.slice(iAlt) *= R;
    refy_Corner.slice(iAlt) *= R;
  }

  report.exit(function);
  return;
}
