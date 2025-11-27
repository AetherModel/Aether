/*
  This is an example of a second order 2D solver for the Euler equations.

  to compile:
  g++ -I/usr/local/include -I/Users/ridley/Software/Json/json/include -o cubesphere2d cubesphere2d.cpp

*/

#include <armadillo>
#include <fstream>

using precision_t = double;

/// Armadillo type vector (single column) with compile-time precision.
using arma_vec = arma::Col<precision_t>;
/// Armadillo type matrix (two dimension) with compile-time precision.
using arma_mat = arma::Mat<precision_t>;
/// Armadillo type cube (three dimension) with compile-time precision.
using arma_cube = arma::Cube<precision_t>;

precision_t cPI = 3.141592653589793;
precision_t cTWOPI = 2.0 * cPI;
precision_t cRtoD = 180.0 / cPI;
precision_t cPIdiv2 = cPI / 2;
precision_t cGamma = 5.0 / 3.0;  // Specific ratio of heat
precision_t cKb = 1.38e-23;
precision_t mmm = 16.0 * 1.67e-27;

// ---------------------------------------------------------
// A couple of global variables
// ---------------------------------------------------------

int64_t verbose = 1;

struct projection_struct {
  arma_mat gradLR;
  arma_mat gradDU;
  arma_mat R;
  arma_mat L;
  arma_mat U;
  arma_mat D;
};

struct grid_struct {

  // sizes:
  int64_t nXt, nYt, nGCs;
  int64_t iXfirst_, iXlast_;
  int64_t iYfirst_, iYlast_;

  // Positions:
  arma_mat lon;
  arma_mat lat;

  // These are for Ronchi et al., JCP 124, 93-114, 1996
  arma_mat X, Y, Z, C, D, d;
  arma_mat dlx, dln, dS;
  // xi is the LR direction
  // nu is the UD direction
  arma_mat xi, nu;
  arma_mat x, y, r;
  arma_mat Apn, Apx, Atn, Atx;
  arma_mat Axt, Axp, Ant, Anp;
  precision_t dxi, dnu, R;
  arma_mat alpha;
  arma_mat sinAlpha;

  arma_mat nXiLon;
  arma_mat nXiLat;
  arma_mat nNuLon;
  arma_mat nNuLat;
  // These are eq28 of Nair (g lower ij):
  arma_mat gl11, gl12, gl21, gl22;
  // These are eq29 of Nair (g upper ij):
  arma_mat sqrtg;
  arma_mat gu11, gu12, gu21, gu22;
  // These are eq32 of Nair (sphere-to-cube):
  arma_mat s2c11, s2c12, s2c21, s2c22;
  arma_mat c2s11, c2s12, c2s21, c2s22;
};

// ---------------------------------------------------------
//
// ---------------------------------------------------------

precision_t calc_dt(arma_mat dx,
                    arma_mat dy,
                    arma_mat &wsLR,
                    arma_mat &wsDU,
                    int64_t nGCs) {

  if (verbose > 2)
    std::cout << "  --> calc_dt\n";

  int64_t nX = wsLR.n_rows;
  int64_t nY = wsLR.n_cols;

  precision_t wsX, wsY, dtX, dtY, dt;

  dt = 1e32;

  for (int64_t j = nGCs; j < nY - nGCs; j++) {
    for (int64_t i = nGCs; i < nX - nGCs; i++) {
      wsX = (wsLR(i + 1, j) + wsLR(i, j)) / 2;
      dtX = dx(i, j) / wsX;
      wsY = (wsDU(i, j + 1) + wsDU(i, j)) / 2;
      dtY = dy(i, j) / wsY;

      if (dtX < dt)
        dt = dtX;

      if (dtY < dt)
        dt = dtY;
    }
  }

  return dt;
}

/**
 * Output function
 *
 * @param values Values
 * @param filename FileName
 * @param DoAppend
 */
void output(arma_mat &values,
            std::string filename,
            bool DoAppend) {

  std::ofstream outfile;

  if (DoAppend)
    outfile.open(filename, std::ios_base::app);
  else {
    outfile.open(filename);
    int64_t nX = values.n_rows;
    int64_t nY = values.n_cols;
    outfile << nX << " " << nY << "\n";
  }

  outfile << values;
  outfile.close();
}

/**
 * Transform spherical coordinates to 3D Cartesian
 *
 * doi: 10.1016/j.jcp.2007.07.022
 * Section 3, Eqn (23)
 *
 * @return dh Great Circle Distance between two points
 */
arma_vec sph2cart(precision_t lon,
                  precision_t lat,
                  precision_t r) {
  arma_vec xyz(3);
  xyz(0) = r * std::cos(lat) * std::cos(lon);
  xyz(1) = r * std::cos(lat) * std::sin(lon);
  xyz(2) = r * std::sin(lat);
  return xyz;
}

grid_struct init_grid_equidistant(int iFace,
                                  int64_t nX,
                                  int64_t nY,
                                  int64_t nGCs,
                                  precision_t R,
                                  precision_t xOff,
                                  precision_t yOff) {

  double a = R / std::sqrt(3);

  grid_struct grid;
  int64_t nXt = nX + 2 * nGCs;
  int64_t nYt = nY + 2 * nGCs;

  grid.nXt = nXt;
  grid.nYt = nYt;
  grid.nGCs = nGCs;
  grid.iXfirst_ = nGCs;
  grid.iYfirst_ = nGCs;
  grid.iXlast_ = nX + nGCs;
  grid.iYlast_ = nY + nGCs;

  // Positions:
  grid.lon.resize(nXt, nYt);
  grid.lat.resize(nXt, nYt);

  grid.x.resize(nXt, nYt);
  grid.y.resize(nXt, nYt);
  grid.r.resize(nXt, nYt);
  grid.xi.resize(nXt, nYt);
  grid.nu.resize(nXt, nYt);
  grid.X.resize(nXt, nYt);
  grid.Y.resize(nXt, nYt);
  grid.Z.resize(nXt, nYt);
  grid.C.resize(nXt, nYt);
  grid.D.resize(nXt, nYt);
  grid.d.resize(nXt, nYt);
  grid.dlx.resize(nXt, nYt);
  grid.dln.resize(nXt, nYt);
  grid.dS.resize(nXt, nYt);

  grid.s2c11.resize(nXt, nYt);
  grid.s2c12.resize(nXt, nYt);
  grid.s2c21.resize(nXt, nYt);
  grid.s2c22.resize(nXt, nYt);
  grid.c2s11.resize(nXt, nYt);
  grid.c2s12.resize(nXt, nYt);
  grid.c2s21.resize(nXt, nYt);
  grid.c2s22.resize(nXt, nYt);

  grid.gl11.resize(nXt, nYt);
  grid.gl12.resize(nXt, nYt);
  grid.gl21.resize(nXt, nYt);
  grid.gl22.resize(nXt, nYt);
  grid.gu11.resize(nXt, nYt);
  grid.gu12.resize(nXt, nYt);
  grid.gu21.resize(nXt, nYt);
  grid.gu22.resize(nXt, nYt);
  grid.sqrtg.resize(nXt, nYt);

  double iD, iL, x, y, theta, phi, r;

  double dx = 2 * a / nX;
  double dy = 2 * a / nY;
  double r3, r4;
  double R2 = R * R;

  // Loop through each point and derive the coordinate
  // DU is y-direction (down-up)
  // LR is x-direction (left-right)
  for (int iDU = 0; iDU < nYt; iDU++) {
    for (int iLR = 0; iLR < nXt; iLR++) {

      // the offsets are so we can find cell centers, edges, and corners
      // Centers assume Off = 0.5, which edges assume Off = 0
      double iD = iDU - nGCs + yOff;
      double iL = iLR - nGCs + xOff;

      x = -a + iL * dx;
      y = -a + iD * dy;
      phi = std::atan(x / a);
      // y = a * tan(theta) * sec(phi) => y * cos(phi) = a * tan(theta)
      theta = std::atan(y * std::cos(phi) / a);
      //std::cout << "Grid creation : "
      //          << iDU << " "
      //          << iLR << " "
      //          << y << " "
      //          << x << " "
      //          << theta * cRtoD << " "
      //          << phi * cRtoD << "\n";
      grid.lon(iLR, iDU) = phi;
      grid.lat(iLR, iDU) = theta;
      grid.X(iLR, iDU) = R * std::cos(theta) * std::cos(phi);
      grid.Y(iLR, iDU) = R * std::cos(theta) * std::sin(phi);
      grid.Z(iLR, iDU) = R * std::sin(theta);
      grid.r(iLR, iDU) = std::sqrt(a * a + x * x + y * y);

      // Equation 28 of Nair:
      r4 = r * r * r * r;
      grid.gl11(iLR, iDU) = R2 / r4 * (a * a + y * y);
      grid.gl12(iLR, iDU) = - R2 / r4 * (x * y);
      grid.gl21(iLR, iDU) = - R2 / r4 * (x * y);
      grid.gl22(iLR, iDU) = R2 / r4 * (a * a + x * x);
      // Equation 29 of Nair:
      r3 = r * r * r;
      grid.sqrtg(iLR, iDU) = R2 * a / r3;
      grid.gu11(iLR, iDU) = grid.gl22(iLR, iDU) / grid.sqrtg(iLR, iDU);
      grid.gu12(iLR, iDU) = - grid.gl12(iLR, iDU) / grid.sqrtg(iLR, iDU);
      grid.gu21(iLR, iDU) = - grid.gl21(iLR, iDU) / grid.sqrtg(iLR, iDU);
      grid.gu22(iLR, iDU) = grid.gl11(iLR, iDU) / grid.sqrtg(iLR, iDU);

      grid.s2c11(iLR, iDU) = a / (R * std::cos(theta) * std::cos(phi)) *
                             (1 / std::cos(theta));
      grid.s2c12(iLR, iDU) = 0.0;
      grid.s2c21(iLR, iDU) = a / (R * std::cos(theta) * std::cos(phi)) *
                             (std::tan(theta) * std::tan(phi));
      grid.s2c22(iLR, iDU) = a / (R * std::cos(theta) * std::cos(phi)) *
                             (1 / std::cos(phi));
      grid.c2s11(iLR, iDU) = (R * std::cos(theta) * std::cos(phi)) / a * std::cos(
                               theta);
      grid.c2s12(iLR, iDU) = 0.0;
      grid.c2s21(iLR, iDU) = -(R * std::cos(theta) * std::cos(phi)) / a *
                             std::sin(theta) * std::sin(phi);
      grid.c2s22(iLR, iDU) = (R * std::cos(theta) * std::cos(phi)) / a * std::cos(
                               phi);

    }
  }

  //xxx
  return grid;
}


grid_struct init_grid(int iFace,
                      int64_t nX, int64_t nY, int64_t nGCs,
                      precision_t R,
                      precision_t xOff, precision_t yOff) {

  grid_struct grid;
  int64_t nXt = nX + 2 * nGCs;
  int64_t nYt = nY + 2 * nGCs;

  grid.nXt = nXt;
  grid.nYt = nYt;
  grid.nGCs = nGCs;
  grid.iXfirst_ = nGCs;
  grid.iYfirst_ = nGCs;
  grid.iXlast_ = nX + nGCs;
  grid.iYlast_ = nY + nGCs;

  // Positions:
  grid.lon.resize(nXt, nYt);
  grid.lat.resize(nXt, nYt);

  grid.xi.resize(nXt, nYt);
  grid.nu.resize(nXt, nYt);
  grid.X.resize(nXt, nYt);
  grid.Y.resize(nXt, nYt);
  grid.C.resize(nXt, nYt);
  grid.D.resize(nXt, nYt);
  grid.d.resize(nXt, nYt);
  grid.dlx.resize(nXt, nYt);
  grid.dln.resize(nXt, nYt);
  grid.dS.resize(nXt, nYt);

  grid.Axt.resize(nXt, nYt);
  grid.Axp.resize(nXt, nYt);
  grid.Ant.resize(nXt, nYt);
  grid.Anp.resize(nXt, nYt);

  grid.Apn.resize(nXt, nYt);
  grid.Apx.resize(nXt, nYt);
  grid.Atn.resize(nXt, nYt);
  grid.Atx.resize(nXt, nYt);

  precision_t fortyfive = cPI / 4.0;
  // Xi is LR (x), Nu is UD (y)
  precision_t dxi = 2.0 * fortyfive / (nX - 1 + xOff * 2);
  precision_t dnu = 2.0 * fortyfive / (nY - 1 + yOff * 2);

  grid.dxi = dxi;
  grid.dnu = dnu;
  grid.R = R;

  precision_t latp, lonp;

  // Loop through each point and derive the coordinate

  precision_t total_area = 0.0, det, dmo;

  for (int iDU = 0; iDU < nYt; iDU++) {
    for (int iLR = 0; iLR < nXt; iLR++) {

      // the offsets are so we can find cell centers, edges, and corners
      double iD = iDU - nGCs + yOff;
      double iL = iLR - nGCs + xOff;

      // Define local coordinates:
      // Xi is LR (x), Nu is UD (y)
      grid.nu(iLR, iDU) = (-fortyfive + dnu * iD);
      grid.xi(iLR, iDU) = (-fortyfive + dxi * iL);

      grid.X(iLR, iDU) = tan(grid.xi(iLR, iDU));
      grid.Y(iLR, iDU) = tan(grid.nu(iLR, iDU));

      // Transformation from 3D Cartesian to LatLong
      // lonp = std::atan2(y_cart, x_cart) + cPI/2.0;
      if (iFace == 0) {
        lonp = std::atan(grid.X(iLR, iDU));
        // Theta in Ronchi is from the north pole, so lat is 90 - theta
        latp = std::atan(1.0 / grid.Y(iLR, iDU) / std::cos(lonp));
      }

      if (iFace == 1) {
        lonp = std::atan(-1.0 / grid.X(iLR, iDU));

        if (lonp < 0)
          lonp = cPI + lonp;

        // Theta in Ronchi is from the north pole, so lat is 90 - theta
        latp = std::atan(1.0 / grid.Y(iLR, iDU) / std::sin(lonp));
      }

      if (iFace == 2) {
        lonp = std::atan(grid.X(iLR, iDU)) + cPI;
        // Theta in Ronchi is from the north pole, so lat is 90 - theta
        latp = std::atan(-1.0 / grid.Y(iLR, iDU) / std::cos(lonp));
      }

      if (iFace == 3) {
        lonp = std::atan(-1.0 / grid.X(iLR, iDU));

        if (lonp > 0)
          lonp = lonp + cPI;
        else
          lonp = 2 * cPI + lonp;

        // Theta in Ronchi is from the north pole, so lat is 90 - theta
        latp = std::atan(-1.0 / grid.Y(iLR, iDU) / std::sin(lonp));
      }

      if (iFace == 4) {
        lonp = std::atan2(grid.X(iLR, iDU), grid.Y(iLR, iDU));
        latp = std::atan2(-grid.Y(iLR, iDU), cos(lonp) );
      }

      if (iFace == 5) {
        lonp = std::atan2(-grid.X(iLR, iDU), grid.Y(iLR, iDU));
        latp = -std::atan2(-grid.Y(iLR, iDU), cos(lonp) );
      }

      if (latp > 0)
        latp = cPI / 2 - latp;
      else
        latp = -(cPI / 2 + latp);

      // Fill Computed coords
      grid.lat(iLR, iDU) = latp;
      grid.lon(iLR, iDU) = lonp;

      grid.d(iLR, iDU) =
        1 +
        grid.X(iLR, iDU) * grid.X(iLR, iDU) +
        grid.Y(iLR, iDU) * grid.Y(iLR, iDU);
      grid.C(iLR, iDU) = sqrt(1 +
                              grid.X(iLR, iDU) * grid.X(iLR, iDU));
      grid.D(iLR, iDU) = sqrt(1 +
                              grid.Y(iLR, iDU) * grid.Y(iLR, iDU));

      if (iFace < 4) {
        grid.Axt(iLR, iDU) = 0.0;
        grid.Axp(iLR, iDU) = grid.C(iLR, iDU) * grid.D(iLR, iDU) /
                             sqrt(grid.d(iLR, iDU));
        grid.Ant(iLR, iDU) = -1.0;
        grid.Anp(iLR, iDU) = grid.X(iLR, iDU) * grid.Y(iLR, iDU) /
                             sqrt(grid.d(iLR, iDU));
      } else {
        dmo = 1.0 / std::sqrt(grid.d(iLR, iDU) - 1);

        //if (dmo > 100.0)
        //  dmo = 100.0;

        if (iFace == 4) {
          grid.Axt(iLR, iDU) = - dmo * grid.D(iLR, iDU) * grid.X(iLR, iDU);
          grid.Axp(iLR, iDU) = dmo * grid.D(iLR, iDU) * grid.Y(iLR, iDU) /
                               sqrt(grid.d(iLR, iDU));
          grid.Ant(iLR, iDU) = - dmo * grid.C(iLR, iDU) * grid.Y(iLR, iDU);
          grid.Anp(iLR, iDU) = - dmo * grid.C(iLR, iDU) * grid.X(iLR, iDU) /
                               sqrt(grid.d(iLR, iDU));

        } else {
          // iFace == 5
          grid.Axt(iLR, iDU) = dmo * grid.D(iLR, iDU) * grid.X(iLR, iDU);
          grid.Axp(iLR, iDU) = - dmo * grid.D(iLR, iDU) * grid.Y(iLR, iDU) /
                               sqrt(grid.d(iLR, iDU));
          grid.Ant(iLR, iDU) = dmo * grid.C(iLR, iDU) * grid.Y(iLR, iDU);
          grid.Anp(iLR, iDU) = dmo * grid.C(iLR, iDU) * grid.X(iLR, iDU) /
                               sqrt(grid.d(iLR, iDU));
        }
      }

      // Calculate inverse of matrix for calculating Ax and An from At and Ap:
      det = 1.0 / (grid.Axt(iLR, iDU) * grid.Anp(iLR, iDU) -
                   grid.Axp(iLR, iDU) * grid.Ant(iLR, iDU));

      grid.Atx(iLR, iDU) = det * grid.Anp(iLR, iDU);
      grid.Atn(iLR, iDU) = - det * grid.Axp(iLR, iDU);
      grid.Apx(iLR, iDU) = - det * grid.Ant(iLR, iDU);
      //grid.Apx(iLR, iDU) = grid.d(iLR, iDU) / (grid.C(iLR, iDU) * grid.D(iLR, iDU));
      grid.Apn(iLR, iDU) = det * grid.Axt(iLR, iDU);

      grid.dlx(iLR, iDU) =
        R * grid.D(iLR, iDU) * dxi /
        grid.d(iLR, iDU) /
        (cos(grid.xi(iLR, iDU)) * cos(grid.xi(iLR, iDU)));
      grid.dln(iLR, iDU) =
        R * grid.C(iLR, iDU) * dnu /
        grid.d(iLR, iDU) /
        (cos(grid.nu(iLR, iDU)) * cos(grid.nu(iLR, iDU)));

      //grid.dS(iLR, iDU) = R * R * dxi * dnu /
      //                    (sqrt(grid.d(iLR, iDU) * grid.d(iLR, iDU) * grid.d(iLR, iDU))) *
      //                    grid.C(iLR, iDU) * grid.C(iLR, iDU) *
      //                    grid.D(iLR, iDU) * grid.D(iLR, iDU);
      grid.dS(iLR, iDU) = R * R * dxi * dnu /
                          (sqrt(grid.d(iLR, iDU) * grid.d(iLR, iDU) * grid.d(iLR, iDU)) *
                           cos(grid.xi(iLR, iDU)) * cos(grid.xi(iLR, iDU)) *
                           cos(grid.nu(iLR, iDU)) * cos(grid.nu(iLR, iDU)));
      //grid.dS(iLR, iDU) = R * R * dxi * dnu *
      //                    grid.C(iLR, iDU) * grid.D(iLR, iDU) /
      //                    (grid.d(iLR, iDU) * grid.d(iLR, iDU) *
      //                     cos(grid.xi(iLR, iDU)) * cos(grid.xi(iLR, iDU)) *
      //                     cos(grid.nu(iLR, iDU)) * cos(grid.nu(iLR, iDU)));

      if (iLR > 2 && iLR < nXt - 3 &&
          iDU > 2 && iDU < nYt - 3)
        total_area = total_area + grid.dS(iLR, iDU);
    }
  }

  std::cout << "Total Area : " << total_area << "; expected : " << 4 * cPI * R *
            R / 6.0 << "\n";

  return grid;
}

/**
 * Calculate Great Circle Distance
 *
 * doi: 10.1016/j.jcp.2007.07.022
 * Section 3, Eqn (23)
 *
 * @return dh Great Circle Distance between two points
 */
precision_t calc_great_circle(precision_t lon1,
                              precision_t lon2,
                              precision_t lat1,
                              precision_t lat2) {

  precision_t dlon_2 = (lon2 - lon1) / 2.0;
  precision_t dlat_2 = (lat2 - lat1) / 2.0;

  precision_t dh = 2.0 * std::asin(std::sqrt(std::sin(dlat_2) * std::sin(dlat_2) +
                                             std::sin(dlon_2) * std::sin(dlon_2) * std::cos(lat1) * std::cos(lat2)));

  return dh;
}

// ---------------------------------------------------------
// Angle between three points on a sphere
// ---------------------------------------------------------
precision_t calc_angle_given_three_lon_lat(precision_t p1_lon,
                                           precision_t p1_lat,
                                           precision_t p2_lon,
                                           precision_t p2_lat,
                                           precision_t p3_lon,
                                           precision_t p3_lat,
                                           precision_t R) {

  arma_vec p1 = sph2cart(p1_lon, p1_lat, R);
  arma_vec p2 = sph2cart(p2_lon, p2_lat, R);
  arma_vec p3 = sph2cart(p3_lon, p3_lat, R);
  arma_vec d1 = p1 - p2;
  arma_vec d2 = p3 - p2;
  precision_t n1 = std::sqrt(d1(0) * d1(0) +
                             d1(1) * d1(1) +
                             d1(2) * d1(2));
  precision_t n2 = std::sqrt(d2(0) * d2(0) +
                             d2(1) * d2(1) +
                             d2(2) * d2(2));
  d1 = d1 / n1;
  d2 = d2 / n2;
  precision_t angle = std::acos(d1(0) * d2(0) +
                                d1(1) * d2(1) +
                                d1(2) * d2(2));

  return angle;
}

// ---------------------------------------------------------
// bin edges
// ---------------------------------------------------------

arma_vec calc_bin_edges(arma_vec centers) {

  int64_t nPts = centers.n_elem;
  arma_vec edges(nPts + 1);

  precision_t dc = centers(1) - centers(0);

  edges(0) = centers(0) - dc / 2.0;
  edges(1) = centers(0) + dc / 2.0;

  for (int64_t i = 2; i < nPts + 1; i++)
    edges(i) = 2 * centers(i - 1) - edges(i - 1);

  return edges;
}

// ---------------------------------------------------------
// bin edges
// ---------------------------------------------------------

arma_mat calc_bin_edges(arma_mat centers, bool DoX) {

  // X is first dimension (row), Y is second dimension (col)

  int64_t nX = centers.n_rows;
  int64_t nY = centers.n_cols;
  arma_mat edges;
  arma_vec centers1d;

  if (DoX) {
    if (verbose > 2)
      std::cout << "  --> x\n";

    edges.resize(nX + 1, nY);

    for (int64_t j = 0; j < nY; j++) {
      centers1d = centers.col(j);
      edges.col(j) = calc_bin_edges(centers1d);
    }
  } else {
    if (verbose > 2)
      std::cout << "  --> y\n";

    edges.resize(nX, nY + 1);

    for (int64_t i = 0; i < nX; i++) {
      centers1d = centers.row(i).as_col();
      edges.row(i) = calc_bin_edges(centers1d).as_row();
    }
  }

  return edges;
}

// ---------------------------------------------------------
// bin widths
// ---------------------------------------------------------

arma_vec calc_bin_widths(arma_vec edges) {

  int64_t nPts = edges.n_elem - 1;
  arma_vec widths(nPts);

  for (int64_t i = 0; i < nPts; i++)
    widths(i) = edges(i + 1) - edges(i);

  return widths;
}

// ---------------------------------------------------------
// bin widths 2d
// ---------------------------------------------------------

arma_mat calc_bin_widths(arma_mat edges, bool DoX) {

  int64_t nX = edges.n_rows;
  int64_t nY = edges.n_cols;

  arma_mat widths;
  arma_vec edges1d;

  if (DoX) {
    if (verbose > 2)
      std::cout << "  --> x\n";

    nX--;
    widths.resize(nX, nY);

    for (int64_t j = 0; j < nY; j++) {
      edges1d = edges.col(j);
      widths.col(j) = calc_bin_widths(edges1d);
    }
  } else {
    if (verbose > 2)
      std::cout << "  --> y\n";

    nY--;
    widths.resize(nX, nY);

    for (int64_t i = 0; i < nX; i++) {
      edges1d = edges.row(i).as_col();
      widths.row(i) = calc_bin_widths(edges1d).as_row();
    }
  }

  return widths;
}

/**SOME PROJECTION AND GRADIENT CODE **/
// ---------------------------------------------------------
//
// ---------------------------------------------------------

arma_vec limiter_mc(arma_vec &left,
                    arma_vec &right,
                    int64_t nPts,
                    int64_t nGCs) {

  precision_t beta = 0.8;

  arma_vec s = left % right;
  arma_vec combined = (left + right) * 0.5;

  left = left * beta;
  right = right * beta;
  arma_vec limited = left;

  for (int64_t i = 1; i < nPts + 2 * nGCs - 1; i++) {
    if (s(i) < 0) {
      // Sign < 0 means opposite signed left and right:
      limited(i) = 0.0;
    } else {
      if (left(i) > 0 && right(i) > 0) {
        if (right(i) < limited(i))
          limited(i) = right(i);

        if (combined(i) < limited(i))
          limited(i) = combined(i);
      } else {
        if (right(i) > limited(i))
          limited(i) = right(i);

        if (combined(i) > limited(i))
          limited(i) = combined(i);
      }
    }
  }

  return limited;
}

void print(arma_vec values) {
  int64_t nP = values.n_elem;

  for (int64_t i = 0; i < nP; i++)
    std::cout << values(i) << " ";

  std::cout << "\n";
}

// ---------------------------------------------------------
// calc gradients at centers
//   - values and x defined at centers
// ---------------------------------------------------------

arma_vec calc_grad_1d(arma_vec &values,
                      arma_vec &x,
                      int64_t nPts,
                      int64_t nGCs) {

  arma_vec gradients = values * 0.0;
  arma_vec gradL = values * 0.0;
  arma_vec gradR = values * 0.0;

  precision_t factor1 = 0.625;
  precision_t factor2 = 0.0416667;
  precision_t h;

  int64_t i;
  arma_vec hv = values * 0.0;

  i = nGCs - 1;
  h = 2.0 / (x(i + 1) - x(i));
  gradR(i) = h * (factor1 * (values(i + 1) - values(i)) -
                  factor2 * (values(i + 2) - values(i - 1)));
  gradL(i) = (values(i) - values(i - 1)) / (x(i) - x(i - 1));

  // This is attempting to vectorize the problem, but it seems to be slower?
  //  int64_t iS = nGCs;
  //  int64_t iE = nPts + nGCs - 1;
  //  hv.rows(iS, iE) = 2.0 / (x.rows(iS, iE) - x.rows(iS-1, iE-1));
  //  gradL.rows(iS, iE) = hv.rows(iS,iE) % (factor1 * (values.rows(iS, iE) -
  //                values.rows(iS-1, iE-1)) -
  //           factor2 * (values.rows(iS+1, iE+1) -
  //                values.rows(iS-2, iE-2)));
  //  hv.rows(iS, iE) = 2.0 / (x.rows(iS+1, iE+1) - x.rows(iS, iE));
  //  gradR.rows(iS, iE) = hv.rows(iS,iE) % (factor1 * (values.rows(iS+1, iE+1) -
  //                values.rows(iS, iE)) -
  //           factor2 * (values.rows(iS+2, iE+2) -
  //                values.rows(iS-1, iE-1)));

  for (i = nGCs; i < nPts + nGCs; i++) {
    h = 2.0 / (x(i) - x(i - 1));
    gradL(i) = h * (factor1 * (values(i) - values(i - 1)) -
                    factor2 * (values(i + 1) - values(i - 2)));
    h = 2.0 / (x(i + 1) - x(i));
    gradR(i) = h * (factor1 * (values(i + 1) - values(i)) -
                    factor2 * (values(i + 2) - values(i - 1)));
  }

  i = nPts + nGCs;
  h = 2.0 / (x(i) - x(i - 1));
  gradL(i) = h * (factor1 * (values(i) - values(i - 1)) -
                  factor2 * (values(i + 1) - values(i - 2)));
  gradR(i) = (values(i + 1) - values(i)) / (x(i + 1) - x(i));

  gradients = limiter_mc(gradL, gradR, nPts, nGCs);

  return gradients;
}

// ---------------------------------------------------------
// calc gradients at centers for 2d matrices
//   - values and x defined at centers
// ---------------------------------------------------------

arma_mat calc_grad(arma_mat values,
                   arma_mat x,
                   int64_t nGCs,
                   bool DoX) {

  arma_mat v2d, x2d;

  if (DoX) {
    v2d = values;
    x2d = x;
  } else {
    v2d = values.t();
    x2d = x.t();
  }

  int64_t nX = v2d.n_rows;
  int64_t nY = v2d.n_cols;
  arma_mat grad2d = v2d * 0.0;

  int64_t nPts = nX - 2 * nGCs;
  arma_vec values1d(nX);
  arma_vec x1d(nX);

  for (int64_t j = 1; j < nY - 1; j++) {
    values1d = v2d.col(j);
    x1d = x2d.col(j);
    grad2d.col(j) = calc_grad_1d(values1d, x1d, nPts, nGCs);
  }

  arma_mat gradients;

  if (DoX)
    gradients = grad2d;
  else
    gradients = grad2d.t();

  return gradients;
}

// ---------------------------------------------------------
// Project gradients + values to the right face, from the left
//   returned values are on the i - 1/2 edges
//     (between i-1 and i cell center)
// ---------------------------------------------------------

arma_mat project_from_left(arma_mat values,
                           arma_mat gradients,
                           arma_mat x_centers,
                           arma_mat x_edges,
                           int64_t nGCs) {

  int64_t nX = values.n_rows;
  int64_t nY = values.n_cols;

  // Define at edges:
  arma_mat projected(nX + 1, nY);
  projected.zeros();

  // no gradient in the 0 or iEnd cells
  for (int64_t j = 0; j < nY; j++) {
    for (int64_t i = 1; i < nX - 1; i++) {
      projected(i + 1, j) = values(i, j) +
                            gradients(i, j) * (x_edges(i + 1, j) - x_centers(i, j));
    }

    projected(1, j) = projected(2, j);
    projected(0, j) = projected(1, j);
    projected(nX, j) = projected(nX - 1, j);
  }

  return projected;
}

// ---------------------------------------------------------
// Project gradients + values to the left face, from the right
//   returned values are on the i - 1 edges
//     (between i-1 and i cell center)
// ---------------------------------------------------------

arma_mat project_from_right(arma_mat values,
                            arma_mat gradients,
                            arma_mat x_centers,
                            arma_mat x_edges,
                            int64_t nGCs) {
  int64_t nX = values.n_rows;
  int64_t nY = values.n_cols;

  // Define at edges:
  arma_mat projected(nX + 1, nY);
  projected.zeros();

  // no gradient in the 0 or iEnd cells
  for (int64_t j = 0; j < nY; j++) {
    for (int64_t i = 1; i < nX - 1; i++) {
      projected(i, j) = values(i, j) +
                        gradients(i, j) * (x_edges(i, j) - x_centers(i, j));
    }

    projected(0, j) = projected(1, j);
    projected(nX - 1, j) = projected(nX - 2, j);
    projected(nX, j) = projected(nX - 1, j);
  }

  return projected;
}

// ---------------------------------------------------------
// Limiter on values
//   projected is assumed to be on the edge between the
//   i-1 and i cell (i-1/2)
//   limited is returned at edges
// ---------------------------------------------------------

arma_vec limiter_value(arma_vec projected,
                       arma_vec values,
                       int64_t nPts,
                       int64_t nGCs) {

  int64_t iStart = 0;
  int64_t iEnd = nPts + 2 * nGCs;

  arma_vec limited = projected;

  precision_t mini, maxi;

  for (int64_t i = iStart + 1; i < iEnd - 1; i++) {

    mini = values(i - 1);

    if (values(i) < mini)
      mini = values(i);

    maxi = values(i - 1);

    if (values(i) > maxi)
      maxi = values(i);

    if (limited(i) < mini)
      limited(i) = mini;

    if (limited(i) > maxi)
      limited(i) = maxi;
  }

  return limited;
}

// ---------------------------------------------------------
// take gradients and project to all edges
// ---------------------------------------------------------

projection_struct project_to_edges(arma_mat &values,
                                   arma_mat &x_centers, arma_mat &x_edges,
                                   arma_mat &y_centers, arma_mat &y_edges,
                                   int64_t nGCs) {

  int64_t nX = values.n_rows;
  int64_t nY = values.n_cols;

  projection_struct proj;

  proj.gradLR = calc_grad(values, x_centers, nGCs, true);
  proj.gradDU = calc_grad(values.t(), y_centers.t(), nGCs, true).t();

  proj.R = project_from_left(values, proj.gradLR,
                             x_centers, x_edges, nGCs);
  // Left side of edge from left
  proj.L = project_from_right(values, proj.gradLR,
                              x_centers, x_edges, nGCs);
  // Up side of edge from down (left)
  proj.U = project_from_left(values.t(), proj.gradDU.t(),
                             y_centers.t(), y_edges.t(), nGCs)
           .t();
  // Down side of edge from up (right)
  proj.D = project_from_right(values.t(), proj.gradDU.t(),
                              y_centers.t(), y_edges.t(), nGCs)
           .t();

  return proj;
}

/**** SOME INITIALIZATION FUNCTION, NOT CORE ****/
// ---------------------------------------------------------
// initial rho: initialize the whole field to be 2.0
// ---------------------------------------------------------
arma_mat init_rho(arma_mat &x,
                  arma_mat &y) {

  int64_t nX = x.n_rows;
  int64_t nY = x.n_cols;

  arma_mat rho(nX, nY);
  arma_mat r;

  r = sqrt((x - 0.0) % (x - 0.0) + (y - 0.0) % (y - 0.0));
  rho.fill(1.0);
  rho.elem( find( r < 0.25)).fill(1.2);

  return rho;
}

// ---------------------------------------------------------
// initial velocity: initialize zero velocity
// ---------------------------------------------------------

arma_mat init_vel(arma_mat &x,
                  arma_mat &y) {
  int64_t nX = x.n_rows;
  int64_t nY = x.n_cols;
  arma_mat vel(nX, nY);
  // all cells positive to right:
  vel.zeros();
  vel.fill(0.5);
  return vel;
}


// ---------------------------------------------------------
// initial values
// ---------------------------------------------------------

arma_mat init_value(arma_mat &x,
                    arma_mat &y,
                    precision_t inVal) {
  int64_t nX = x.n_rows;
  int64_t nY = x.n_cols;
  arma_mat val(nX, nY);
  val.fill(inVal);
  return val;
}

arma_mat init_vel2(arma_mat &x,
                   arma_mat &y) {
  int64_t nX = x.n_rows;
  int64_t nY = x.n_cols;
  arma_mat vel(nX, nY);
  // all cells positive to right:
  vel.zeros();
  //vel.fill(1.0);
  return vel;
}

// ---------------------------------------------------------
// initial temp (E): constant total energy
// THIS IS NOT e but E, the total energy
// ---------------------------------------------------------

arma_mat init_temp(arma_mat &x,
                   arma_mat &y) {
  int64_t nX = x.n_rows;
  int64_t nY = x.n_cols;

  arma_mat temp(nX, nY);
  temp.fill(100.0);
  return temp;
}

// ---------------------------------------------------------
// Calculate the max speed in the x and y directions
// ---------------------------------------------------------

void calc_cmax(arma_mat &xVel,
               arma_mat &yVel,
               arma_mat &temp,
               arma_mat &xMax,
               arma_mat &yMax) {

  if (verbose > 2)
    std::cout << "  --> calc_max\n";

  arma_mat xVel2, yVel2;

  xVel2 = xVel % xVel;
  yVel2 = yVel % yVel;
  xMax = sqrt(xVel2) + sqrt(cKb / mmm * temp);
  yMax = sqrt(yVel2) + sqrt(cKb / mmm * temp);

  return;
}


// ---------------------------------------------------------
// Set Boundary Conditions
// ---------------------------------------------------------

void set_bcs(arma_mat &rho,
             arma_mat &xVel,
             arma_mat &yVel,
             arma_mat &temp,
             grid_struct gridC) {

  if (verbose > 2)
    std::cout << "  --> set_bcs\n";

  // ------------------------------------------------
  // Exchange messages (set BCs, really):
  for (int64_t i = gridC.iXfirst_; i < gridC.iXlast_; i++) {
    for (int64_t j = 0; j < gridC.nGCs; j++) {
      // bottom bc:
      rho(i, gridC.iYfirst_ - 1 - j) = rho(i, gridC.iYlast_ - 1 - j);
      // top bc:
      rho(i, gridC.iYlast_ + j) = rho(i, gridC.iYfirst_ + j);
    }
  }

  for (int64_t j = gridC.iYfirst_; j < gridC.iYlast_; j++) {
    for (int64_t i = 0; i < gridC.nGCs; i++) {
      // left bc:
      rho(gridC.iXfirst_ - 1 - i, j) = rho(gridC.iXlast_ - 1 - i, j);
      // right bc:
      rho(gridC.iXlast_ + i, j) = rho(gridC.iXfirst_ + i, j);
    }
  }

  return;
}

// ---------------------------------------------------------
// Convert vector from Alat, Alon to Axi, Anu
//   -> Using equation (7) of Ronchi et al:
// ---------------------------------------------------------

void convert_vector_ll_to_xn(arma_mat aLon,
                             arma_mat aLat,
                             arma_mat &aXi,
                             arma_mat &aNu,
                             grid_struct grid) {

  // Ronchi defines aPhi = aLon, aTheta = -aLat
  aXi = -grid.Axt % aLat + grid.Axp % aLon;
  aNu = -grid.Ant % aLat + grid.Anp % aLon;
  return;
}

// ---------------------------------------------------------
// Convert vector from sphere to cube
//   -> Using equation (32) of Nair:
// ---------------------------------------------------------

void convert_vector_sphere_to_cube(arma_mat u,
                                   arma_mat v,
                                   arma_mat &u1,
                                   arma_mat &u2,
                                   grid_struct grid) {

  // Ronchi defines aPhi = aLon, aTheta = -aLat
  u1 = grid.s2c11 % u + grid.s2c12 % v;
  u2 = grid.s2c21 % u + grid.s2c22 % v;
  return;
}

// ---------------------------------------------------------
// Convert vector cube to sphere
//   -> Using equation (32) of Nair:
// ---------------------------------------------------------

void convert_vector_cube_to_sphere(arma_mat u1,
                                   arma_mat u2,
                                   arma_mat &u,
                                   arma_mat &v,
                                   grid_struct grid) {

  // Ronchi defines aPhi = aLon, aTheta = -aLat
  u = grid.c2s11 % u1 + grid.c2s12 % u2;
  v = grid.c2s21 % u1 + grid.c2s22 % u2;
  return;
}

// ---------------------------------------------------------
// Convert vector from Alat, Alon to Axi, Anu
//   -> Using equation (7) of Ronchi et al:
// ---------------------------------------------------------

void convert_vector_xn_to_ll(arma_mat aXi,
                             arma_mat aNu,
                             arma_mat &aLon,
                             arma_mat &aLat,
                             grid_struct grid) {

  // Ronchi defines aPhi = aLon, aTheta = -aLat
  aLat = -(grid.Atx % aXi + grid.Atn % aNu);
  aLon = grid.Apx % aXi + grid.Apn % aNu;

  return;
}

// ---------------------------------------------------------
// Convert vector from Alat, Alon to Axi, Anu
//   -> Using equation (7) of Ronchi et al:
// ---------------------------------------------------------

arma_mat calc_angle_between_coords(grid_struct grid) {

  arma_mat e1Lat, e1Lon, e2Lat, e2Lon, m, one, zero;
  arma_mat e1dote2, alpha;

  m.resize(grid.nXt, grid.nYt);
  one.resize(grid.nXt, grid.nYt);
  one.fill(1.0);
  zero.resize(grid.nXt, grid.nYt);
  zero.fill(0.0);

  // define e1 as the LR (xi) direction:
  e1Lat.resize(grid.nXt, grid.nYt);
  e1Lon.resize(grid.nXt, grid.nYt);
  convert_vector_xn_to_ll(one, zero, e1Lon, e1Lat, grid);
  m = sqrt(e1Lon % e1Lon + e1Lat % e1Lat);
  e1Lon = e1Lon / m;
  e1Lat = e1Lat / m;

  // define e2 as the DU (nu) direction:
  e2Lat.resize(grid.nXt, grid.nYt);
  e2Lon.resize(grid.nXt, grid.nYt);
  convert_vector_xn_to_ll(zero, one, e2Lon, e2Lat, grid);
  m = sqrt(e2Lon % e2Lon + e2Lat % e2Lat);
  e2Lon = e2Lon / m;
  e2Lat = e2Lat / m;

  alpha = acos(e1Lat % e2Lat + e1Lon % e2Lon);

  return alpha;
}

// ---------------------------------------------------------
// Convert vector from Alat, Alon to Axi, Anu
//   -> Using equation (7) of Ronchi et al:
// ---------------------------------------------------------

void calc_norms(grid_struct &grid) {

  arma_mat e1Lat, e1Lon, e2Lat, e2Lon, m, one, zero;

  grid.nXiLon.resize(grid.nXt, grid.nYt);
  grid.nXiLat.resize(grid.nXt, grid.nYt);
  grid.nNuLon.resize(grid.nXt, grid.nYt);
  grid.nNuLat.resize(grid.nXt, grid.nYt);

  m.resize(grid.nXt, grid.nYt);
  one.resize(grid.nXt, grid.nYt);
  one.fill(1.0);
  zero.resize(grid.nXt, grid.nYt);
  zero.fill(0.0);

  // define e1 as the LR (xi) direction:
  e1Lat.resize(grid.nXt, grid.nYt);
  e1Lon.resize(grid.nXt, grid.nYt);
  convert_vector_xn_to_ll(one, zero, e1Lon, e1Lat, grid);
  m = sqrt(e1Lon % e1Lon + e1Lat % e1Lat);

  // Rotate by 90 deg (CCW) to get the norm:
  grid.nNuLon = -e1Lat / m;
  grid.nNuLat = e1Lon / m;

  // define e2 as the DU (nu) direction:
  e2Lat.resize(grid.nXt, grid.nYt);
  e2Lon.resize(grid.nXt, grid.nYt);
  convert_vector_xn_to_ll(zero, one, e2Lon, e2Lat, grid);
  m = sqrt(e2Lon % e2Lon + e2Lat % e2Lat);
  // Rotate by 90 deg (CW) to get the norm:
  grid.nXiLon = e2Lat / m;
  grid.nXiLat = -e2Lon / m;

  return;
}


// ---------------------------------------------------------
// Update States
// ---------------------------------------------------------

void update_states(arma_mat rho,
                   arma_mat &xVel,
                   arma_mat &yVel,
                   arma_mat &temp,
                   arma_mat &drhodt,
                   arma_mat &dlonVeldt,
                   arma_mat &dlatVeldt,
                   arma_mat &dtempdt,
                   grid_struct gridC,
                   grid_struct gridL,
                   grid_struct gridD,
                   precision_t dt) {

  arma_mat xMomentum, yMomentum;
  arma_mat rhoE, energy, vel2;

  precision_t cv = 1500.0;

  if (verbose > 2)
    std::cout << "  --> update_states\n";

  // Derived variables:
  xMomentum = rho % xVel; // x1momentum, pure scalar field
  yMomentum = rho % yVel; // y1momentum, pure scalar field
  rhoE = rho % temp;

  vel2 = xVel % xVel + yVel % yVel;
  //energy = rho % (0.5 * vel2 + cv * temp);
  energy = cv * rho % temp;

  /** Initialize projection constructs */
  static projection_struct rhoP;
  static projection_struct xMomentumP, xVelP;
  static projection_struct yMomentumP, yVelP;
  static projection_struct energyP;
  static projection_struct tempP;

  // They are all pure scalar fields without sqrt(g)
  static arma_mat totaleL, totaleR, totaleD, totaleU;
  static arma_mat velL2, velR2, velD2, velU2;
  static arma_mat pressureL, pressureR, pressureD, pressureU;

  arma_mat dxVeldt = xVel * 0.0;
  arma_mat dyVeldt = yVel * 0.0;

  dlonVeldt = dxVeldt * 0.0 + 1;
  dlatVeldt = dyVeldt * 0.0 + 1;

  static arma_mat velNormL, velNormR, velNormU, velNormD;

  /** Initialize Flux and Wave Speed Storages */
  static arma_mat eq1FluxLR, eq1FluxDU;
  static arma_mat eq1FluxL, eq1FluxR, eq1FluxD, eq1FluxU;
  static arma_mat eq2FluxLR, eq2FluxDU;
  static arma_mat eq2FluxL, eq2FluxR, eq2FluxD, eq2FluxU;
  static arma_mat eq3FluxLR, eq3FluxDU;
  static arma_mat eq3FluxL, eq3FluxR, eq3FluxD, eq3FluxU;
  static arma_mat eq4FluxLR, eq4FluxDU;
  static arma_mat eq4FluxL, eq4FluxR, eq4FluxD, eq4FluxU;

  arma_mat wsL, wsR, wsD, wsU, wsLR, wsDU;

  arma_mat diff; // for Riemann Solver

  if (verbose > 3)
    std::cout << "  ---> Projecting\n";

  rhoP = project_to_edges(rho, gridC.xi, gridL.xi, gridC.nu, gridD.nu,
                          gridC.nGCs);
  // project the lon / lat velocities to the edges:
  xVelP = project_to_edges(xVel, gridC.xi, gridL.xi, gridC.nu, gridD.nu,
                           gridC.nGCs);
  yVelP = project_to_edges(yVel, gridC.xi, gridL.xi, gridC.nu, gridD.nu,
                           gridC.nGCs);
  xMomentumP = project_to_edges(xMomentum, gridC.xi, gridL.xi, gridC.nu, gridD.nu,
                                gridC.nGCs);
  yMomentumP = project_to_edges(yMomentum, gridC.xi, gridL.xi, gridC.nu, gridD.nu,
                                gridC.nGCs);
  energyP = project_to_edges(energy, gridC.xi, gridL.xi, gridC.nu, gridD.nu,
                             gridC.nGCs);
  tempP = project_to_edges(temp, gridC.xi, gridL.xi, gridC.nu, gridD.nu,
                           gridC.nGCs);

  if (verbose > 3)
    std::cout << "  ---> Derived values\n";

  velL2 = (xVelP.L % xVelP.L + yVelP.L % yVelP.L);
  velR2 = (xVelP.R % xVelP.R + yVelP.R % yVelP.R);
  velD2 = (xVelP.D % xVelP.D + yVelP.D % yVelP.D);
  velU2 = (xVelP.U % xVelP.U + yVelP.U % yVelP.U);

  precision_t k = 1.38e-23;
  // let's be Oxygen:
  precision_t mass = 16.0 * 1.67e-27;
  pressureL = k / mass * (rhoP.L % tempP.L);
  pressureR = k / mass * (rhoP.R % tempP.R);
  pressureD = k / mass * (rhoP.D % tempP.D);
  pressureU = k / mass * (rhoP.U % tempP.U);

  arma_mat pressureLR = (pressureL + pressureR) / 2;
  arma_mat pressureDU = (pressureD + pressureU) / 2;

  if (verbose > 3)
    std::cout << "  ---> Normal Velocities\n";

  // Calculate the normal velocity at the boundaries:
  velNormL = xVelP.L % gridL.nXiLon + yVelP.L % gridL.nXiLat;
  velNormR = xVelP.R % gridL.nXiLon + yVelP.R % gridL.nXiLat;
  velNormU = xVelP.U % gridD.nNuLon + yVelP.U % gridD.nNuLat;
  velNormD = xVelP.D % gridD.nNuLon + yVelP.D % gridD.nNuLat;

  if (verbose > 3)
    std::cout << "  ---> Fluxes eq 1\n";

  // Flux calculated from the left of the edge
  eq1FluxL = rhoP.L % velNormL;
  // Flux calculated from the right of the edge
  eq1FluxR = rhoP.R % velNormR;
  // Flux calculated from the down of the edge
  eq1FluxD = rhoP.D % velNormD;
  // Flux calculated from the up of the edge
  eq1FluxU = rhoP.U % velNormU;

  if (verbose > 3)
    std::cout << "  ---> Fluxes eq 2\n";

  eq2FluxL = (xMomentumP.L % velNormL);
  eq2FluxR = (xMomentumP.R % velNormR);
  eq2FluxD = (xMomentumP.D % velNormD);
  eq2FluxU = (xMomentumP.U % velNormU);

  if (verbose > 3)
    std::cout << "  ---> Fluxes eq 3\n";

  eq3FluxL = (yMomentumP.L % velNormL);
  eq3FluxR = (yMomentumP.R % velNormR);
  eq3FluxD = (yMomentumP.D % velNormD);
  eq3FluxU = (yMomentumP.U % velNormU);

  eq4FluxL = energyP.L % velNormL;
  eq4FluxR = energyP.R % velNormR;
  eq4FluxD = energyP.D % velNormD;
  eq4FluxU = energyP.U % velNormU;

  // ------------------------------------------------
  // Calculate the wave speed for the diffusive flux:
  // In Reference velocities
  if (verbose > 3)
    std::cout << "  ---> Diffusive Fluxes\n";

  wsL = sqrt(velL2) + sqrt(cGamma * (cGamma - 1) * tempP.L);
  wsR = sqrt(velR2) + sqrt(cGamma * (cGamma - 1) * tempP.R);
  wsD = sqrt(velD2) + sqrt(cGamma * (cGamma - 1) * tempP.D);
  wsU = sqrt(velU2) + sqrt(cGamma * (cGamma - 1) * tempP.U);

  wsLR = wsR;

  for (int64_t i = 0; i < gridC.nXt + 1; i++) {
    for (int64_t j = 0; j < gridC.nYt; j++) {
      if (wsL(i, j) > wsLR(i, j))
        wsLR(i, j) = wsL(i, j);
    }
  }

  wsDU = wsD;

  for (int64_t i = 0; i < gridC.nXt; i++) {
    for (int64_t j = 0; j < gridC.nYt + 1; j++) {
      if (wsU(i, j) > wsDU(i, j))
        wsDU(i, j) = wsU(i, j);
    }
  }

  // ------------------------------------------------
  // Calculate average flux at the edges (Rusanov Flux):

  if (verbose > 3)
    std::cout << "  ---> Averaging fluxes at edges\n";

  diff = (rhoP.R - rhoP.L);
  eq1FluxLR = (eq1FluxL + eq1FluxR) / 2 + 0.5 * wsLR % diff;
  diff = (rhoP.U - rhoP.D);
  eq1FluxDU = (eq1FluxD + eq1FluxU) / 2 + 0.5 * wsDU % diff;

  diff = (xMomentumP.R - xMomentumP.L);
  eq2FluxLR = (eq2FluxL + eq2FluxR) / 2 + 0.5 * wsLR % diff;
  diff = (xMomentumP.U - xMomentumP.D);
  eq2FluxDU = (eq2FluxD + eq2FluxU) / 2 + 0.5 * wsDU % diff;

  diff = (yMomentumP.R - yMomentumP.L);
  eq3FluxLR = (eq3FluxL + eq3FluxR) / 2 + 0.5 * wsLR % diff;
  diff = (yMomentumP.U - yMomentumP.D);
  eq3FluxDU = (eq3FluxD + eq3FluxU) / 2 + 0.5 * wsDU % diff;

  diff = (energyP.R - energyP.L);
  eq4FluxLR = (eq4FluxL + eq4FluxR) / 2 + 0.5 * wsLR % diff;
  diff = (energyP.U - energyP.D);
  eq4FluxDU = (eq4FluxD + eq4FluxU) / 2 + 0.5 * wsDU % diff;

  // ------------------------------------------------
  // Update values:
  if (verbose > 3)
    std::cout << "  ---> Updating equations of state\n";

  precision_t dpdx, dpdn, pp, pm;

  arma_mat ax, an;

  ax = xVel * 0.0;
  an = yVel * 0.0;
  arma_mat dedt = xVel * 0.0;

  arma_mat rhoNew = rho;

  // Only deal with inner cell
  for (int64_t j = gridC.iYfirst_; j < gridC.iYlast_; j++) {
    for (int64_t i = gridC.iXfirst_; i < gridC.iXlast_; i++) {
      precision_t rhoResidual_ij = (gridL.dln(i + 1, j) * eq1FluxLR(i + 1, j) -
                                    gridL.dln(i, j) * eq1FluxLR(i, j) +
                                    gridD.dlx(i, j + 1) * eq1FluxDU(i, j + 1) -
                                    gridD.dlx(i, j) * eq1FluxDU(i, j));
      drhodt(i, j) = rhoResidual_ij / gridC.dS(i, j);

      rhoNew(i, j) = rho(i, j) + dt * drhodt(i, j);

      precision_t xMomentumResidual_ij = (gridL.dln(i + 1, j) * eq2FluxLR(i + 1, j) -
                                          gridL.dln(i, j) * eq2FluxLR(i, j) +
                                          gridD.dlx(i, j + 1) * eq2FluxDU(i, j + 1) -
                                          gridD.dlx(i, j) * eq2FluxDU(i, j));
      dxVeldt(i, j) = xMomentumResidual_ij / gridC.dS(i, j) / rhoNew(i, j);

      precision_t yMomentumResidual_ij = (gridL.dln(i + 1, j) * eq3FluxLR(i + 1, j) -
                                          gridL.dln(i, j) * eq3FluxLR(i, j) +
                                          gridD.dlx(i, j + 1) * eq3FluxDU(i, j + 1) -
                                          gridD.dlx(i, j) * eq3FluxDU(i, j));
      dyVeldt(i, j) = yMomentumResidual_ij / gridC.dS(i, j) / rhoNew(i, j);

      // Calculate the gradient in the potential in the cubesphere
      // coordinate system:
      dpdx = 1 / gridC.R * gridC.D(i, j) *
             (pressureLR(i + 1, j) - pressureLR(i, j)) / gridC.dxi;
      dpdn = 1 / gridC.R * gridC.X(i, j) * gridC.Y(i, j) /
             gridC.D(i, j) *
             (pressureDU(i, j + 1) - pressureDU(i, j)) / gridC.dnu;
      ax(i, j) = (dpdx + dpdn) / rhoNew(i, j);

      dpdx = 1 / gridC.R * gridC.X(i, j) * gridC.Y(i, j) /
             gridC.C(i, j) * (pressureLR(i + 1, j) - pressureLR(i, j)) / gridC.dxi;
      dpdn = 1 / gridC.R * gridC.C(i, j) *
             (pressureDU(i, j + 1) - pressureDU(i, j)) / gridC.dnu;
      an(i, j) = (dpdx + dpdn) / rhoNew(i, j);

      precision_t energyResidual_ij = (gridL.dln(i + 1, j) * eq4FluxLR(i + 1, j) -
                                       gridL.dln(i, j) * eq4FluxLR(i, j) +
                                       gridD.dlx(i, j + 1) * eq4FluxDU(i, j + 1) -
                                       gridD.dlx(i, j) * eq4FluxDU(i, j));
      dedt(i, j) = energyResidual_ij / gridC.dS(i, j);

    }
  }

  dlatVeldt = dyVeldt - (ax % gridC.Atx + an % gridC.Atn);
  dlonVeldt = dxVeldt + ax % gridC.Apx + an % gridC.Apn;
  dtempdt = dedt / rhoNew / cv;

  return;
}


// ---------------------------------------------------------
// Main Code!
// ---------------------------------------------------------

int main() {
  precision_t dt;        // Time Step
  precision_t current_time = 0.0; // Initial Time 0
  precision_t total_time = 100.0;   // Total simulation time
  precision_t cfl = 0.75;
  int iFace = 5;

  precision_t dtOut = total_time / 50.0; // Output Interval
  precision_t dtReport = total_time / 200.0; // Output Interval

  int64_t iStep;    // Iterator of Time Step

  int64_t nX = 50;  // Number of x grid cells
  int64_t nY = 50; // Number of y grid cells
  int64_t nGCs = 2; // Number of ghost cells

  // Radius of Sphere
  precision_t R = 10000.;

  if (verbose > 0)
    std::cout << "> generating cubesphere cell center and metrics\n";

  grid_struct gridC_eq = init_grid_equidistant(iFace, nX, nY, nGCs, R, 0.5, 0.5);

  output(gridC_eq.lat, "eq_lat.txt", false);
  output(gridC_eq.lon, "eq_lon.txt", false);
  output(gridC_eq.r, "eq_r.txt", false);

  grid_struct gridC = init_grid(iFace, nX, nY, nGCs, R, 0.5, 0.5);
  grid_struct gridL = init_grid(iFace, nX + 1, nY, nGCs, R, 0.0, 0.5);
  grid_struct gridD = init_grid(iFace, nX, nY + 1, nGCs, R, 0.5, 0.0);

  gridC.alpha = calc_angle_between_coords(gridC);
  gridL.alpha = calc_angle_between_coords(gridL);
  gridD.alpha = calc_angle_between_coords(gridD);

  gridC.sinAlpha = sin(gridC.alpha);
  gridL.sinAlpha = sin(gridL.alpha);
  gridD.sinAlpha = sin(gridD.alpha);

  calc_norms(gridC);
  calc_norms(gridL);
  calc_norms(gridD);

  /** State Initialization */
  /// Initialize Density
  if (verbose > 0)
    std::cout << "> initializing rho\n";

  arma_mat rho = init_rho(gridC.xi, gridC.nu); // rho, pure scalar field

  /// Initialize Velocity and Momentum
  if (verbose > 0)
    std::cout << "> initializing vel\n";

  // Initialize spherical velocity
  // Supposed to be Longitudinal Velocity:
  arma_mat vLon = init_value(gridC.xi, gridC.nu, 0.0);
  // Supposed to be Latitudinal Velocity:
  arma_mat vLat = init_value(gridC.xi, gridC.nu, 0.0);

  arma_mat vXi, vNu;
  convert_vector_ll_to_xn(vLon, vLat, vXi, vNu, gridC);


  /// Initialize total energy
  if (verbose > 0)
    std::cout << "> initializing energy\n";

  arma_mat temp = init_temp(gridC.xi, gridC.nu);

  /** Output some pre-simulation results */
  if (verbose > 0)
    std::cout << "-> outputting\n";

  output(gridC.xi, "xi.txt", false);
  output(gridC.dS, "dS.txt", false);
  output(gridC.d, "d.txt", false);
  output(gridC.nu, "nu.txt", false);
  output(gridC.lat, "lat.txt", false);
  output(gridC.lon, "lon.txt", false);
  output(gridC.alpha, "alpha.txt", false);
  output(rho, "rho.txt", false);
  output(vLon, "vLon.txt", false);
  output(vLat, "vLat.txt", false);
  output(vXi, "vXi.txt", false);
  output(vNu, "vNu.txt", false);
  output(temp, "temp.txt", false);
  iStep = 0;

  arma_mat xMax, yMax;

  arma_mat drhodt;
  arma_mat dxVeldt, dyVeldt;
  arma_mat dtempdt;

  drhodt.resize(gridC.nXt, gridC.nYt);
  dxVeldt.resize(gridC.nXt, gridC.nYt);
  dyVeldt.resize(gridC.nXt, gridC.nYt);
  dtempdt.resize(gridC.nXt, gridC.nYt);

  arma_mat k1rho, k2rho, k3rho, k4rho, rhoInterK1, rhoInterK2, rhoInterK3;
  k1rho.resize(gridC.nXt, gridC.nYt);
  k2rho.resize(gridC.nXt, gridC.nYt);
  k3rho.resize(gridC.nXt, gridC.nYt);
  k4rho.resize(gridC.nXt, gridC.nYt);
  rhoInterK1.resize(gridC.nXt, gridC.nYt);
  rhoInterK2.resize(gridC.nXt, gridC.nYt);
  rhoInterK3.resize(gridC.nXt, gridC.nYt);

  arma_mat k1vLon, k2vLon, k3vLon, k4vLon, vLonInterK1, vLonInterK2, vLonInterK3;
  k1vLon.resize(gridC.nXt, gridC.nYt);
  k2vLon.resize(gridC.nXt, gridC.nYt);
  k3vLon.resize(gridC.nXt, gridC.nYt);
  k4vLon.resize(gridC.nXt, gridC.nYt);
  vLonInterK1.resize(gridC.nXt, gridC.nYt);
  vLonInterK2.resize(gridC.nXt, gridC.nYt);
  vLonInterK3.resize(gridC.nXt, gridC.nYt);

  arma_mat k1vLat, k2vLat, k3vLat, k4vLat, vLatInterK1, vLatInterK2, vLatInterK3;
  k1vLat.resize(gridC.nXt, gridC.nYt);
  k2vLat.resize(gridC.nXt, gridC.nYt);
  k3vLat.resize(gridC.nXt, gridC.nYt);
  k4vLat.resize(gridC.nXt, gridC.nYt);
  vLatInterK1.resize(gridC.nXt, gridC.nYt);
  vLatInterK2.resize(gridC.nXt, gridC.nYt);
  vLatInterK3.resize(gridC.nXt, gridC.nYt);

  arma_mat k1temp, k2temp, k3temp, k4temp, tempInterK1, tempInterK2, tempInterK3;
  k1temp.resize(gridC.nXt, gridC.nYt);
  k2temp.resize(gridC.nXt, gridC.nYt);
  k3temp.resize(gridC.nXt, gridC.nYt);
  k4temp.resize(gridC.nXt, gridC.nYt);
  tempInterK1.resize(gridC.nXt, gridC.nYt);
  tempInterK2.resize(gridC.nXt, gridC.nYt);
  tempInterK3.resize(gridC.nXt, gridC.nYt);

  while (current_time < total_time) {

    if (int((current_time - dt) / dtReport) != int((current_time ) / dtReport)) {
      std::cout << "step : " << iStep << "; time : " << current_time << "\n";
      arma_vec amin, amax;
      precision_t mini, maxi;
      amin = rho.min();
      mini = amin.min();
      amax = max(rho, 1);
      maxi = arma::max(amax);
      std::cout << "  -> min/max (rho) : " << mini << " / " << maxi << "\n";
      amin = temp.min();
      mini = amin.min();
      amax = max(temp, 1);
      maxi = arma::max(amax);
      std::cout << "  -> min/max (temp) : " << mini << " / " << maxi << "\n";
    }

    calc_cmax(vLon, vLat, temp, xMax, yMax);
    dt = calc_dt(gridC.dlx, gridC.dln, xMax, yMax, gridC.nGCs);
    dt = cfl * dt;
    std::cout << " dt: " << dt << "\n";

    // k1 - start at t0, go to t+1/2 to figure out slope at t0 (k1)
    update_states(rho, vLon, vLat, temp,
                  k1rho, k1vLon, k1vLat, k1temp,
                  gridC, gridL, gridD, dt / 2);
    // Take 1/2 step to figure out t+1/2 values, using k1:
    rhoInterK1 = rho + k1rho * dt / 2;
    vLonInterK1 = vLon + k1vLon * dt / 2;
    vLatInterK1 = vLat + k1vLat * dt / 2;
    tempInterK1 = temp + k1temp * dt / 2;
    set_bcs(rhoInterK1, vLonInterK1, vLatInterK1, tempInterK1, gridC);

    // k2 - start at t+1/2, go to t+1 to figure out slope at t+1/2 (k2)
    update_states(rhoInterK1, vLonInterK1, vLatInterK1, tempInterK1,
                  k2rho, k2vLon, k2vLat, k2temp,
                  gridC, gridL, gridD, dt / 2);
    // Take 1/2 step to figure out t+1/2 values, using k2:
    rhoInterK2 = rho + k2rho * dt / 2;
    vLonInterK2 = vLon + k2vLon * dt / 2;
    vLatInterK2 = vLat + k2vLat * dt / 2;
    tempInterK2 = temp + k2temp * dt / 2;
    set_bcs(rhoInterK2, vLonInterK2, vLatInterK2, tempInterK2, gridC);

    // k3 - start at t+1/2, go to t+1 to figure out slope at t+1/2 (k3)
    update_states(rhoInterK2, vLonInterK2, vLatInterK2, tempInterK2,
                  k3rho, k3vLon, k3vLat, k3temp,
                  gridC, gridL, gridD, dt / 2);
    // Take full step to figure out k4, using k3 slope:
    rhoInterK3 = rho + k3rho * dt;
    vLonInterK3 = vLon + k3vLon * dt;
    vLatInterK3 = vLat + k3vLat * dt;
    tempInterK3 = temp + k3rho * dt;
    set_bcs(rhoInterK3, vLonInterK3, vLatInterK3, tempInterK3, gridC);

    // k4 - start at t+1, go to t+2 to figure out slope at t+1 (k4)
    update_states(rhoInterK3, vLonInterK3, vLatInterK3, tempInterK3,
                  k4rho, k4vLon, k4vLat, k4temp,
                  gridC, gridL, gridD, dt);

    rho = rho - dt / 6 * (k1rho + 2 * k2rho + 2 * k3rho + k4rho);
    vLon = vLon - dt / 6 * (k1vLon + 2 * k2vLon + 2 * k3vLon + k4vLon);
    vLat = vLat - dt / 6 * (k1vLat + 2 * k2vLat + 2 * k3vLat + k4vLat);
    temp = temp - dt / 6 * (k1temp + 2 * k2temp + 2 * k3temp + k4temp);
    set_bcs(rho, vLon, vLat, temp, gridC);

    iStep++;
    current_time += dt;

    if (verbose > 3)
      std::cout << "  ---> Outputing\n";

    if (int((current_time - dt) / dtOut) != int((current_time ) / dtOut)) {
      std::cout << "> Outputing at time : " << current_time << "\n";
      output(rho, "rho.txt", true);
      output(vLon, "vLon.txt", true);
      output(vLat, "vLat.txt", true);
      output(temp, "temp.txt", true);
    }
  }

  return 0;
}
