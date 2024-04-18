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
precision_t cPIdiv2 = cPI / 2;

// ---------------------------------------------------------
// A couple of global variables
// ---------------------------------------------------------

int64_t verbose = 1;

struct projection_struct
{

  arma_mat gradLR;
  arma_mat gradDU;
  arma_mat R;
  arma_mat L;
  arma_mat U;
  arma_mat D;
};

// ---------------------------------------------------------
//
// ---------------------------------------------------------

precision_t calc_dt(double dx,
                    double dy,
                    arma_mat &wsLR,
                    arma_mat &wsDU,
                    int64_t nGCs)
{

  int64_t nX = wsLR.n_rows;
  int64_t nY = wsLR.n_cols;

  precision_t wsX, wsY, dtX, dtY, dt;

  dt = 1e32;

  for (int64_t j = nGCs; j < nY - nGCs; j++)
  {
    for (int64_t i = nGCs; i < nX - nGCs; i++)
    {
      wsX = (wsLR(i + 1, j) + wsLR(i, j)) / 2;
      dtX = dx / wsX;
      wsY = (wsDU(i, j + 1) + wsDU(i, j)) / 2;
      dtY = dy / wsY;
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
            bool DoAppend)
{

  std::ofstream outfile;
  if (DoAppend)
    outfile.open(filename, std::ios_base::app);
  else
  {
    outfile.open(filename);
    int64_t nX = values.n_rows;
    int64_t nY = values.n_cols;
    outfile << nX << " " << nY << "\n";
  }
  outfile << values;
  outfile << "----";
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
                  precision_t r)
{
  arma_vec xyz(3);
  xyz(0) = r * std::cos(lat) * std::cos(lon);
  xyz(1) = r * std::cos(lat) * std::sin(lon);
  xyz(2) = r * std::sin(lat);
  return xyz;
}

/**
 * Edge metrics generation
 *
 * Params:
 * @param x x coordinates
 * @param y y coordinates
 * @param R Radius of Earth
 * @param DoX Do Xedge or Yedge
 *
 * "Return" in reference:
 * @param edge_jacobian Edge Jacobian
 * @param edge_g11_inv
 * @param edge_g12_inv
 * @param edge_g21_inv
 * @param edge_g22_inv
 */

void edge_metrics(int64_t face, arma_mat x, arma_mat y, precision_t R,
                  arma_mat &edge_jacobian, arma_mat &edge_g11_inv, arma_mat &edge_g12_inv,
                  arma_mat &edge_g21_inv, arma_mat &edge_g22_inv)
{

  double a = std::sqrt(3) / 3 * R;

  double x_cart, y_cart, z_cart;
  double xp, yp, zp, rp, latp, lonp, p;

  int64_t nX = edge_jacobian.n_rows;
  int64_t nY = edge_jacobian.n_cols;
  for (int64_t j = 0; j < nY; j++)
  {
    for (int64_t i = 0; i < nX; i++)
    {
      // Extract Coordinates
      xp = x(i, 0);
      yp = y(0, j);
      rp = std::sqrt(xp * xp + yp * yp + a * a);

      // Transformation from 2D reference to LatLong
      lonp = std::atan(xp / a) + 3 * cPI / 2.;
      latp = std::atan(yp / a * std::cos(lonp - 3 * cPI / 2.));

      // jacobian = R^2 * a / r^3
      double g = (R * R * a / (rp * rp * rp));
      edge_jacobian(i, j) = std::sqrt(g);

      // Lower indices, Eqn(B1), Nair
      p = R * R / (rp * rp * rp * rp);
      double g11 = p * (a * a + yp * yp);
      double g12 = -p * xp * yp;
      double g21 = -p * xp * yp;
      double g22 = p * (a * a + xp * xp);

      // Now upper indicies, Eqn(B1), Nair
      edge_g11_inv(i, j) = g22 / g;
      edge_g12_inv(i, j) = -g12 / g;
      edge_g21_inv(i, j) = -g21 / g;
      edge_g22_inv(i, j) = g11 / g;
    }
  }
}

/**
 * CubeSphere Grid Generation
 * Applies to both corner and edge generation as the mapping 
 * is the same.
 *
 * Below are input parameters
 * @param face Face Number
 * @param nX Number of cells along local x direction.
 * @param nY Number of cells along local y direction.
 * @param R Radius of Earth (of global spherical coordinates)
 * @param left_off, down_off offset value to define edge or cell center
 * @param lon2d, lat2d global longitude and latitude coordinates
 * @param DoMetric computes metrics or not
 *
 * Below are "return" values
 * @param jacobian Jacobian of the grid transformation
 * @param A11, A12, A21, A22, A11_inv, A12_inv, A21_inv, A22_inv
 * are transformation matrix
 * @param x, y returned local x and y coordinates
 *
 */
void generate_x_y_cubesphere(int64_t face, int64_t nX, int64_t nY, int64_t nGCs,
                             precision_t R,
                             precision_t left_off,
                             precision_t down_off,
                             arma_mat &lon2d, arma_mat &lat2d,
                             arma_mat &jacobian,
                             arma_mat &A11_inv, arma_mat &A12_inv,
                             arma_mat &A21_inv, arma_mat &A22_inv,
                             arma_mat &A11, arma_mat &A12,
                             arma_mat &A21, arma_mat &A22,
                             arma_mat &x, arma_mat &y, bool DoMetric)
{

  // Normalized global coordinates
  // (x_cart, y_cart, z_cart) in 3D Cartesian with spherical coord R
  double x_cart, y_cart, z_cart;
  // Coordinates at each desired point
  // (xp, yp, zp) are local coordinates
  // rp is a metric for xp yp zp
  // latp and longp are latitude and longitude for cell center
  double xp, yp, zp, rp, latp, lonp, p, p2;
  double dr, du;

  dr = 2.0 / (nX - 2.0 * nGCs);
  du = 2.0 / (nY - 2.0 * nGCs);

  // Range factor (Local face of CubeSphere ranges [-a, a])
  // doi: 10.1016/j.jcp.2007.07.022
  // a = sqrt(3)/3*R
  double a = std::sqrt(3) / 3 * R;

  // Loop through each point and derive the coordinate
  for (int iDU = 0; iDU < nY; iDU++)
  {
    for (int iLR = 0; iLR < nX; iLR++)
    {
      // the offsets are so we can find cell centers, edges, and corners
      double iD = iDU - nGCs + down_off;
      double iL = iLR - nGCs + left_off;

      /** Equidistant projection
       *
       * doi: 10.1016/j.jcp.2007.07.022
       * Section 3.1
       *
       * (X, Y, Z) = R/r(x, -a, y)
       * r = sqrt(a^2+x^2+y^2)
       *
       */

      // Define local coordinates based on a
      xp = (-1.0 + dr * iL) * a;
      yp = (-1.0 + du * iD) * a;
      rp = std::sqrt(xp * xp + yp * yp + a * a);

      // Normalize local variable
      // In fact, transforming local variable to global Cartesian
      // x_cart = xp / rp * R;
      // y_cart = -a / rp * R;
      // z_cart = yp / rp * R;

      // Transformation from 3D Cartesian to LatLong
      // lonp = std::atan2(y_cart, x_cart) + cPI/2.0;
      lonp = std::atan(xp / a) + 3 * cPI / 2.;
      latp = std::atan(yp / a * std::cos(lonp - 3 * cPI / 2.));

      // Fill Computed coords
      lat2d(iLR, iDU) = latp;
      lon2d(iLR, iDU) = lonp;
      x(iLR, iDU) = xp;
      y(iLR, iDU) = yp;

      if (DoMetric)
      {
        /** Jacobian and Metric Tensor
         *
         * doi.org/10.1175/MWR2890.1
         * Section 4b
         *
         */

        // jacobian = R^2 * a / r^3
        jacobian(iLR, iDU) = R * R * a / (rp * rp * rp);

        // Now A inv, Eqn(B10), Nair
        p = a / std::cos(latp) / std::cos(lonp- 3 * cPI / 2.) / R;
        A11_inv(iLR, iDU) = p / std::cos(lonp- 3 * cPI / 2.);
        A12_inv(iLR, iDU) = 0;
        A21_inv(iLR, iDU) = p * std::tan(latp) * std::tan(lonp- 3 * cPI / 2.);
        A22_inv(iLR, iDU) = p / std::cos(latp);

        // Just A
        p2 = R * std::cos(latp) * std::cos(lonp- 3 * cPI / 2.) / a;
        A11(iLR, iDU) = p2 * std::cos(lonp- 3 * cPI / 2.);
        A12(iLR, iDU) = 0;
        A21(iLR, iDU) = -p2 * std::sin(latp) * std::sin(lonp- 3 * cPI / 2.);
        A22(iLR, iDU) = p2 * std::cos(latp);
      }
    }
  }
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
                              precision_t lat2)
{

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
                                           precision_t R)
{

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

arma_vec calc_bin_edges(arma_vec centers)
{

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

arma_mat calc_bin_edges(arma_mat centers, bool DoX)
{

  // X is first dimension (row), Y is second dimension (col)

  int64_t nX = centers.n_rows;
  int64_t nY = centers.n_cols;
  arma_mat edges;
  arma_vec centers1d;
  if (DoX)
  {
    if (verbose > 2)
      std::cout << "  --> x\n";
    edges.resize(nX + 1, nY);
    for (int64_t j = 0; j < nY; j++)
    {
      centers1d = centers.col(j);
      edges.col(j) = calc_bin_edges(centers1d);
    }
  }
  else
  {
    if (verbose > 2)
      std::cout << "  --> y\n";
    edges.resize(nX, nY + 1);
    for (int64_t i = 0; i < nX; i++)
    {
      centers1d = centers.row(i).as_col();
      edges.row(i) = calc_bin_edges(centers1d).as_row();
    }
  }
  return edges;
}

// ---------------------------------------------------------
// bin widths
// ---------------------------------------------------------

arma_vec calc_bin_widths(arma_vec edges)
{

  int64_t nPts = edges.n_elem - 1;
  arma_vec widths(nPts);

  for (int64_t i = 0; i < nPts; i++)
    widths(i) = edges(i + 1) - edges(i);

  return widths;
}

// ---------------------------------------------------------
// bin widths 2d
// ---------------------------------------------------------

arma_mat calc_bin_widths(arma_mat edges, bool DoX)
{

  int64_t nX = edges.n_rows;
  int64_t nY = edges.n_cols;

  arma_mat widths;
  arma_vec edges1d;

  if (DoX)
  {
    if (verbose > 2)
      std::cout << "  --> x\n";
    nX--;
    widths.resize(nX, nY);
    for (int64_t j = 0; j < nY; j++)
    {
      edges1d = edges.col(j);
      widths.col(j) = calc_bin_widths(edges1d);
    }
  }
  else
  {
    if (verbose > 2)
      std::cout << "  --> y\n";
    nY--;
    widths.resize(nX, nY);
    for (int64_t i = 0; i < nX; i++)
    {
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
                    int64_t nGCs)
{

  precision_t beta = 0.8;

  arma_vec s = left % right;
  arma_vec combined = (left + right) * 0.5;

  left = left * beta;
  right = right * beta;
  arma_vec limited = left;

  for (int64_t i = 1; i < nPts + 2 * nGCs - 1; i++)
  {
    if (s(i) < 0)
    {
      // Sign < 0 means opposite signed left and right:
      limited(i) = 0.0;
    }
    else
    {
      if (left(i) > 0 && right(i) > 0)
      {
        if (right(i) < limited(i))
          limited(i) = right(i);
        if (combined(i) < limited(i))
          limited(i) = combined(i);
      }
      else
      {
        if (right(i) > limited(i))
          limited(i) = right(i);
        if (combined(i) > limited(i))
          limited(i) = combined(i);
      }
    }
  }
  return limited;
}

void print(arma_vec values)
{
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
                      int64_t nGCs)
{

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
  //						    values.rows(iS-1, iE-1)) -
  //					 factor2 * (values.rows(iS+1, iE+1) -
  //						    values.rows(iS-2, iE-2)));
  //  hv.rows(iS, iE) = 2.0 / (x.rows(iS+1, iE+1) - x.rows(iS, iE));
  //  gradR.rows(iS, iE) = hv.rows(iS,iE) % (factor1 * (values.rows(iS+1, iE+1) -
  //						    values.rows(iS, iE)) -
  //					 factor2 * (values.rows(iS+2, iE+2) -
  //						    values.rows(iS-1, iE-1)));

  for (i = nGCs; i < nPts + nGCs; i++)
  {
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
                   bool DoX)
{

  arma_mat v2d, x2d;

  if (DoX)
  {
    v2d = values;
    x2d = x;
  }
  else
  {
    v2d = values.t();
    x2d = x.t();
  }

  int64_t nX = v2d.n_rows;
  int64_t nY = v2d.n_cols;
  arma_mat grad2d = v2d * 0.0;

  int64_t nPts = nX - 2 * nGCs;
  arma_vec values1d(nX);
  arma_vec x1d(nX);
  for (int64_t j = 1; j < nY - 1; j++)
  {
    values1d = v2d.col(j);
    x1d = x2d.col(j);
    grad2d.col(j) = calc_grad_1d(values1d, x1d, nPts, nGCs);
  }

  arma_mat gradients;

  if (DoX)
  {
    gradients = grad2d;
  }
  else
  {
    gradients = grad2d.t();
  }
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
                           int64_t nGCs)
{

  int64_t nX = values.n_rows;
  int64_t nY = values.n_cols;

  // Define at edges:
  arma_mat projected(nX + 1, nY);
  projected.zeros();

  // no gradient in the 0 or iEnd cells
  for (int64_t j = 0; j < nY; j++)
  {
    for (int64_t i = 1; i < nX - 1; i++)
    {
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
                            int64_t nGCs)
{
  int64_t nX = values.n_rows;
  int64_t nY = values.n_cols;

  // Define at edges:
  arma_mat projected(nX + 1, nY);
  projected.zeros();

  // no gradient in the 0 or iEnd cells
  for (int64_t j = 0; j < nY; j++)
  {
    for (int64_t i = 1; i < nX - 1; i++)
    {
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
                       int64_t nGCs)
{

  int64_t iStart = 0;
  int64_t iEnd = nPts + 2 * nGCs;

  arma_vec limited = projected;

  precision_t mini, maxi;

  for (int64_t i = iStart + 1; i < iEnd - 1; i++)
  {

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
                                   int64_t nGCs)
{

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
                  arma_mat &y)
{

  int64_t nX = x.n_rows;
  int64_t nY = x.n_cols;

  arma_mat rho(nX, nY);
  arma_mat r;

  r = sqrt((x - 0.0) % (x - 0.0) + (y - 0.0) % (y - 0.0));
  rho.fill(1.0);
  // rho.elem( find( r < 0.25)).fill(2.2);
  // rho.elem( find( r < 0.25)) = 2.25 - r.elem( find( r < 0.25));

  return rho;
}

// ---------------------------------------------------------
// initial velocity: initialize zero velocity
// ---------------------------------------------------------

arma_mat init_vel(arma_mat &x,
                  arma_mat &y)
{
  int64_t nX = x.n_rows;
  int64_t nY = x.n_cols;
  arma_mat vel(nX, nY);
  // all cells positive to right:
  vel.zeros();
  //vel.fill(0.1);
  return vel;
}

arma_mat init_vel2(arma_mat &x,
                   arma_mat &y)
{
  int64_t nX = x.n_rows;
  int64_t nY = x.n_cols;
  arma_mat vel(nX, nY);
  // all cells positive to right:
  //vel.zeros();
  vel.fill(0.1);
  return vel;
}

// ---------------------------------------------------------
// initial temp (E): constant total energy
// THIS IS NOT e but E, the total energy
// ---------------------------------------------------------

arma_mat init_temp(arma_mat &x,
                   arma_mat &y)
{
  int64_t nX = x.n_rows;
  int64_t nY = x.n_cols;

  arma_mat temp(nX, nY);
  double gamma = 5.0 / 3.0;
  temp.fill(1. / ((gamma - 1) * gamma));
  //temp.fill(1. / ((gamma - 1) * gamma));
  return temp;
}

int main()
{
  precision_t dt = 0.0001;        // Time Step
  precision_t current_time = 0.0; // Initial Time 0
  precision_t total_time = 2.0;   // Total simulation time
  precision_t cfl = 0.1;         // CFL Number
  precision_t gamma = 5.0 / 3.0;  // Specific ratio of heat

  precision_t dtOut = 0.01; // Output Interval

  int64_t nSteps = 50; // Number of Total Time steps
  int64_t iStep;       // Iterator of Time Step

  int64_t nX = 48;  // Number of x grid cells
  int64_t nY = nX;  // Number of y grid cells
  int64_t nGCs = 2; // Number of ghost cells

  // Local Coordinates are nX+nGCs*2 by nY+nGCs*2
  arma_mat x(nX + 2 * nGCs, nY + 2 * nGCs);
  arma_mat y(nX + 2 * nGCs, nY + 2 * nGCs);

  // Local Latitude and Longitude are nX+nGCs*2 by nY+nGCs*2
  arma_mat lat(nX + 2 * nGCs, nY + 2 * nGCs);
  arma_mat lon(nX + 2 * nGCs, nY + 2 * nGCs);

  // Metrics regarding each cell are the same size as the number of cells
  arma_mat jacobian(nX + 2 * nGCs, nY + 2 * nGCs);
  arma_mat A11_inv(nX + 2 * nGCs, nY + 2 * nGCs);
  arma_mat A12_inv(nX + 2 * nGCs, nY + 2 * nGCs);
  arma_mat A21_inv(nX + 2 * nGCs, nY + 2 * nGCs);
  arma_mat A22_inv(nX + 2 * nGCs, nY + 2 * nGCs);
  arma_mat A11(nX + 2 * nGCs, nY + 2 * nGCs);
  arma_mat A12(nX + 2 * nGCs, nY + 2 * nGCs);
  arma_mat A21(nX + 2 * nGCs, nY + 2 * nGCs);
  arma_mat A22(nX + 2 * nGCs, nY + 2 * nGCs);

  // Radius of Sphere
  precision_t R = 1.;

  if (verbose > 0)
    std::cout << " --> generating cubesphere cell center and metrics\n";

  /** Generate cell centers */
  generate_x_y_cubesphere(4, nX + 2 * nGCs, nY + 2 * nGCs, nGCs, R, 0.5, 0.5,
                          lon, lat, jacobian, A11_inv, A12_inv, A21_inv, A22_inv, A11, A12, A21, A22, x, y, true);

  /** Generate widths */
  if (verbose > 0)
    std::cout << " --> calculating dx dy\n";
  // Since it is equidistant
  precision_t dx = std::abs(x(0, 0) - x(1, 0));
  precision_t dy = std::abs(y(0, 0) - y(0, 1));

  /** Generate cell area (in reference space) */
  if (verbose > 0)
    std::cout << " --> calculating area\n";
  precision_t area = dx * dy;

  /** Calculate edge midpoint locations in reference space */
  if (verbose > 0)
    std::cout << " --> calculating edges\n";
  arma_mat xEdges = calc_bin_edges(x, true);
  arma_mat yEdges = calc_bin_edges(y, false);

  if (verbose > 0)
    std::cout << " --> calculating edge metrics\n";
  // Metrics regarding each cell edges are the same size as the number of edges
  arma_mat jacobian_xEdges(nX + 2 * nGCs + 1, nY + 2 * nGCs);

  // Remember, they are in fact metric tensor with upper indices.
  arma_mat g11_xEdges(nX + 2 * nGCs + 1, nY + 2 * nGCs);
  arma_mat g12_xEdges(nX + 2 * nGCs + 1, nY + 2 * nGCs);
  arma_mat g21_xEdges(nX + 2 * nGCs + 1, nY + 2 * nGCs);
  arma_mat g22_xEdges(nX + 2 * nGCs + 1, nY + 2 * nGCs);
  edge_metrics(4, xEdges, y, R, jacobian_xEdges, g11_xEdges,
               g12_xEdges, g21_xEdges, g22_xEdges);

  arma_mat jacobian_yEdges(nX + 2 * nGCs, nY + 2 * nGCs + 1);
  arma_mat g11_yEdges(nX + 2 * nGCs, nY + 2 * nGCs + 1);
  arma_mat g12_yEdges(nX + 2 * nGCs, nY + 2 * nGCs + 1);
  arma_mat g21_yEdges(nX + 2 * nGCs, nY + 2 * nGCs + 1);
  arma_mat g22_yEdges(nX + 2 * nGCs, nY + 2 * nGCs + 1);
  edge_metrics(4, x, yEdges, R, jacobian_yEdges, g11_yEdges,
               g12_yEdges, g21_yEdges, g22_yEdges);

  /** State Initialization */
  /// Initialize Density
  if (verbose > 0)
    std::cout << "---> initializing rho\n";
  arma_mat rho = init_rho(x, y); // rho, pure scalar field

  /// Initialize Velocity and Momentum
  if (verbose > 0)
    std::cout << "---> initializing vel\n";
  // Initialize spherical velocity
  arma_mat xVelSph = init_vel(x, y);
  arma_mat yVelSph = init_vel2(x, y);
  // Transform to reference/local velocity (contravariant velocity)
  arma_mat xVel = xVelSph % A11_inv + yVelSph % A12_inv; // u^1
  arma_mat yVel = xVelSph % A21_inv + yVelSph % A22_inv; // u^2

  if (verbose > 0)
    std::cout << "---> initializing momentum\n";
  // The following two are used for iteration
  arma_mat xMomentum = rho % xVel; // x1momentum, pure scalar field
  arma_mat yMomentum = rho % yVel; // x2momentum, pure scalar field

  /// Initialize total energy
  if (verbose > 0)
    std::cout << "---> initializing energy\n";
  arma_mat rhoE = rho % init_temp(x, y); // rhoE, pure scalar field

  /** Initialize projection constructs */
  projection_struct rhoP;
  projection_struct xMomentumP;
  projection_struct yMomentumP;
  projection_struct rhoEP;

  // They are all pure scalar fields without sqrt(g)
  arma_mat rhoL, rhoR, rhoD, rhoU;
  arma_mat xVelL, xVelR, xVelD, xVelU;
  arma_mat yVelL, yVelR, yVelD, yVelU;
  arma_mat totaleL, totaleR, totaleD, totaleU;

  arma_mat velL2, velR2, velD2, velU2;
  arma_mat tempL, tempR, tempD, tempU;
  arma_mat pressureL, pressureR, pressureD, pressureU;
  
  /** Initialize Flux and Wave Speed Storages */
  arma_mat eq1FluxLR, eq1FluxDU;
  arma_mat eq1FluxL, eq1FluxR, eq1FluxD, eq1FluxU;

  arma_mat eq2FluxLR, eq2FluxDU;
  arma_mat eq2FluxL, eq2FluxR, eq2FluxD, eq2FluxU;

  arma_mat eq3FluxLR, eq3FluxDU;
  arma_mat eq3FluxL, eq3FluxR, eq3FluxD, eq3FluxU;

  arma_mat eq4FluxLR, eq4FluxDU;
  arma_mat eq4FluxL, eq4FluxR, eq4FluxD, eq4FluxU;

  arma_mat wsL, wsR, wsD, wsU, wsLR, wsDU;

  arma_mat diff; // for Riemann Solver

  /** SOME BC EXCHANGE CODE THAT DOES NOT EXIST (Dirichilet BC) */

  /** Output some pre-simulation results */
  if (verbose > 0)
    std::cout << "---> outputting\n";
  output(x, "x.txt", false);
  output(y, "y.txt", false);
  output(xEdges, "xEdges.txt", false);
  output(yEdges, "yEdges.txt", false);
  output(lat, "lat.txt", false);
  output(lon, "lon.txt", false);
  output(rho, "rho.txt", false);
  output(xMomentum, "xMomentum.txt", false);
  output(yMomentum, "yMomentum.txt", false);
  output(xVel, "xVel.txt", false);
  output(yVel, "yVel.txt", false);
  output(rhoE, "rhoE.txt", false);

  iStep = 0;

  while (current_time < total_time)
  {

    if (verbose > 0)
      std::cout << "step : " << iStep << "; time : " << current_time << "\n";

    iStep++;

    // -----------------------------------
    /** Project State Variables (pure scalar fields with pulled sqrt(g)) to Edges: */

    if (verbose > 3)
      std::cout << "Projecting\n";

    rhoP = project_to_edges(rho, x, xEdges, y, yEdges, nGCs);
    xMomentumP = project_to_edges(xMomentum, x, xEdges, y, yEdges, nGCs);
    yMomentumP = project_to_edges(yMomentum, x, xEdges, y, yEdges, nGCs);
    rhoEP = project_to_edges(rhoE, x, xEdges, y, yEdges, nGCs);

    // Resolve Scalar Fields into rho, xVel, yVel, and totalE (without rho)
    rhoL = rhoP.L;
    rhoR = rhoP.R;
    rhoD = rhoP.D;
    rhoU = rhoP.U;

    xVelL = xMomentumP.L / rhoL;
    xVelR = xMomentumP.R / rhoR;
    xVelD = xMomentumP.D / rhoD;
    xVelU = xMomentumP.U / rhoU;

    yVelL = yMomentumP.L / rhoL;
    yVelR = yMomentumP.R / rhoR;
    yVelD = yMomentumP.D / rhoD;
    yVelU = yMomentumP.U / rhoU;

    totaleL = rhoEP.L / rhoL;
    totaleR = rhoEP.R / rhoR;
    totaleD = rhoEP.D / rhoD;
    totaleU = rhoEP.U / rhoU;

    velL2 = (xVelL % xVelL + yVelL % yVelL);
    velR2 = (xVelR % xVelR + yVelR % yVelR);
    velD2 = (xVelD % xVelD + yVelD % yVelD);
    velU2 = (xVelU % xVelU + yVelU % yVelU);

    tempL = totaleL - 0.5 * velL2;
    tempR = totaleR - 0.5 * velR2;
    tempD = totaleD - 0.5 * velD2;
    tempU = totaleU - 0.5 * velU2;

    pressureL = (gamma - 1) * (rhoP.L % tempL);
    pressureR = (gamma - 1) * (rhoP.R % tempR);
    pressureD = (gamma - 1) * (rhoP.D % tempD);
    pressureU = (gamma - 1) * (rhoP.U % tempU);

    // ------------------------------------------------
    // Calculate fluxes of different terms at the edges:
    // All fluxes are in spherical sense, local coords

    if (verbose > 3)
      std::cout << "Calculating Fluxes\n";
    
    // Note that dot product between normal vector at edge and flux vector
    // resolves into a pure one component flux or either hat{x} or hat{y}

    // Flux calculated from the left of the edge
    eq1FluxL = rhoL % xVelL % jacobian_xEdges;
    // Flux calculated from the right of the edge
    eq1FluxR = rhoR % xVelR % jacobian_xEdges;
    // Flux calculated from the down of the edge
    eq1FluxD = rhoD % yVelD % jacobian_yEdges;
    // Flux calculated from the up of the edge
    eq1FluxU = rhoU % yVelU % jacobian_yEdges;

    eq2FluxL = (rhoL % xVelL % xVelL + pressureL % g11_xEdges) % jacobian_xEdges;
    eq2FluxR = (rhoR % xVelR % xVelR + pressureR % g11_xEdges) % jacobian_xEdges;
    eq2FluxD = (rhoD % yVelD % xVelD + pressureD % g12_yEdges) % jacobian_yEdges;
    eq2FluxU = (rhoU % yVelU % xVelU + pressureU % g12_yEdges) % jacobian_yEdges;

    eq3FluxL = (rhoL % xVelL % yVelL + pressureL % g21_xEdges) % jacobian_xEdges;
    eq3FluxR = (rhoR % xVelR % yVelR + pressureR % g21_xEdges) % jacobian_xEdges;
    eq3FluxD = (rhoD % yVelD % yVelD + pressureD % g22_yEdges) % jacobian_yEdges;
    eq3FluxU = (rhoU % yVelU % yVelU + pressureU % g22_yEdges) % jacobian_yEdges;

    eq4FluxL = (totaleL % rhoL + pressureL) % xVelL % jacobian_xEdges;
    eq4FluxR = (totaleR % rhoR + pressureR) % xVelR % jacobian_xEdges;
    eq4FluxD = (totaleD % rhoD + pressureD) % yVelD % jacobian_yEdges;
    eq4FluxU = (totaleU % rhoU + pressureU) % yVelU % jacobian_yEdges;

    // ------------------------------------------------
    // Calculate the wave speed for the diffusive flux:
    // In Reference velocities

    if (verbose > 3)
      std::cout << "Diffusive Fluxes\n";
    wsL = sqrt(velL2) + sqrt(gamma * (gamma - 1) * tempL);
    wsR = sqrt(velR2) + sqrt(gamma * (gamma - 1) * tempR);
    wsD = sqrt(velD2) + sqrt(gamma * (gamma - 1) * tempD);
    wsU = sqrt(velU2) + sqrt(gamma * (gamma - 1) * tempU);

    wsLR = wsR;
    for (int64_t i = 0; i < nX + 1; i++)
    {
      for (int64_t j = 0; j < nY; j++)
      {
        if (wsL(i, j) > wsLR(i, j))
          wsLR(i, j) = wsL(i, j);
      }
    }

    wsDU = wsD;
    for (int64_t i = 0; i < nX; i++)
    {
      for (int64_t j = 0; j < nY + 1; j++)
      {
        if (wsU(i, j) > wsDU(i, j))
          wsDU(i, j) = wsU(i, j);
      }
    }

    // ------------------------------------------------
    // Calculate dt based on max waves speeds and cell sizes

    if (verbose > 3)
      std::cout << "Calculating dt\n";

    dt = calc_dt(dx, dy, wsLR, wsDU, nGCs);
    dt = cfl * dt;
    current_time += dt;

    // ------------------------------------------------
    // Calculate average flux at the edges (Rusanov Flux):

    if (verbose > 3)
      std::cout << "Averaging fluxes at edges\n";

    diff = (rhoR - rhoL) % jacobian_xEdges; // State difference, need to add sqrt(g)
    eq1FluxLR = (eq1FluxL + eq1FluxR) / 2 + 0.5 * wsLR % diff;
    diff = (rhoU - rhoD) % jacobian_yEdges;
    eq1FluxDU = (eq1FluxD + eq1FluxU) / 2 + 0.5 * wsDU % diff;

    diff = (rhoR % xVelR - rhoL % xVelL) % jacobian_xEdges;
    eq2FluxLR = (eq2FluxL + eq2FluxR) / 2 + 0.5 * wsLR % diff;
    diff = (rhoU % xVelU - rhoD % xVelD) % jacobian_yEdges;
    eq2FluxDU = (eq2FluxD + eq2FluxU) / 2 + 0.5 * wsDU % diff;

    diff = (rhoR % yVelR - rhoL % yVelL) % jacobian_xEdges;
    eq3FluxLR = (eq3FluxL + eq3FluxR) / 2 + 0.5 * wsLR % diff;
    diff = (rhoU % yVelU - rhoD % yVelD) % jacobian_yEdges;
    eq3FluxDU = (eq3FluxD + eq3FluxU) / 2 + 0.5 * wsDU % diff;

    diff = (rhoR % totaleR - rhoL % totaleL) % jacobian_xEdges;
    eq4FluxLR = (eq4FluxL + eq4FluxR) / 2 + 0.5 * wsLR % diff;
    diff = (rhoU % totaleU - rhoD % totaleD) % jacobian_yEdges;
    eq4FluxDU = (eq4FluxD + eq4FluxU) / 2 + 0.5 * wsDU % diff;

    // ------------------------------------------------
    // Update values:
    if (verbose > 3)
      std::cout << "Updating equations of state\n";

    // Total residual norm for this iteration (for L1 Residua; Norm Calc)
    precision_t residual_ij = 0;

    // Only deal with inner cell
    for (int64_t j = nGCs; j < nY + nGCs; j++)
    {
      for (int64_t i = nGCs; i < nX + nGCs; i++)
      {
        precision_t rhoResidual_ij = (dy * eq1FluxLR(i + 1, j) -
                                      dy * eq1FluxLR(i, j) +
                                      dx * eq1FluxDU(i, j + 1) -
                                      dx * eq1FluxDU(i, j));
        rho(i, j) = rho(i, j) - dt / area / jacobian(i, j) * rhoResidual_ij;

        if (i == 25 && j == 2)
        {
          std::cout << (dy * eq1FluxLR(i + 1, j) -
                        dy * eq1FluxLR(i, j) +
                        dx * eq1FluxDU(i, j + 1) -
                        dx * eq1FluxDU(i, j))
                    << std::endl;
        }


        precision_t xMomentumResidual_ij = dy * eq2FluxLR(i + 1, j) -
                                           dy * eq2FluxLR(i, j) +
                                           dx * eq2FluxDU(i, j + 1) -
                                           dx * eq2FluxDU(i, j);
        xMomentum(i, j) = xMomentum(i, j) - dt / area / jacobian(i, j) * xMomentumResidual_ij;
        if (i == 25 && j == 2)
        {
          std::cout << (dy * eq2FluxLR(i + 1, j) -
                        dy * eq2FluxLR(i, j) +
                        dx * eq2FluxDU(i, j + 1) -
                        dx * eq2FluxDU(i, j))
                    << std::endl;
        }

        precision_t yMomentumResidual_ij = (dy * eq3FluxLR(i + 1, j) -
                                            dy * eq3FluxLR(i, j) +
                                            dx * eq3FluxDU(i, j + 1) -
                                            dx * eq3FluxDU(i, j));
        if (i == 25 && j == 2)
        {
          std::cout << (dy * eq3FluxLR(i + 1, j) -
                        dy * eq3FluxLR(i, j) +
                        dx * eq3FluxDU(i, j + 1) -
                        dx * eq3FluxDU(i, j))
                    << std::endl;
        }
        yMomentum(i, j) = yMomentum(i, j) - dt / area / jacobian(i, j) * yMomentumResidual_ij;

        precision_t rhoEResidual_ij = (dy * eq4FluxLR(i + 1, j) -
                                         dy * eq4FluxLR(i, j) +
                                         dx * eq4FluxDU(i, j + 1) -
                                         dx * eq4FluxDU(i, j));
        rhoE(i, j) = rhoE(i,j) - dt / area / jacobian(i, j) * rhoEResidual_ij;

        if (i == 25 && j == 2)
        {
          std::cout << (dy * eq4FluxLR(i + 1, j) -
                        dy * eq4FluxLR(i, j) +
                        dx * eq4FluxDU(i, j + 1) -
                        dx * eq4FluxDU(i, j))
                    << std::endl;
        }

        residual_ij += std::abs(rhoResidual_ij) + std::abs(xMomentumResidual_ij) +
                       std::abs(yMomentumResidual_ij) + std::abs(rhoEResidual_ij);
      }
    }
    std::cout << "L1 Residual Norm : " << residual_ij << "\n";

    // ------------------------------------------------
    // Exchange messages (set BCs, really):

    if (verbose > 3)
      std::cout << "Outputing\n";

    // if (int((current_time - dt)/dtOut) != int((current_time )/dtOut)) {
    std::cout << "Outputing at time : " << current_time << "\n";
    output(rho, "rho.txt", true);
    output(xVel, "xvel.txt", true);
    output(yVel, "yvel.txt", true);
    output(rhoE, "rhoe.txt", true);
    // output(rhoL, "rhor.txt", false);
    // output(rhoR, "rhol.txt", false);
    // }
  }

  arma_mat xVelSph_output = xVel%A11 + yVel%A12;
  arma_mat yVelSph_output = xVel%A21 + yVel%A22;
  output(xVelSph_output, "xvel_sph.txt", true);
  output(yVelSph_output, "yvel_sph.txt", true);
  return 0;
}