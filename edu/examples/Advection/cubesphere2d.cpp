
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
precision_t cTWOPI = 2.0*cPI;
precision_t cPIdiv2 = cPI/2;

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

// ---------------------------------------------------------
// grid creation
// ---------------------------------------------------------

void generate_x_y_cubesphere(int64_t nX, int64_t nY, int64_t nGCs,
			     precision_t left_off,
			     precision_t down_off,
			     arma_mat &lon2d, arma_mat &lat2d,
			     arma_mat &jacobian,
			     arma_mat &gP11,
			     arma_mat &gP12,
			     arma_mat &gP21,
			     arma_mat &gP22,
			     arma_mat &x,
			     arma_mat &y) {

  double xn, yn, zn, rn;
  double xp, yp, zp, rp, latp, lonp, p;
  double dr, du;
  
  double a = sqrt(3);

  arma_vec xyz, xyzn, xyz_wrapped;

  dr = 2.0 / (nX - 2.0 * nGCs);
  du = 2.0 / (nY - 2.0 * nGCs);

  double g11, g12, g21, g22, detg, R, g;

  R = sqrt(3) * a;
  
  // Loop through each point and derive the coordinate
  for (int iDU = 0; iDU < nY; iDU++) {
    for (int iLR = 0; iLR < nX; iLR++) {

      // the offsets are so we can find cell centers, edges, and corners
      double iD = iDU - nGCs + down_off;
      double iL = iLR - nGCs + left_off;

      xp = (-1.0 + dr * iL) * a;
      yp = -a;
      zp = (-1.0 + du * iD) * a;
      rp = sqrt(xp * xp + yp * yp + zp * zp);
      xn = xp / rp;
      yn = yp / rp;
      zn = zp / rp;
      latp = asin(zn);
      lonp = atan2(yn, xn) + cPI/2.0;

      lat2d(iLR, iDU) = latp;
      lon2d(iLR, iDU) = lonp;
      x(iLR, iDU) = xp;
      y(iLR, iDU) = zp;

      /*
	This stuff is from:
	A Discontinuous Galerkin Transport Scheme on the Cubed Sphere
	RAMACHANDRAN D. NAIR, STEPHEN J. THOMAS, AND RICHARD D. LOFT
	MONTHLY WEATHER REVIEW - VOLUME 133
	2004
      */
      
      // jacobian = R^2 * a / r^3
      jacobian(iLR, iDU) = R * R * a / (rp * rp * rp);

      // (28) - in this case, z = y
      p = R * R / (rp * rp * rp * rp);
      g11 = p * (a * a + zp * zp);
      g12 = -p * xp * zp;
      g21 = g12;
      g22 = p * (a * a + xp * xp);

      // (29) - in this case, z = y

      //g = g11 * g22 - g21 * g12;
      //
      //p = 1 / g;
      //gP11(iLR, iDU) = p * g22;
      //gP12(iLR, iDU) = -p * g12;
      //gP21(iLR, iDU) = -p * g21;
      //gP22(iLR, iDU) = p * g11;

      p = a / (R * cos(latp) * cos(lonp));
      gP11(iLR, iDU) = p * (1 / cos(lonp));
      gP12(iLR, iDU) = 0.0;
      gP21(iLR, iDU) = p * (tan(latp) * tan(lonp));
      gP22(iLR, iDU) = p * (1 / cos(latp));
      /*
      std::cout << iDU << " "
		<< iLR << " "
		<< lonp << " "
		<< xp << " "
		<< zp << " "
		<< g << " "
		<< jacobian(iLR, iDU) * jacobian(iLR, iDU) << "\n";
      */
    }
  }
  
  return;
}

// ---------------------------------------------------------
// calculate angles between coordinates
// ---------------------------------------------------------

arma_mat calc_grid_angles(arma_mat x, arma_mat y) {

  int64_t nX = x.n_rows;
  int64_t nY = x.n_cols;
  arma_mat angles(nX, nY);

  precision_t dx_ud, dy_ud;
  precision_t dx_lr, dy_lr;
  precision_t m;
  
  for (int iDU = 1; iDU < nY-1; iDU++) {
    for (int iLR = 1; iLR < nX-1; iLR++) {
      dx_lr = (x(iLR+1, iDU) - x(iLR-1, iDU))/2.0;
      dy_lr = (y(iLR+1, iDU) - y(iLR-1, iDU))/2.0;
      m = sqrt(dx_lr*dx_lr + dy_lr*dy_lr);
      dx_lr = dx_lr/m;
      dy_lr = dy_lr/m;

      dx_ud = (x(iLR, iDU+1) - x(iLR, iDU-1))/2.0;
      dy_ud = (y(iLR, iDU+1) - y(iLR, iDU-1))/2.0;
      m = sqrt(dx_ud*dx_ud + dy_ud*dy_ud);
      dx_ud = dx_ud/m;
      dy_ud = dy_ud/m;

      angles(iLR, iDU) = acos(dx_lr * dx_ud + dy_lr * dy_ud);
      std::cout << iDU << " "
		<< iLR << " "
		<< sin(angles(iLR, iDU)) << "\n";
      
    }
  }

  angles.row(0) = angles.row(1);
  angles.row(nX-1) = angles.row(nX-2);
  angles.col(0) = angles.col(1);
  angles.col(nY-1) = angles.col(nY-2);
  
  return angles;
  
}



// ---------------------------------------------------------
// bin edges
// ---------------------------------------------------------

arma_vec calc_bin_edges(arma_vec centers) {

  int64_t nPts = centers.n_elem;
  arma_vec edges(nPts+1);

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
    if (verbose > 2) std::cout << "  --> x\n";
    edges.resize(nX+1, nY);
    for (int64_t j = 0; j < nY; j++) {
      centers1d = centers.col(j);
      edges.col(j) = calc_bin_edges(centers1d);
    }
  } else {
    if (verbose > 2) std::cout << "  --> y\n";
    edges.resize(nX, nY+1);
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
    if (verbose > 2) std::cout << "  --> x\n";
    nX--;
    widths.resize(nX, nY);
    for (int64_t j = 0; j < nY; j++) {
      edges1d = edges.col(j);
      widths.col(j) = calc_bin_widths(edges1d);
    }
  } else {
    if (verbose > 2) std::cout << "  --> y\n";
    nY--;
    widths.resize(nX, nY);
    for (int64_t i = 0; i < nX; i++) {
      edges1d = edges.row(i).as_col();
      widths.row(i) = calc_bin_widths(edges1d).as_row();
    }
  }
  return widths;
}

// ---------------------------------------------------------
// initial rho
// ---------------------------------------------------------

arma_mat init_rho(arma_mat &x,
		  arma_mat &y) {

  int64_t nX = x.n_rows;
  int64_t nY = x.n_cols;

  arma_mat rho(nX, nY);
  arma_mat r;

  r = sqrt( (x - 0.0) % (x - 0.0) + (y - 0.0) % (y - 0.0));
  rho.fill(2.0);
  rho.elem( find( r < 0.25)).fill(12.2);
  //rho.elem( find( r < 0.25)) = 2.25 - r.elem( find( r < 0.25));
  
  return rho;
}

// ---------------------------------------------------------
// initial velocity
// ---------------------------------------------------------

arma_mat init_vel(arma_mat &x,
		  arma_mat &y) {
  int64_t nX = x.n_rows;
  int64_t nY = x.n_cols;
  arma_mat vel(nX, nY);
  // all cells positive to right:
  vel.zeros();
  return vel;
}

// ---------------------------------------------------------
// initial temp (e)
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
// exchange messages
// ---------------------------------------------------------

void exchange(arma_mat &values, int64_t nGCs) {

  int64_t nX = values.n_rows;
  int64_t nY = values.n_cols;

  // this is a periodic BC:
  if (verbose > 2) std::cout << "  --> x\n";
  for (int64_t j = 0; j < nY; j++) {
    for (int64_t i = 0; i < nGCs; i++) {
      values(i, j) = values(nX - 2 * nGCs + i, j);
      values(nX - nGCs + i, j) = values(nGCs + i, j);
    }
  }
  if (verbose > 2) std::cout << "  --> y\n";
  for (int64_t i = 0; i < nX; i++) {
    for (int64_t j = 0; j < nGCs; j++) {
      values(i, j) = values(i, nY - 2 * nGCs + j);
      values(i, nY - nGCs + j) = values(i, nGCs + j);
    }
  }
}

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
  for (int64_t i = 0 ; i < nP ; i++)
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
  h = 2.0 / (x(i+1) - x(i));
  gradR(i) = h * (factor1 * (values(i+1) - values(i)) -
		  factor2 * (values(i+2) - values(i-1)));
  gradL(i) = (values(i) - values(i-1)) / (x(i) - x(i-1));

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
  
  for (i = nGCs; i < nPts + nGCs; i++) {
    h = 2.0 / (x(i) - x(i-1));
    gradL(i) = h * (factor1 * (values(i) - values(i-1)) -
		    factor2 * (values(i+1) - values(i-2)));
    h = 2.0 / (x(i+1) - x(i));
    gradR(i) = h * (factor1 * (values(i+1) - values(i)) -
		    factor2 * (values(i+2) - values(i-1)));
  }
  i = nPts + nGCs;
  h = 2.0 / (x(i) - x(i-1));
  gradL(i) = h * (factor1 * (values(i) - values(i-1)) -
		  factor2 * (values(i+1) - values(i-2)));
  gradR(i) = (values(i+1) - values(i)) / (x(i+1) - x(i));
  
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
  for (int64_t j = 1; j < nY-1; j++) {
    values1d = v2d.col(j);
    x1d = x2d.col(j);
    grad2d.col(j) = calc_grad_1d(values1d, x1d, nPts, nGCs);
  }

  arma_mat gradients;

  if (DoX) {
    gradients = grad2d;
  } else {
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

    mini = values(i-1);
    if (values(i) < mini)
      mini = values(i);
    maxi = values(i-1);
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
			     y_centers.t(), y_edges.t(), nGCs).t();
  // Down side of edge from up (right)
  proj.D = project_from_right(values.t(), proj.gradDU.t(),
			      y_centers.t(), y_edges.t(), nGCs).t();

  return proj;
}

// ---------------------------------------------------------
//
// ---------------------------------------------------------

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

// ---------------------------------------------------------
//
// ---------------------------------------------------------

precision_t calc_dt(arma_mat &xWidth,
		    arma_mat &yWidth,
		    arma_mat &wsLR,
		    arma_mat &wsDU,
		    int64_t nGCs) {

  int64_t nX = xWidth.n_rows;
  int64_t nY = yWidth.n_cols;

  precision_t wsX, wsY, dtX, dtY, dt;

  dt = 1e32;
  
  for (int64_t j = nGCs; j < nY - nGCs; j++) {
    for (int64_t i = nGCs; i < nX - nGCs; i++) {
      wsX = (wsLR(i+1, j) + wsLR(i, j))/2;
      dtX = xWidth(i, j) / wsX;
      wsY = (wsDU(i, j+1) + wsDU(i, j))/2;
      dtY = yWidth(i, j) / wsY;
      if (dtX < dt) dt = dtX;
      if (dtY < dt) dt = dtY;
    }
  }
  return dt;
}

// ---------------------------------------------------------
// Calculate the normal to a line segment
// ---------------------------------------------------------

precision_t calc_norm(precision_t lon1,
		      precision_t lon2,
		      precision_t lat1,
		      precision_t lat2) {

  precision_t dlon = lon2 * cos(lat2) - lon1 * cos(lat1);
  precision_t dlat = lat2 - lat1;
  precision_t len = sqrt(dlon * dlon + dlat * dlat);
  
  arma_vec norm(2);
  // 90-deg CCW rotation:
  norm(0) = dlat / len;
  norm(1) = -dlon / len;
  
  return norm;
}

// ---------------------------------------------------------
// Calculate great circle distance
// ---------------------------------------------------------

precision_t calc_great_circle(precision_t lon1,
			      precision_t lon2,
			      precision_t lat1,
			      precision_t lat2) {

  precision_t dlon_2 = (lon2 - lon1)/2.0;
  precision_t dlat_2 = (lat2 - lat1)/2.0;
  
  precision_t dh = 2.0 *
    asin(sqrt(sin(dlat_2) * sin(dlat_2) +
	      sin(dlon_2) * sin(dlon_2) * cos(lat1) * cos(lat2)));

  return dh;
}

// ---------------------------------------------------------
// x, y, z from r, lat, lon
// ---------------------------------------------------------

arma_vec calc_xyz_from_lon_lat_r(precision_t lon,
				  precision_t lat,
				  precision_t r) {
  arma_vec xyz(3);
  xyz(0) = r * cos(lat) * cos(lon);
  xyz(1) = r * cos(lat) * sin(lon);
  xyz(2) = r * sin(lat);
  return xyz;
}

// ---------------------------------------------------------
// cross product
// ---------------------------------------------------------

arma_vec calc_cross(arma_vec v1,
		    arma_vec v2) {

  arma_vec c(3);
  c(0) = v1(1) * v2(2) - v2(1) * v1(2);
  c(1) = - (v1(0) * v2(2) - v2(0) * v1(2));
  c(2) = v1(0) * v2(1) - v2(0) * v1(1);
  return c;
}

// ---------------------------------------------------------
// Angle between three points on a sphere
// ---------------------------------------------------------

precision_t calc_angle_given_three_lon_lat(precision_t p1_lon,
					   precision_t p1_lat,
					   precision_t p2_lon,
					   precision_t p2_lat,
					   precision_t p3_lon,
					   precision_t p3_lat) {

  arma_vec p1 = calc_xyz_from_lon_lat_r(p1_lon, p1_lat, 1.0);
  arma_vec p2 = calc_xyz_from_lon_lat_r(p2_lon, p2_lat, 1.0);
  arma_vec p3 = calc_xyz_from_lon_lat_r(p3_lon, p3_lat, 1.0);
  arma_vec d1 = p1 - p2;
  arma_vec d2 = p3 - p2;
  precision_t n1 = sqrt(d1(0)*d1(0) +
			d1(1)*d1(1) +
			d1(2)*d1(2));
  precision_t n2 = sqrt(d2(0)*d2(0) +
			d2(1)*d2(1) +
			d2(2)*d2(2));
  d1 = d1 / n1;
  d2 = d2 / n2;
  precision_t angle = acos(d1(0) * d2(0) +
			   d1(1) * d2(1) +
			   d1(2) * d2(2));
  
  return angle;
}

// ---------------------------------------------------------
// main code
// ---------------------------------------------------------

int main() {

  precision_t dt = 0.0005/4;
  precision_t current_time = 0.0;
  precision_t total_time = 1.0;
  precision_t cfl = 0.1;
  precision_t gamma = 5.0/3.0;

  precision_t dtOut = 0.01;
  
  int64_t nSteps = 50;
  int64_t iStep;
  
  int64_t nX = 45;
  int64_t nY = nX;
  int64_t nGCs = 2;

  if (verbose > 0) std::cout << "---> initializing grid\n";
  arma_mat x(nX + 2*nGCs, nY + 2*nGCs);
  arma_mat y(nX + 2*nGCs, nY + 2*nGCs);
  arma_mat xCorner(nX + 2*nGCs + 1, nY + 2*nGCs + 1);
  arma_mat yCorner(nX + 2*nGCs + 1, nY + 2*nGCs + 1);
  arma_mat lat(nX + 2*nGCs, nY + 2*nGCs);
  arma_mat lon(nX + 2*nGCs, nY + 2*nGCs);
  arma_mat latCorner(nX + 2*nGCs + 1, nY + 2*nGCs + 1);
  arma_mat lonCorner(nX + 2*nGCs + 1, nY + 2*nGCs + 1);

  arma_mat jacobian(nX + 2*nGCs, nY + 2*nGCs);
  arma_mat g11(nX + 2*nGCs, nY + 2*nGCs);
  arma_mat g12(nX + 2*nGCs, nY + 2*nGCs);
  arma_mat g21(nX + 2*nGCs, nY + 2*nGCs);
  arma_mat g22(nX + 2*nGCs, nY + 2*nGCs);

  arma_mat DummyJacobian(nX + 2*nGCs + 1, nY + 2*nGCs + 1);
  arma_mat DummyG11(nX + 2*nGCs + 1, nY + 2*nGCs + 1);
  arma_mat DummyG12(nX + 2*nGCs + 1, nY + 2*nGCs + 1);
  arma_mat DummyG21(nX + 2*nGCs + 1, nY + 2*nGCs + 1);
  arma_mat DummyG22(nX + 2*nGCs + 1, nY + 2*nGCs + 1);

  arma_mat dxDU(nX + 2*nGCs, nY + 2*nGCs + 1);
  arma_mat dyLR(nX + 2*nGCs + 1, nY + 2*nGCs);

  // ---------------------------------------------
  // Cell Centers
  // ---------------------------------------------
  generate_x_y_cubesphere(nX + 2*nGCs, nY+ 2*nGCs, nGCs, 0.5, 0.5,
			  lon, lat, jacobian, g11, g12, g21, g22, x, y);

  generate_x_y_cubesphere(nX + 2*nGCs + 1, nY+ 2*nGCs + 1, nGCs,
			  0.0, 0.0,
			  lonCorner, latCorner,
			  DummyJacobian,
			  DummyG11, DummyG12, DummyG21, DummyG22,
			  xCorner, yCorner);

  // calc dx, defined as:
  // distance between lower-left corner and lower-right corner
  for (int64_t j = 0; j < nY + 1; j++) {
    for (int64_t i = 0; i < nX; i++) {
      dxDU(i, j) = 3 * calc_great_circle(lonCorner(i, j),
					 lonCorner(i+1, j),
					 latCorner(i, j),
					 latCorner(i+1, j));
    }
  }

  // calc dy, defined as:
  // distance between lower-left corner and upper-left corner
  for (int64_t j = 0; j < nY + 1; j++) {
    for (int64_t i = 0; i < nX; i++) {
      dyLR(i, j) = 3 * calc_great_circle(lonCorner(i, j),
					 lonCorner(i, j+1),
					 latCorner(i, j),
					 latCorner(i, j+1));
    }
  }


  // arma_mat angles = calc_grid_angles(x, y);

  if (verbose > 0) std::cout << " --> calculating edges\n";
  arma_mat xEdges = calc_bin_edges(x, true);
  arma_mat yEdges = calc_bin_edges(y, false);

  if (verbose > 0) std::cout << " --> calculating widths\n";
  arma_mat xWidth = calc_bin_widths(xEdges, true);
  arma_mat yWidth = calc_bin_widths(yEdges, false);

  arma_mat xDt, yDt;
  
  arma_mat area = xWidth % yWidth;

  precision_t a1, a2, a3, a4;
  for (int64_t j = 0; j < nY; j++) {
    for (int64_t i = 0; i < nX; i++) {

      a1 =calc_angle_given_three_lon_lat(lonCorner(i, j+1),
						    latCorner(i, j+1),
						    lonCorner(i, j),
						    latCorner(i, j),
						    lonCorner(i+1, j),
						    latCorner(i+1, j));
//      a1 = cPIdiv2 - calc_angle_given_three_lon_lat(lonCorner(0, nY-1),
//						    latCorner(0, nY-1),
//						    lonCorner(0, 0),
//						    latCorner(0, 0),
//						    lonCorner(nX-1, 0),
//						    latCorner(nX-1, 0));
      a2 = calc_angle_given_three_lon_lat(lonCorner(i, j),
					  latCorner(i, j),
					  lonCorner(i+1, j),
					  latCorner(i+1, j),
					  lonCorner(i+1, j+1),
					  latCorner(i+1, j+1));
      a3 = calc_angle_given_three_lon_lat(lonCorner(i+1, j),
					  latCorner(i+1, j),
					  lonCorner(i+1, j+1),
					  latCorner(i+1, j+1),
					  lonCorner(i, j+1),
					  latCorner(i, j+1));
      a4 = calc_angle_given_three_lon_lat(lonCorner(i+1, j+1),
					  latCorner(i+1, j+1),
					  lonCorner(i, j+1),
					  latCorner(i, j+1),
					  lonCorner(i, j),
					  latCorner(i, j));
      area(i,j) = 3 * 3 * (a1 + a2 + a3 + a4 - cTWOPI);
      area(i,j) = area(i,j) + 3 * 3 *
	(dxDU(i,j) + dxDU(i,j+1))/2 * (dyLR(i,j) + dyLR(i+1,j))/2;
      std::cout << "angles : " << j << " " << i << " "
		<< a1 << "  "
		<< a2 << "  "
		<< a3 << "  "
		<< a4 << " "
		<< dxDU(i,j) << " "
		<< dyLR(i,j) << " "
		<< 3*3 * dxDU(i,j) * dyLR(i,j) << " "
		<< area(i,j) << "\n";
      
    }
  }
	
  
  // state variables:
  if (verbose > 0) std::cout << "---> initializing rho\n";
  arma_mat rho = init_rho(x, y);

  if (verbose > 0) std::cout << "---> initializing vel\n";
  arma_mat xVel = init_vel(x, y) + 5.0;
  arma_mat yVel = init_vel(x, y);
  arma_mat vel = sqrt(xVel % xVel + yVel % yVel);
  arma_mat velL2, velR2, velU2, velD2;
  arma_mat xVelPrime, yVelPrime;

  // temp is "e" (not E):
  if (verbose > 0) std::cout << "---> initializing temp\n";
  arma_mat temp = init_temp(x, y);
  
  arma_mat eq1FluxLR, eq1FluxDU;
  arma_mat eq1FluxL, eq1FluxR, eq1FluxD, eq1FluxU;

  arma_mat eq2FluxLR, eq2FluxDU;
  arma_mat eq2FluxL, eq2FluxR, eq2FluxD, eq2FluxU;

  arma_mat eq3FluxLR, eq3FluxDU;
  arma_mat eq3FluxL, eq3FluxR, eq3FluxD, eq3FluxU;

  arma_mat eq4FluxLR, eq4FluxDU;
  arma_mat eq4FluxL, eq4FluxR, eq4FluxD, eq4FluxU;

  arma_mat wsL, wsR, wsD, wsU, wsLR, wsDU;

  if (verbose > 0) std::cout << "---> exchanging\n";
  exchange(rho, nGCs);
  exchange(xVel, nGCs);
  exchange(yVel, nGCs);
  exchange(temp, nGCs);
  
  if (verbose > 0) std::cout << "---> initializing momentum\n";
  arma_mat xMomentum = rho % xVel;
  arma_mat yMomentum = rho % yVel;
  arma_mat grad_xMomenum, xMomentumL, xMomentumR, xMomentumD, xMomentumU;
  arma_mat grad_yMomenum, yMomentumL, yMomentumR, yMomentumD, yMomentumU;

  if (verbose > 0) std::cout << "---> initializing totale\n";
  arma_mat totalE = rho % temp + 0.5 * rho % vel % vel;
  arma_mat grad_totalE, totaleL, totaleR, totaleD, totaleU;
  
  if (verbose > 0) std::cout << "---> outputting\n";
  
  output(x, "x.txt", false);
  output(y, "y.txt", false);
  output(jacobian, "jacobian.txt", false);
  output(g11, "g11.txt", false);
  output(g12, "g12.txt", false);
  output(g21, "g21.txt", false);
  output(g22, "g22.txt", false);
  output(rho, "rho.txt", false);
  output(xVel, "xVel.txt", false);
  output(yVel, "yVel.txt", false);
  output(temp, "temp.txt", false);
  output(totalE, "totale.txt", false);

  arma_mat diff;

  projection_struct rhoP;
  projection_struct xVelP, xVelPrimeP;
  projection_struct yVelP, yVelPrimeP;
  projection_struct tempP;
  projection_struct jacobianP;

  jacobianP = project_to_edges(jacobian, x, xEdges, y, yEdges, nGCs);
  

  iStep = 0;
  while (current_time < total_time) {

    if (verbose > 0)
      std::cout << "step : " << iStep << "; time : " << current_time << "\n";

    iStep++;

    // ------------------------------------------------
    // Calculate derived equations (at cell centers - these will be updated)
    // eq 1 = rho
    // eq 2 = rho * xVel (momentum)
    // eq 3 = rho * yVel (momentum)
    // eq 4 = E --> rho * (temp + 0.5 * vel^2) (totalE) 
    
    vel = sqrt(xVel % xVel + yVel % yVel);
    totalE = rho % temp + 0.5 * rho % vel % vel;
    xMomentum = rho % xVel;
    yMomentum = rho % yVel;

    // Get the velocities in the local coordinate system:
    xVelPrime = xVel % g11 + yVel % g12;
    yVelPrime = xVel % g21 + yVel % g22;
        
    // -----------------------------------
    // Project State Variables to Edges:

    if (verbose > 3) std::cout << "Projecting\n";
    
    rhoP = project_to_edges(rho, x, xEdges, y, yEdges, nGCs);
    xVelP = project_to_edges(xVel, x, xEdges, y, yEdges, nGCs);
    yVelP = project_to_edges(yVel, x, xEdges, y, yEdges, nGCs);
    tempP = project_to_edges(temp, x, xEdges, y, yEdges, nGCs);

    xVelPrimeP = project_to_edges(xVelPrime, x, xEdges, y, yEdges, nGCs);
    yVelPrimeP = project_to_edges(yVelPrime, x, xEdges, y, yEdges, nGCs);
    
    // ------------------------------------------------
    // Calculate derived equations (at edges)
    // eq 1 = rho
    // eq 2 = rho * xVel (momentum)
    // eq 3 = rho * yVel (momentum)
    // eq 4 = E --> rho * (temp + 0.5 * vel^2) (totalE) 
    
    if (verbose > 3) std::cout << "Deriving\n";

    xMomentumL = rhoP.L % xVelP.L;
    xMomentumR = rhoP.R % xVelP.R;
    xMomentumD = rhoP.D % xVelP.D;
    xMomentumU = rhoP.U % xVelP.U;
    
    yMomentumL = rhoP.L % yVelP.L;
    yMomentumR = rhoP.R % yVelP.R;
    yMomentumD = rhoP.D % yVelP.D;
    yMomentumU = rhoP.U % yVelP.U;

    velL2 = xVelP.L % xVelP.L + yVelP.L % yVelP.L;
    velR2 = xVelP.R % xVelP.R + yVelP.R % yVelP.R;
    velD2 = xVelP.D % xVelP.D + yVelP.D % yVelP.D;
    velU2 = xVelP.U % xVelP.U + yVelP.U % yVelP.U;
        
    totaleL = rhoP.L % tempP.L + 0.5 * rhoP.L % velL2;
    totaleR = rhoP.R % tempP.R + 0.5 * rhoP.R % velR2;
    totaleD = rhoP.D % tempP.D + 0.5 * rhoP.D % velD2;
    totaleU = rhoP.U % tempP.U + 0.5 * rhoP.U % velU2;

    // ------------------------------------------------
    // Calculate fluxes of different terms at the edges:
    
    if (verbose > 3) std::cout << "Calculating Fluxes\n";

    eq1FluxL = jacobianP.L % rhoP.L % xVelPrimeP.L;
    eq1FluxR = jacobianP.R % rhoP.R % xVelPrimeP.R;
    eq1FluxD = jacobianP.D % rhoP.D % yVelPrimeP.D;
    eq1FluxU = jacobianP.U % rhoP.U % yVelPrimeP.U;

    eq2FluxL = jacobianP.L % rhoP.L % (xVelP.L % xVelPrimeP.L + (gamma-1) * tempP.L);
    eq2FluxR = jacobianP.R % rhoP.R % (xVelP.R % xVelPrimeP.R + (gamma-1) * tempP.R);
    eq2FluxD = jacobianP.D % rhoP.D % xVelP.D % yVelPrimeP.D;
    eq2FluxU = jacobianP.U % rhoP.U % xVelP.U % yVelPrimeP.U;

    eq3FluxR = jacobianP.L % rhoP.R % xVelPrimeP.R % yVelP.R;
    eq3FluxL = jacobianP.R % rhoP.L % xVelPrimeP.L % yVelP.L;
    eq3FluxD = jacobianP.D % rhoP.D % (yVelPrimeP.D % yVelP.D + (gamma-1) * tempP.D);
    eq3FluxU = jacobianP.U % rhoP.U % (yVelPrimeP.U % yVelP.U + (gamma-1) * tempP.U);

    eq4FluxL = jacobianP.L % rhoP.L % xVelPrimeP.L % (0.5 * velL2 + gamma * tempP.L);
    eq4FluxR = jacobianP.R % rhoP.R % xVelPrimeP.R % (0.5 * velR2 + gamma * tempP.R);
    eq4FluxD = jacobianP.D % rhoP.D % yVelPrimeP.D % (0.5 * velD2 + gamma * tempP.D);
    eq4FluxU = jacobianP.U % rhoP.U % yVelPrimeP.U % (0.5 * velU2 + gamma * tempP.U);
			    
    // ------------------------------------------------
    // Calculate the wave speed for the diffusive flux:

    if (verbose > 3) std::cout << "Diffusive Fluxes\n";
    
    wsL = sqrt(velL2) + sqrt(gamma * (gamma-1) * tempP.L);
    wsR = sqrt(velR2) + sqrt(gamma * (gamma-1) * tempP.R);
    wsD = sqrt(velD2) + sqrt(gamma * (gamma-1) * tempP.D);
    wsU = sqrt(velU2) + sqrt(gamma * (gamma-1) * tempP.U);

    wsLR = wsR;
    for (int64_t i = 0; i < nX+1; i++) {
      for (int64_t j = 0; j < nY; j++) {
	if (wsL(i, j) > wsLR(i, j)) wsLR(i, j) = wsL(i, j);
      }
    }

    wsDU = wsD;
    for (int64_t i = 0; i < nX; i++) {
      for (int64_t j = 0; j < nY+1; j++) {
	if (wsU(i, j) > wsDU(i, j)) wsDU(i, j) = wsU(i, j);
      }
    }

    // ------------------------------------------------
    // Calculate dt based on max waves speeds and cell sizes
    
    if (verbose > 3) std::cout << "Calculating dt\n";

    dt = calc_dt(dxUD, dyLR, wsLR, wsDU, nGCs);
    dt = cfl * dt;
    current_time += dt;

    // ------------------------------------------------
    // Calculate average flux at the edges:

    if (verbose > 3) std::cout << "Averaging fluxes at edges\n";
    
    diff = rhoP.R - rhoP.L;
    eq1FluxLR = (eq1FluxL + eq1FluxR) / 2 + 0.5 * wsLR % diff;
    diff = rhoP.U - rhoP.D;
    eq1FluxDU = (eq1FluxD + eq1FluxU) / 2 + 0.5 * wsDU % diff;

    diff = xMomentumR - xMomentumL;
    eq2FluxLR = (eq2FluxL + eq2FluxR) / 2 + 0.5 * wsLR % diff;
    diff = xMomentumU - xMomentumD;
    eq2FluxDU = (eq2FluxD + eq2FluxU) / 2 + 0.5 * wsDU % diff;

    diff = yMomentumR - yMomentumL;
    eq3FluxLR = (eq3FluxL + eq3FluxR) / 2 + 0.5 * wsLR % diff;
    diff = yMomentumU - yMomentumD;
    eq3FluxDU = (eq3FluxD + eq3FluxU) / 2 + 0.5 * wsDU % diff;

    diff = totaleR - totaleL;
    eq4FluxLR = (eq4FluxL + eq4FluxR) / 2 + 0.5 * wsLR % diff;
    diff = totaleU - totaleD;
    eq4FluxDU = (eq4FluxD + eq4FluxU) / 2 + 0.5 * wsDU % diff;

    // ------------------------------------------------
    // Update values:
    if (verbose > 3) std::cout << "Updating equations of state\n";

    for (int64_t j = nGCs; j < nY + nGCs; j++) {
      for (int64_t i = nGCs; i < nX + nGCs; i++) {
	rho(i,j) = rho(i,j) - dt / area(i,j) *
	  (dyLR(i+1,j) * eq1FluxLR(i+1,j) -
	   dyLR(i,j) * eq1FluxLR(i,j) +
	   dxUD(i,j+1) * eq1FluxDU(i,j+1) -
	   dxUD(i,j) * eq1FluxDU(i,j))/jacobian(i,j);
	xMomentum(i,j) = xMomentum(i,j) - dt / area(i,j) *
	  (dyLR(i+1,j) * eq2FluxLR(i+1,j) -
	   dyLR(i,j) * eq2FluxLR(i,j) +
	   dxUD(i,j+1) * eq2FluxDU(i,j+1) -
	   dxUD(i,j) * eq2FluxDU(i,j))/jacobian(i,j);
	yMomentum(i,j) = yMomentum(i,j) - dt / area(i,j) *
	  (dyLR(i+1,j) * eq3FluxLR(i+1,j) -
	   dyLR(i,j) * eq3FluxLR(i,j) +
	   dxUD(i,j+1) * eq3FluxDU(i,j+1) -
	   dxUD(i,j) * eq3FluxDU(i,j))/jacobian(i,j);
	totalE(i,j) = totalE(i,j) - dt / area(i,j) *
	  (dyLR(i+1,j) * eq4FluxLR(i+1,j) -
	   dyLR(i,j) * eq4FluxLR(i,j) +
	   dxUD(i,j+1) * eq4FluxDU(i,j+1) -
	   dxUD(i,j) * eq4FluxDU(i,j))/jacobian(i,j);
      }
    }

    // ------------------------------------------------
    // Rederive state variables from the equations of state:

    if (verbose > 3) std::cout << "Deriving State Variables\n";

    //xVel = xMomentum / rho;
    //yVel = yMomentum / rho;
    //temp = totalE / rho - 0.5 * (xVel % xVel + yVel % yVel);
      
    // ------------------------------------------------
    // Exchange messages (set BCs, really):
    if (verbose > 3) std::cout << "Exchanging messages\n";

    exchange(rho, nGCs);
    exchange(xVel, nGCs);
    exchange(yVel, nGCs);
    exchange(temp, nGCs);

    if (verbose > 3) std::cout << "Outputing\n";

    if (int((current_time - dt)/dtOut) != int((current_time )/dtOut)) {
      std::cout << "Outputing at time : " << current_time << "\n";
      output(rho, "rho.txt", true);
      output(xVel, "xvel.txt", true);
      output(yVel, "yvel.txt", true);
      output(temp, "temp.txt", true);
      output(totalE, "totale_a.txt", true);
      //output(rhoL, "rhor.txt", false);
      //output(rhoR, "rhol.txt", false);
    }
  }   
  return 0;
}
