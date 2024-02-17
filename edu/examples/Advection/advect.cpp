
// g++ -I/usr/local/include -I/Users/ridley/Software/Json/json/include main.cpp
// g++ -I/usr/local/include -o advect1d advect.cpp

/// The armadillo library is to allow the use of 3d cubes and other
/// array types, with array math built in. This eliminates loops!
#include <armadillo>

/// This is used for timing and the random seed generator:
#include <chrono>

// Types
// Precision compile-time aliasing
#ifdef AETHER_USE_PRECISION_DOUBLE
/// Precision type chosen to be `double` through `AETHER_USE_PRECISION_DOUBLE`
using precision_t = double;
#else
/// Precision type compile-time default to float.
using precision_t = float;
#endif

/// Armadillo type vector (single column) with compile-time precision.
using arma_vec = arma::Col<precision_t>;
/// Armadillo type matrix (two dimension) with compile-time precision.
using arma_mat = arma::Mat<precision_t>;
/// Armadillo type cube (three dimension) with compile-time precision.
using arma_cube = arma::Cube<precision_t>;

#include <fstream>

// ---------------------------------------------------------
// grid creation
// ---------------------------------------------------------

arma_vec init_grid(int64_t nPts, int64_t nGCs) {

  precision_t dx = 1.0 / nPts;
  arma_vec x(nPts + nGCs * 2);

  // uniform grid:
  for (int64_t i = -nGCs; i < nPts + nGCs; i++) {
    x(i + nGCs) = i * dx;
  }
  precision_t maxX = x(nPts + nGCs - 1);
  x = 100.0 * x / maxX;

  return x;
}

// ---------------------------------------------------------
// grid stretched creation
// ---------------------------------------------------------

arma_vec init_stretched_grid(int64_t nPts, int64_t nGCs) {

  precision_t dx = 1.0;
  arma_vec x(nPts + nGCs * 2);

  precision_t factor = 1.0;
  precision_t i2pi = 2.0 * 3.1415927 / (nPts-1);

  x(nGCs) = 0.0;
  
  for (int64_t i = 1; i < nPts + nGCs; i++) {
    x(i + nGCs) = x(i - 1 + nGCs) + dx  + factor * (1 + cos(i * i2pi));
    std::cout << "i : " << i << " " << cos(i * i2pi) << "\n";
  }
  for (int64_t i = -1; i >= -nGCs; i--) {
    x(i + nGCs) = x(i + 1 + nGCs) - dx  - factor * (1 + cos(i * i2pi));
    std::cout << "i : " << i << " " << cos(i * i2pi) << "\n";
  }
  precision_t maxX = x(nPts + nGCs - 1);
  x = 100.0 * x / maxX;
  
  return x;
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
// initial density
// ---------------------------------------------------------

arma_vec init_den(int64_t nPts) {

  arma_vec den(nPts);
  for (int64_t i = 0; i < nPts; i++) {
    if (i < nPts/2)
      den(i) = 2.0;
    else
      den(i) = 1.0;
  }

  return den;
}

// ---------------------------------------------------------
// initial velocity
// ---------------------------------------------------------

arma_vec init_vel(int64_t nPts) {
  arma_vec vel(nPts);
  // all cells positive to right:
  vel.ones();
  vel = -1.0 * vel;
  return vel;
}

// ---------------------------------------------------------
// exchange messages
// ---------------------------------------------------------

void exchange(arma_vec &values, int64_t nPts, int64_t nGCs) {

  int64_t iEnd = nPts + 2 * nGCs;
  // this is a periodic BC:
  for (int64_t i = 0; i < nGCs; i++) {
    values(i) = values(iEnd - 2 * nGCs + i);
    values(iEnd - nGCs + i) = values(nGCs + i);
  }
}

// ---------------------------------------------------------
//
// ---------------------------------------------------------

arma_vec limiter_mc(arma_vec left,
		    arma_vec right,
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

// ---------------------------------------------------------
// calc gradients at centers
//   - values and x defined at centers
// ---------------------------------------------------------

arma_vec calc_grad(arma_vec values,
		   arma_vec x,
		   int64_t nPts,
		   int64_t nGCs) {

  arma_vec gradients = values * 0.0;
  arma_vec gradL = values * 0.0;
  arma_vec gradR = values * 0.0;

  precision_t factor1 = 0.625;
  precision_t factor2 = 0.0416667;
  precision_t h;

  int64_t i;
  
  i = nGCs - 1;
  h = 2.0 / (x(i+1) - x(i));
  gradR(i) = h * (factor1 * (values(i+1) - values(i)) -
		  factor2 * (values(i+2) - values(i-1)));
  gradL(i) = (values(i) - values(i-1)) / (x(i) - x(i-1));
  
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
// Project gradients + values to the right face, from the left
//   returned values are on the i - 1/2 edges
//     (between i-1 and i cell center)
// ---------------------------------------------------------

arma_vec project_from_left(arma_vec values,
			   arma_vec gradients,
			   arma_vec x_centers,
			   arma_vec x_edges,
			   int64_t nPts,
			   int64_t nGCs) {
  int64_t iStart = 0;
  int64_t iEnd = nPts + 2 * nGCs;

  // Define at edges:
  arma_vec projected(nPts + 2 * nGCs + 1);
  projected.zeros();

  // no gradient in the 0 or iEnd cells
  for (int64_t i = iStart + 1; i < iEnd - 1; i++)
    projected(i + 1) = values(i) +
      gradients(i) * (x_edges(i + 1) - x_centers(i));

  return projected;
}

// ---------------------------------------------------------
// Project gradients + values to the right face, from the left
//   returned values are on the i - 1/2 edges
//     (between i-1 and i cell center)
// ---------------------------------------------------------

arma_vec project_from_left_new(arma_vec values,
			       arma_vec x_centers,
			       arma_vec x_edges,
			       int64_t nPts,
			       int64_t nGCs) {
  int64_t iStart = 1;
  int64_t iEnd = nPts + 2 * nGCs - 1;

  // Define at edges:
  arma_vec projected(nPts + 2 * nGCs + 1);
  projected.zeros();

  precision_t dxei, dxci, dxcip1, r;
  
  // no gradient in the 0 or iEnd cells
  for (int64_t i = iStart; i < iEnd; i++) {
    dxei = x_edges(i + 1) - x_edges(i);
    dxci = x_centers(i) - x_centers(i - 1);
    dxcip1 = x_centers(i + 1) - x_centers(i);
    r = dxcip1 / dxci;
    projected(i + 1) = values(i) + 
      0.5 * dxei * (values(i) - values(i - 1)) / dxci + 
      0.125 * dxei * dxei * (values(i + 1) + r * values(i - 1) - (1 + r) * values(i)) / (dxci * dxcip1);
  }

  projected = limiter_value(projected, values, nPts, nGCs);
  
  return projected;
}


// ---------------------------------------------------------
// Project gradients + values to the left face, from the right
//   returned values are on the i - 1 edges
//     (between i-1 and i cell center)
// ---------------------------------------------------------

arma_vec project_from_right(arma_vec values,
			    arma_vec gradients,
			    arma_vec x_centers,
			    arma_vec x_edges,
			    int64_t nPts,
			    int64_t nGCs) {
  int64_t iStart = 0;
  int64_t iEnd = nPts + 2 * nGCs;

  // Define at edges:
  arma_vec projected(nPts + 2 * nGCs + 1);
  projected.zeros();

  // no gradient in the 0 or iEnd cells
  for (int64_t i = iStart + 1; i < iEnd - 1; i++)
    projected(i) = values(i) +
      gradients(i) * (x_edges(i) - x_centers(i));

  return projected;
}

// ---------------------------------------------------------
// Project gradients + values to the left face, from the right
//   returned values are on the i - 1 edges
//     (between i-1 and i cell center)
// ---------------------------------------------------------

arma_vec project_from_right_new(arma_vec values,
				arma_vec x_centers,
				arma_vec x_edges,
				int64_t nPts,
				int64_t nGCs) {
  int64_t iStart = 1;
  int64_t iEnd = nPts + 2 * nGCs - 1;

  // Define at edges:
  arma_vec projected(nPts + 2 * nGCs + 1);
  precision_t dxei, dxci, dxcip1, r;
  
  projected.zeros();

  // no gradient in the 0 or iEnd cells
  for (int64_t i = iStart; i < iEnd; i++) {
    dxei = x_edges(i + 1) - x_edges(i);
    dxci = x_centers(i) - x_centers(i - 1);
    dxcip1 = x_centers(i + 1) - x_centers(i);
    r = dxcip1 / dxci;
    projected(i) = values(i) - 
      0.5 * dxei * (values(i + 1) - values(i)) / dxcip1 +
      0.125 * dxei * dxei * (values(i + 1) + r * values(i - 1) - (1 + r) * values(i)) / (dxci * dxcip1);
  }

  projected = limiter_value(projected, values, nPts, nGCs);
  
  return projected;
}

// ---------------------------------------------------------
// gudonov upwind scheme
// ---------------------------------------------------------

arma_vec gudonov(arma_vec valL,
		 arma_vec valR,
		 arma_vec velL,
		 arma_vec velR,
		 int64_t nPts,
		 int64_t nGCs) {
  
  int64_t iStart = 0;
  int64_t iEnd = nPts + 2 * nGCs;

  arma_vec flux = velL * 0.0;
  arma_vec vel = (velL + velR)/2.0;
  
  for (int64_t i = iStart + 1; i < iEnd - 1; i++) {
    if (vel(i) > 0)
      flux(i) = valR(i) * vel(i);
    else
      flux(i) = valL(i) * vel(i);
  }
  return flux;
}

// ---------------------------------------------------------
// gudonov upwind scheme
// ---------------------------------------------------------

arma_vec rusanov(arma_vec valL,
		 arma_vec valR,
		 arma_vec velL,
		 arma_vec velR,
		 arma_vec widths,
		 int64_t nPts,
		 int64_t nGCs) {
  
  int64_t iStart = 0;
  int64_t iEnd = nPts + 2 * nGCs;

  arma_vec ws = abs((velL + velR)/2.0) + 1.;
  arma_vec fluxL = valL % velL;
  arma_vec fluxR = valR % velR;
  arma_vec valDiff = valL - valR;
  arma_vec flux = (fluxL + fluxR) / 2.0;
  for (int64_t i = iStart + 1; i < iEnd - 1; i++)
    flux(i) = flux(i) - ws(i)/2 * valDiff(i);
    
  return flux;
}

// ---------------------------------------------------------
//
// ---------------------------------------------------------

void output(arma_vec values,
	    std::string filename,
	    bool DoAppend,
	    int64_t nPts,
	    int64_t nGCs) {

  std::ofstream outfile;
  if (DoAppend)
    outfile.open(filename, std::ios_base::app);
  else
    outfile.open(filename);
  
  int64_t i;
  for (i = 0; i < nPts + 2 * nGCs ; i++) {
    outfile << values(i) << " ";
  }
  outfile << "\n";
  outfile.close();
    
}


// ---------------------------------------------------------
//
// ---------------------------------------------------------


// ---------------------------------------------------------
// main code
// ---------------------------------------------------------

int main() {

  int64_t nPts = 200;
  int64_t nGCs = 2;
  int64_t nPtsTotal = nGCs + nPts + nGCs;

  arma_vec x = init_grid(nPts, nGCs);
  //arma_vec x = init_stretched_grid(nPts, nGCs);
  arma_vec edges = calc_bin_edges(x);
  arma_vec widths = calc_bin_widths(edges);

  precision_t dt = 0.1 * x(nPts + nGCs - 1) / nPts;
  precision_t time = 0.0;

  int64_t nSteps = x(nPts + nGCs - 1) / dt;
  int64_t iStep;
    
  arma_vec grad_rho;
  arma_vec rhoL;
  arma_vec rhoR;
  arma_vec grad_vel;
  arma_vec velL;
  arma_vec velR;

  arma_vec flux;

  arma_vec rho = init_den(nPtsTotal);
  arma_vec vel = init_vel(nPtsTotal);

  exchange(rho, nPts, nGCs);
  exchange(vel, nPts, nGCs);

  output(rho, "rho.txt", false, nPts, nGCs);
  output(x, "x.txt", false, nPts, nGCs);
  
  for (iStep = 0; iStep < nSteps; iStep++) {

    std::cout << "iStep = " << iStep << "; time =  " << time << "\n";
    time = time + dt;
  
    grad_rho = calc_grad(rho, x, nPts, nGCs);

    // Right side of edge from left
//rhoR = project_from_left(rho, grad_rho,
//			     x, edges,
//			     nPts, nGCs);
    rhoR = project_from_left_new(rho,
			     x, edges,
			     nPts, nGCs);
    //rhoR = limiter_value(rhoR, rho, nPts, nGCs);

    // Left side of edge from left
//    rhoL = project_from_right(rho, grad_rho,
//			      x, edges,
//			      nPts, nGCs);
    rhoL = project_from_right_new(rho,
			      x, edges,
			      nPts, nGCs);
    //rhoL = limiter_value(rhoL, rho, nPts, nGCs);


    grad_vel = calc_grad(vel, x, nPts, nGCs);

    // Right side of edge from left
//    velR = project_from_left(vel, grad_vel,
//			     x, edges,
//			     nPts, nGCs);
    velR = project_from_left_new(vel,
			     x, edges,
			     nPts, nGCs);
    //velR = limiter_value(velR, vel, nPts, nGCs);

    // Left side of edge from left
//    velL = project_from_right(vel, grad_vel,
//			      x, edges,
//			      nPts, nGCs);
    velL = project_from_right_new(vel,
			      x, edges,
			      nPts, nGCs);
    //velL = limiter_value(velL, vel, nPts, nGCs);


    flux = gudonov(rhoL, rhoR, velL, velR, nPts, nGCs);
    //flux = rusanov(rhoL, rhoR, velL, velR, widths, nPts, nGCs);

    for (int64_t i = nGCs; i < nPts + nGCs*2 - 1; i++) {
      rho(i) = rho(i) - dt / widths(i) *
	(flux(i+1) - flux(i));
    }
    exchange(rho, nPts, nGCs);
    exchange(vel, nPts, nGCs);

    output(rho, "rho.txt", true, nPts, nGCs);
    //output(flux, "flux.txt", false, nPts, nGCs);
    //output(rho, "test.txt", true, nPts, nGCs);
    //output(rhoL, "rhor.txt", false, nPts, nGCs);
    //output(rhoR, "rhol.txt", false, nPts, nGCs);
  }
  
  return 0;
}
