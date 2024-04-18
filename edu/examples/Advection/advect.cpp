
// g++ -I/usr/local/include -I/Users/ridley/Software/Json/json/include main.cpp

#include "../../../include/aether.h"
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

  precision_t dt = 0.01;

  int64_t nSteps = 3;
  int64_t iStep;
  
  int64_t nPts = 60;
  int64_t nGCs = 2;
  int64_t nPtsTotal = nGCs + nPts + nGCs;

  arma_vec x = init_grid(nPts, nGCs);
  arma_vec edges = calc_bin_edges(x);
  arma_vec widths = calc_bin_widths(edges);

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
  
  for (iStep = 0; iStep < nSteps; iStep++) {
  
    grad_rho = calc_grad(rho, x, nPts, nGCs);

    // Right side of edge from left
    rhoR = project_from_left(rho, grad_rho,
			     x, edges,
			     nPts, nGCs);
    //rhoR = limiter_value(rhoR, rho, nPts, nGCs);

    // Left side of edge from left
    rhoL = project_from_right(rho, grad_rho,
			      x, edges,
			      nPts, nGCs);
    //rhoL = limiter_value(rhoL, rho, nPts, nGCs);


    grad_vel = calc_grad(vel, x, nPts, nGCs);

    // Right side of edge from left
    velR = project_from_left(vel, grad_vel,
			     x, edges,
			     nPts, nGCs);
    //velR = limiter_value(velR, vel, nPts, nGCs);

    // Left side of edge from left
    velL = project_from_right(vel, grad_vel,
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
