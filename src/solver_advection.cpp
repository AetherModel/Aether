// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md


#include "aether.h"

// --------------------------------------------------------------------------
//
// --------------------------------------------------------------------------

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
//
// ---------------------------------------------------------

void advect(Grid &grid,
            Times &time,
            Neutrals &neutrals) {

  std::string function = "advect";
  static int iFunction = -1;
  report.enter(function, iFunction);

  projection_struct rhoP;
  projection_struct xVelP;
  projection_struct yVelP;
  projection_struct tempP;
  projection_struct gammaP;

  precision_t gamma = 5.0/3.0;
  precision_t dt = time.get_dt();

  int64_t nGCs = grid.get_nGCs();
  int64_t nX = grid.get_nX();
  int64_t nY = grid.get_nY();
  int64_t nAlts = grid.get_nAlts(), iAlt;

  arma_mat x, xEdges;
  arma_mat y, yEdges;

  arma_mat xVel, yVel, vel;
  arma_mat velL2, velR2, velU2, velD2;
  arma_mat rho, temp;
  arma_mat xMomentum, yMomentum;
  arma_mat grad_xMomenum, xMomentumL, xMomentumR, xMomentumD, xMomentumU;
  arma_mat grad_yMomenum, yMomentumL, yMomentumR, yMomentumD, yMomentumU;

  arma_mat eq1FluxLR, eq1FluxDU;
  arma_mat eq1FluxL, eq1FluxR, eq1FluxD, eq1FluxU;

  arma_mat eq2FluxLR, eq2FluxDU;
  arma_mat eq2FluxL, eq2FluxR, eq2FluxD, eq2FluxU, eq2Flux;

  arma_mat eq3FluxLR, eq3FluxDU;
  arma_mat eq3FluxL, eq3FluxR, eq3FluxD, eq3FluxU, eq3Flux;

  arma_mat eq4FluxLR, eq4FluxDU;
  arma_mat eq4FluxL, eq4FluxR, eq4FluxD, eq4FluxU;

  arma_mat wsL, wsR, wsD, wsU, wsLR, wsDU;

  arma_mat totalE;
  arma_mat grad_totalE, totaleL, totaleR, totaleD, totaleU;
  arma_mat diff;

  arma_mat area, xWidth, yWidth, geometry;

  arma_mat gamma2d;

  // These are all needed by the solver:
  neutrals.calc_mass_density();
  neutrals.calc_mean_major_mass();
  neutrals.calc_specific_heat();

  arma_mat t_to_e;

  for (iAlt = nGCs; iAlt < nAlts - nGCs; iAlt++) {

    if (report.test_verbose(3))
        std::cout << "Advection: Working with iAlt: " << iAlt << "\n";

    xVel = neutrals.velocity_vcgc[0].slice(iAlt);
    yVel = neutrals.velocity_vcgc[1].slice(iAlt);
    rho = neutrals.rho_scgc.slice(iAlt);
    // this is "e", or temperature expressed as an energy
    gamma2d = neutrals.gamma_scgc.slice(iAlt);
    t_to_e = 1.0 / (gamma2d - 1.0) * cKB / neutrals.mean_major_mass_scgc.slice(iAlt);
    temp = t_to_e % neutrals.temperature_scgc.slice(iAlt);

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

    x = grid.x_Center.slice(iAlt) * grid.radius_scgc(1,1,iAlt);
    y = grid.y_Center.slice(iAlt) * grid.radius_scgc(1,1,iAlt);
    xEdges = grid.x_Left.slice(iAlt) * grid.radius_scgc(1,1,iAlt);
    yEdges = grid.y_Down.slice(iAlt) * grid.radius_scgc(1,1,iAlt);

    rhoP = project_to_edges(rho, x, xEdges, y, yEdges, nGCs);
    xVelP = project_to_edges(xVel, x, xEdges, y, yEdges, nGCs);
    yVelP = project_to_edges(yVel, x, xEdges, y, yEdges, nGCs);
    tempP = project_to_edges(temp, x, xEdges, y, yEdges, nGCs);
    gammaP = project_to_edges(gamma2d, x, xEdges, y, yEdges, nGCs);

    // ------------------------------------------------
    // Calculate derived equations (at edges)
    // eq 1 = rho
    // eq 2 = rho * xVel (momentum)
    // eq 3 = rho * yVel (momentum)
    // eq 4 = E --> rho * (temp + 0.5 * vel^2) (totalE)

    report.print(3, "Advection: Deriving State Equations");

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

    report.print(3, "Advection: Calculating Fluxes");

    eq1FluxL = rhoP.L % xVelP.L;
    eq1FluxR = rhoP.R % xVelP.R;
    eq1FluxD = rhoP.D % yVelP.D;
    eq1FluxU = rhoP.U % yVelP.U;

    eq2FluxL = rhoP.L % (xVelP.L % xVelP.L + (gammaP.L - 1) % tempP.L);
    eq2FluxR = rhoP.R % (xVelP.R % xVelP.R + (gammaP.R - 1) % tempP.R);
    eq2FluxD = rhoP.D % xVelP.D % yVelP.D;
    eq2FluxU = rhoP.U % xVelP.U % yVelP.U;
    eq2Flux = rho % xVel % yVel;

    eq3FluxR = rhoP.R % xVelP.R % yVelP.R;
    eq3FluxL = rhoP.L % xVelP.L % yVelP.L;
    eq3FluxD = rhoP.D % (yVelP.D % yVelP.D + (gammaP.D - 1) % tempP.D);
    eq3FluxU = rhoP.U % (yVelP.U % yVelP.U + (gammaP.U - 1) % tempP.U);
    eq3Flux = rho % (yVel % yVel + (gamma2d - 1) % temp);

    eq4FluxL = rhoP.L % xVelP.L % (0.5 * velL2 + gammaP.L % tempP.L);
    eq4FluxR = rhoP.R % xVelP.R % (0.5 * velR2 + gammaP.R % tempP.R);
    eq4FluxD = rhoP.D % yVelP.D % (0.5 * velD2 + gammaP.D % tempP.D);
    eq4FluxU = rhoP.U % yVelP.U % (0.5 * velU2 + gammaP.U % tempP.U);

    // ------------------------------------------------
    // Calculate the wave speed for the diffusive flux:

    report.print(3, "Advection: Diffusive Fluxes");

    wsL = sqrt(velL2) + sqrt(gammaP.L % (gammaP.L - 1) % tempP.L);
    wsR = sqrt(velR2) + sqrt(gammaP.R % (gammaP.R - 1) % tempP.R);
    wsD = sqrt(velD2) + sqrt(gammaP.D % (gammaP.D - 1) % tempP.D);
    wsU = sqrt(velU2) + sqrt(gammaP.U % (gammaP.U - 1) % tempP.U);

    wsLR = wsR;
    for (int64_t i = 0; i < nX + 1; i++) {
      for (int64_t j = 0; j < nY; j++) {
	       if (wsL(i, j) > wsLR(i, j)) wsLR(i, j) = wsL(i, j);
      }
    }

    wsDU = wsD;
    for (int64_t i = 0; i < nX; i++) {
      for (int64_t j = 0; j < nY + 1; j++) {
	       if (wsU(i, j) > wsDU(i, j)) wsDU(i, j) = wsU(i, j);
      }
    }

    // ------------------------------------------------
    // Calculate average flux at the edges:

    report.print(3, "Advection: Averaging fluxes at edges");

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
    report.print(3, "Advection: Updating equations of state");

    area = grid.cell_area.slice(iAlt) * grid.radius2_scgc(1,1,iAlt);
    yWidth = grid.dy_Left.slice(iAlt) * grid.radius_scgc(1,1,iAlt);
    xWidth = grid.dx_Down.slice(iAlt) * grid.radius_scgc(1,1,iAlt);
    
    geometry = 
      sin(grid.geoLat_scgc.slice(iAlt)) /
      cos(grid.geoLat_scgc.slice(iAlt)) /
      grid.radius_scgc(1,1,iAlt);
    
    for (int64_t j = nGCs; j < nY - nGCs; j++) {
      for (int64_t i = nGCs; i < nX - nGCs; i++) {
	//if (i == nGCs) cout << "j = " << j << " " << xWidth(i,j) << "\n";
	rho(i,j) = rho(i,j) - dt *
	  (yWidth(i+1,j) * eq1FluxLR(i+1,j) -
	   yWidth(i,j) * eq1FluxLR(i,j) +
	   xWidth(i,j+1) * eq1FluxDU(i,j+1) -
	   xWidth(i,j) * eq1FluxDU(i,j)) / area(i,j);
	xMomentum(i,j) = xMomentum(i,j) - dt *
	  ((yWidth(i+1,j) * eq2FluxLR(i+1,j) -
	    yWidth(i,j) * eq2FluxLR(i,j) +
	    xWidth(i,j+1) * eq2FluxDU(i,j+1) -
	    xWidth(i,j) * eq2FluxDU(i,j)) / area(i,j) -
	   geometry(i,j) * eq2Flux(i,j));
	yMomentum(i,j) = yMomentum(i,j) - dt *
	  ((yWidth(i+1,j) * eq3FluxLR(i+1,j) -
	    yWidth(i,j) * eq3FluxLR(i,j) +
	    xWidth(i,j+1) * eq3FluxDU(i,j+1) -
	    xWidth(i,j) * eq3FluxDU(i,j)) / area(i,j) +
	   geometry(i,j) * eq3Flux(i,j));
	totalE(i,j) = totalE(i,j) - dt *
	  (yWidth(i+1,j) * eq4FluxLR(i+1,j) -
	   yWidth(i,j) * eq4FluxLR(i,j) +
	   xWidth(i,j+1) * eq4FluxDU(i,j+1) -
	   xWidth(i,j) * eq4FluxDU(i,j)) / area(i,j);
      }
    }

    xVel = xMomentum / rho;
    yVel = yMomentum / rho;

    neutrals.velocity_vcgc[0].slice(iAlt) = xVel;
    neutrals.velocity_vcgc[1].slice(iAlt) = yVel;
    temp = (totalE / rho - 0.5 * (xVel % xVel + yVel % yVel)) / t_to_e;
    precision_t fac, dm, dp;
    
    for (int64_t j = nGCs; j < nY - nGCs; j++) {
      for (int64_t i = nGCs; i < nX - nGCs; i++) {
	fac = 1.0;
	if (cos(grid.geoLat_scgc(i,j,iAlt)) < 0.2) {
	  fac = fac * (0.2 - cos(grid.geoLat_scgc(i,j,iAlt)));
	}
	dm = (1.0 - fac) * neutrals.temperature_scgc(i,j,iAlt);
	dp = (1.0 + fac) * neutrals.temperature_scgc(i,j,iAlt);
	if (temp(i,j) < dm) temp(i,j) = dm;
	if (temp(i,j) > dp) temp(i,j) = dp;
	neutrals.temperature_scgc(i,j,iAlt) = temp(i,j);

	dm = (1.0 - fac) * neutrals.rho_scgc(i,j,iAlt);
	dp = (1.0 + fac) * neutrals.rho_scgc(i,j,iAlt);
	if (rho(i,j) < dm) rho(i,j) = dm;
	if (rho(i,j) > dp) rho(i,j) = dp;
	neutrals.rho_scgc(i,j,iAlt) = rho(i,j);
      }
    }
    if (report.test_verbose(3) && iAlt == 8) {
      std::cout << "end t : " << neutrals.temperature_scgc.slice(iAlt).min() << " " << neutrals.temperature_scgc.slice(iAlt).max() << "\n";
      std::cout << "end temp : " << temp.min() << " " << temp.max() << "\n";
      std::cout << "end xVel : " << xVel.min() << " " << xVel.max() << "\n";
      std::cout << "end yVel : " << yVel.min() << " " << yVel.max() << "\n";
    }
  }
  neutrals.calc_density_from_mass_concentration();

  report.exit(function);
  return;
}
