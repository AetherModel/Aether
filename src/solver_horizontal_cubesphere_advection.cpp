// Copyright 2023, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

// Initial version: F. Cheng, July 2023

#include "aether.h"

using namespace Cubesphere_tools;

// DOES NOT WORK WELL
std::vector<arma_mat> Neutrals::residual_horizontal_hlle_advection(
  std::vector<arma_mat>& states, Grid& grid, Times& time) {
  // Dimensions of Spatial Discretization
  int64_t nXs = grid.get_nX();
  int64_t nYs = grid.get_nY();
  int64_t nGCs = grid.get_nGCs();
  int64_t nAlts = grid.get_nAlts();
  int iAlt, iSpec;

  /** Extract Grid Features **/
  arma_mat x = grid.refx_scgc.slice(iAlt);
  arma_mat xEdges = grid.refx_Left.slice(iAlt);
  arma_mat y = grid.refy_scgc.slice(iAlt);
  arma_mat yEdges = grid.refy_Down.slice(iAlt);

  // Get reference grid dimensions (Assume dx = dy and equidistant)
  arma_vec x_vec = x.col(0);
  precision_t dx = x_vec(1) - x_vec(0);
  precision_t area = dx * dx;
  arma_mat jacobian = grid.sqrt_g_scgc.slice(iAlt);

  /** State/Velocity extraction **/
  /* MASS DENSITY */
  arma_mat rho = states[0];

  /* VELOCITY */
  // Convert to contravariant (reference) velocity
  arma_mat xVel = states[1]; // u^1
  arma_mat yVel = states[2]; // u^2

  // Generate velocity magnitude squared
  arma_mat vel2 = xVel % xVel + yVel % yVel;

  /** Advancing **/
  /* Initialize projection constructs and storages */
  projection_struct rhoP;
  projection_struct xVelP;
  projection_struct yVelP;

  // They are all pure scalar fields without sqrt(g)
  arma_mat rhoL, rhoR, rhoD, rhoU;
  arma_mat xVelL, xVelR, xVelD, xVelU;
  arma_mat yVelL, yVelR, yVelD, yVelU;

  arma_mat velL2, velR2, velD2, velU2;

  /** Initialize Flux and Wave Speed Storages */
  arma_mat eq1FluxLR_left, eq1FluxDU_down;
  arma_mat eq1FluxLR_right, eq1FluxDU_upper;
  arma_mat eq1FluxL, eq1FluxR, eq1FluxD, eq1FluxU;

  arma_mat wsL, wsR, wsD, wsU;
  arma_mat wsL_min, wsL_max, wsR_min, wsR_max;
  arma_mat wsD_min, wsD_max, wsU_min, wsU_max;
  arma_mat wsLR_max, wsDU_max, wsLR_min, wsDU_min;

  arma_mat diff; // for Riemann Solver

  /* Projection */
  rhoP = project_to_edges(rho, x, xEdges, y, yEdges, nGCs);
  xVelP = project_to_edges(xVel, x, xEdges, y, yEdges, nGCs);
  yVelP = project_to_edges(yVel, x, xEdges, y, yEdges, nGCs);

  // Resolve Scalar Fields into rho, xVel, yVel, and totalE (without rho)
  rhoL = rhoP.L;
  rhoR = rhoP.R;
  rhoD = rhoP.D;
  rhoU = rhoP.U;

  xVelL = xVelP.L;
  xVelR = xVelP.R;
  xVelD = xVelP.D;
  xVelU = xVelP.U;

  yVelL = yVelP.L;
  yVelR = yVelP.R;
  yVelD = yVelP.D;
  yVelU = yVelP.U;

  //velL2 = xVelL % xVelL + yVelL % yVelL;
  //velR2 = xVelR % xVelR + yVelR % yVelR;
  //velD2 = xVelD % xVelD + yVelD % yVelD;
  //velU2 = xVelU % xVelU + yVelU % yVelU;

  /* Calculate Edge Fluxes */
  // Note that dot product between normal vector at edge and flux vector
  // resolves into a pure one component flux or either hat{x} or hat{y}

  // Flux calculated from the left of the edge
  eq1FluxL = rhoL % xVelL % grid.sqrt_g_Left.slice(iAlt);
  // Flux calculated from the right of the edge
  eq1FluxR = rhoR % xVelR % grid.sqrt_g_Left.slice(iAlt);
  // Flux calculated from the down of the edge
  eq1FluxD = rhoD % yVelD % grid.sqrt_g_Down.slice(iAlt);
  // Flux calculated from the up of the edge
  eq1FluxU = rhoU % yVelU % grid.sqrt_g_Down.slice(iAlt);

  /* Wave Speed Calculation (Left/Down) */
  wsL = xVelL;
  wsR = xVelR;
  wsD = yVelD;
  wsU = yVelU;

  wsL_max = wsL;
  wsL_min = wsL;
  wsR_max = wsR;
  wsR_min = wsR;
  wsD_max = wsD;
  wsD_min = wsD;
  wsU_max = wsU;
  wsU_min = wsU;

  // Process wave speeds from each direction first
  for (int i = 0; i < nXs + 1; i++) {
    for (int j = 0; j < nYs; j++) {

      if (wsL(i, j) > 0.)
        wsL_min(i, j) = 0.;

      else
        wsL_max(i, j) = 0.;

      if (wsR(i, j) > 0.)
        wsR_min(i, j) = 0.;

      else
        wsR_max(i, j) = 0.;
    }
  }

  for (int i = 0; i < nXs; i++) {
    for (int j = 0; j < nYs + 1; j++) {
      if (wsD(i, j) > 0.)
        wsD_min(i, j) = 0.;

      else
        wsD_max(i, j) = 0.;

      if (wsU(i, j) > 0.)
        wsU_min(i, j) = 0.;

      else
        wsU_max(i, j) = 0.;
    }
  }

  // Process edge wave speeds
  wsLR_max = wsR_max;

  for (int i = 0; i < nXs + 1; i++) {
    for (int j = 0; j < nYs; j++) {
      if (wsL_max(i, j) > wsLR_max(i, j))
        wsLR_max(i, j) = wsL_max(i, j);
    }
  }

  wsDU_max = wsD_max;

  for (int i = 0; i < nXs; i++) {
    for (int j = 0; j < nYs + 1; j++) {
      if (wsU_max(i, j) > wsDU_max(i, j))
        wsDU_max(i, j) = wsU_max(i, j);
    }
  }

  wsLR_min = wsR_min;

  for (int i = 0; i < nXs + 1; i++) {
    for (int j = 0; j < nYs; j++) {
      if (wsL_min(i, j) < wsLR_min(i, j))
        wsLR_min(i, j) = wsL_min(i, j);
    }
  }

  wsDU_min = wsD_min;

  for (int i = 0; i < nXs; i++) {
    for (int j = 0; j < nYs + 1; j++) {
      if (wsU_min(i, j) < wsDU_min(i, j))
        wsDU_min(i, j) = wsU_min(i, j);
    }
  }

  /* Calculate average flux at the edges (HLLE Flux) */
  arma_mat wsLR_sum = wsLR_max + wsLR_min;
  arma_mat wsLR_diff = wsLR_max - wsLR_min;
  diff = (rhoR - rhoL) % grid.sqrt_g_Left.slice(
           iAlt); // State difference, need to add sqrt(g)
  eq1FluxLR_left = 0.5 * (eq1FluxL + eq1FluxR) + 0.5 * (wsLR_sum / wsLR_diff) %
                   (eq1FluxR - eq1FluxL) - (wsLR_max % wsLR_min) / wsLR_diff % diff;

  arma_mat wsDU_sum = wsDU_max + wsDU_min;
  arma_mat wsDU_diff = wsDU_max - wsDU_min;
  diff = (rhoU - rhoD) % grid.sqrt_g_Down.slice(iAlt);
  eq1FluxDU_down = 0.5 * (eq1FluxU + eq1FluxD) + 0.5 * (wsDU_sum / wsDU_diff) %
                   (eq1FluxD - eq1FluxU) - (wsDU_max % wsDU_min) / wsDU_diff % diff;

  /* Wave Speed Calculation (Right/Up) */
  wsL = -xVelL;
  wsR = -xVelR;
  wsD = -yVelD;
  wsU = -yVelU;

  wsL_max = wsL;
  wsL_min = wsL;
  wsR_max = wsR;
  wsR_min = wsR;
  wsD_max = wsD;
  wsD_min = wsD;
  wsU_max = wsU;
  wsU_min = wsU;

  // Process wave speeds from each direction first
  for (int i = 0; i < nXs + 1; i++) {
    for (int j = 0; j < nYs; j++) {

      if (wsL(i, j) > 0.)
        wsL_min(i, j) = 0.;

      else
        wsL_max(i, j) = 0.;

      if (wsR(i, j) > 0.)
        wsR_min(i, j) = 0.;

      else
        wsR_max(i, j) = 0.;
    }
  }

  for (int i = 0; i < nXs; i++) {
    for (int j = 0; j < nYs + 1; j++) {
      if (wsD(i, j) > 0.)
        wsD_min(i, j) = 0.;

      else
        wsD_max(i, j) = 0.;

      if (wsU(i, j) > 0.)
        wsU_min(i, j) = 0.;

      else
        wsU_max(i, j) = 0.;
    }
  }

  // Process edge wave speeds
  wsLR_max = wsR_max;

  for (int i = 0; i < nXs + 1; i++) {
    for (int j = 0; j < nYs; j++) {
      if (wsL_max(i, j) > wsLR_max(i, j))
        wsLR_max(i, j) = wsL_max(i, j);
    }
  }

  wsDU_max = wsD_max;

  for (int i = 0; i < nXs; i++) {
    for (int j = 0; j < nYs + 1; j++) {
      if (wsU_max(i, j) > wsDU_max(i, j))
        wsDU_max(i, j) = wsU_max(i, j);
    }
  }

  wsLR_min = wsR_min;

  for (int i = 0; i < nXs + 1; i++) {
    for (int j = 0; j < nYs; j++) {
      if (wsL_min(i, j) < wsLR_min(i, j))
        wsLR_min(i, j) = wsL_min(i, j);
    }
  }

  wsDU_min = wsD_min;

  for (int i = 0; i < nXs; i++) {
    for (int j = 0; j < nYs + 1; j++) {
      if (wsU_min(i, j) < wsDU_min(i, j))
        wsDU_min(i, j) = wsU_min(i, j);
    }
  }

  /* Calculate average flux at the edges (HLLE Flux) */
  wsLR_sum = wsLR_max + wsLR_min;
  wsLR_diff = wsLR_max - wsLR_min;
  diff = (rhoR - rhoL) % grid.sqrt_g_Left.slice(
           iAlt); // State difference, need to add sqrt(g)
  eq1FluxLR_right = 0.5 * (eq1FluxL + eq1FluxR) + 0.5 * (wsLR_sum / wsLR_diff) %
                    (eq1FluxR - eq1FluxL) - (wsLR_max % wsLR_min) / wsLR_diff % diff;

  wsDU_sum = wsDU_max + wsDU_min;
  wsDU_diff = wsDU_max - wsDU_min;
  diff = (rhoU - rhoD) % grid.sqrt_g_Down.slice(iAlt);
  eq1FluxDU_upper = 0.5 * (eq1FluxU + eq1FluxD) + 0.5 * (wsDU_sum / wsDU_diff) %
                    (eq1FluxD - eq1FluxU) - (wsDU_max % wsDU_min) / wsDU_diff % diff;


  // Setup residual storage for return
  arma_mat eq1_residual(nXs, nYs, fill::zeros);

  // State Update
  // Note the ghost cells WILL NOT BE UPDATED
  for (int j = nGCs; j < nYs - nGCs; j++) {
    for (int i = nGCs; i < nXs - nGCs; i++) {
      precision_t rhoResidual_ij = dx * eq1FluxLR_right(i + 1, j) -
                                   dx * eq1FluxLR_left(i, j) +
                                   dx * eq1FluxDU_upper(i, j + 1) -
                                   dx * eq1FluxDU_down(i, j);
      eq1_residual(i, j) = -1 / area * rhoResidual_ij;
    }
  }

  // Setup return vector
  std::vector<arma_mat> return_vector;
  return_vector.push_back(eq1_residual);

  return return_vector;
}

// WORKS, but diffusive
std::vector<arma_mat> Neutrals::residual_horizontal_rusanov_advection(
  std::vector<arma_mat>& states, Grid& grid, Times& time) {
  // Dimensions of Spatial Discretization
  int64_t nXs = grid.get_nX();
  int64_t nYs = grid.get_nY();
  int64_t nGCs = grid.get_nGCs();
  int64_t nAlts = grid.get_nAlts();
  int iAlt, iSpec;

  /** Extract Grid Features **/
  arma_mat x = grid.refx_scgc.slice(iAlt);
  arma_mat xEdges = grid.refx_Left.slice(iAlt);
  arma_mat y = grid.refy_scgc.slice(iAlt);
  arma_mat yEdges = grid.refy_Down.slice(iAlt);

  // Get reference grid dimensions (Assume dx = dy and equidistant)
  arma_vec x_vec = x.col(0);
  precision_t dx = x_vec(1) - x_vec(0);
  precision_t area = dx * dx;
  arma_mat jacobian = grid.sqrt_g_scgc.slice(iAlt);

  /** State/Velocity extraction **/
  /* MASS DENSITY */
  arma_mat rho = states[0];

  /* VELOCITY */
  // Convert to contravariant (reference) velocity
  arma_mat xVel = states[1]; // u^1
  arma_mat yVel = states[2]; // u^2

  // Generate velocity magnitude squared
  arma_mat vel2 = xVel % xVel + yVel % yVel;

  /** Advancing **/
  /* Initialize projection constructs and storages */
  projection_struct rhoP;
  projection_struct xVelP;
  projection_struct yVelP;

  // They are all pure scalar fields without sqrt(g)
  arma_mat rhoL, rhoR, rhoD, rhoU;
  arma_mat xVelL, xVelR, xVelD, xVelU;
  arma_mat yVelL, yVelR, yVelD, yVelU;

  arma_mat velL2, velR2, velD2, velU2;

  /** Initialize Flux and Wave Speed Storages */
  arma_mat eq1FluxLR, eq1FluxDU;
  arma_mat eq1FluxL, eq1FluxR, eq1FluxD, eq1FluxU;

  arma_mat wsL, wsR, wsD, wsU, wsLR, wsDU;

  arma_mat diff; // for Riemann Solver

  /* Projection */
  rhoP = project_to_edges(rho, x, xEdges, y, yEdges, nGCs);
  xVelP = project_to_edges(xVel, x, xEdges, y, yEdges, nGCs);
  yVelP = project_to_edges(yVel, x, xEdges, y, yEdges, nGCs);

  // Resolve Scalar Fields into rho, xVel, yVel, and totalE (without rho)
  rhoL = rhoP.L;
  rhoR = rhoP.R;
  rhoD = rhoP.D;
  rhoU = rhoP.U;

  xVelL = xVelP.L;
  xVelR = xVelP.R;
  xVelD = xVelP.D;
  xVelU = xVelP.U;

  yVelL = yVelP.L;
  yVelR = yVelP.R;
  yVelD = yVelP.D;
  yVelU = yVelP.U;

  velL2 = xVelL % xVelL + yVelL % yVelL;
  velR2 = xVelR % xVelR + yVelR % yVelR;
  velD2 = xVelD % xVelD + yVelD % yVelD;
  velU2 = xVelU % xVelU + yVelU % yVelU;

  /* Calculate Edge Fluxes */
  // Note that dot product between normal vector at edge and flux vector
  // resolves into a pure one component flux or either hat{x} or hat{y}

  // Flux calculated from the left of the edge
  eq1FluxL = rhoL % xVelL % grid.sqrt_g_Left.slice(iAlt);
  // Flux calculated from the right of the edge
  eq1FluxR = rhoR % xVelR % grid.sqrt_g_Left.slice(iAlt);
  // Flux calculated from the down of the edge
  eq1FluxD = rhoD % yVelD % grid.sqrt_g_Down.slice(iAlt);
  // Flux calculated from the up of the edge
  eq1FluxU = rhoU % yVelU % grid.sqrt_g_Down.slice(iAlt);

  /* Wave Speed Calculation */
  wsL = sqrt(velL2);
  wsR = sqrt(velR2);
  wsD = sqrt(velD2);
  wsU = sqrt(velU2);

  wsLR = wsR;

  for (int i = 0; i < nXs + 1; i++) {
    for (int j = 0; j < nYs; j++) {
      if (wsL(i, j) > wsLR(i, j))
        wsLR(i, j) = wsL(i, j);
    }
  }

  wsDU = wsD;

  for (int i = 0; i < nXs; i++) {
    for (int j = 0; j < nYs + 1; j++) {
      if (wsU(i, j) > wsDU(i, j))
        wsDU(i, j) = wsU(i, j);
    }
  }

  /* Calculate average flux at the edges (Rusanov Flux) */
  diff = (rhoR - rhoL) % grid.sqrt_g_Left.slice(
           iAlt); // State difference, need to add sqrt(g)
  eq1FluxLR = (eq1FluxL + eq1FluxR) / 2 + 0.5 * wsLR % diff;
  diff = (rhoU - rhoD) % grid.sqrt_g_Down.slice(iAlt);
  eq1FluxDU = (eq1FluxD + eq1FluxU) / 2 + 0.5 * wsDU % diff;

  // Setup residual storage for return
  arma_mat eq1_residual(nXs, nYs, fill::zeros);

  // State Update
  // Note the ghost cells WILL NOT BE UPDATED
  for (int j = nGCs; j < nYs - nGCs; j++) {
    for (int i = nGCs; i < nXs - nGCs; i++) {
      precision_t rhoResidual_ij = dx * eq1FluxLR(i + 1, j) -
                                   dx * eq1FluxLR(i, j) +
                                   dx * eq1FluxDU(i, j + 1) -
                                   dx * eq1FluxDU(i, j);
      eq1_residual(i, j) = -1 / area * rhoResidual_ij;
    }
  }

  // Setup return vector
  std::vector<arma_mat> return_vector;
  return_vector.push_back(eq1_residual);

  return return_vector;
}

void Neutrals::solver_horizontal_rusanov_advection(Grid& grid, Times& time) {
  // Function Reporting
  std::string function = "Neutrals::solver_horizontal_rusanov_advection";
  static int iFunction = -1;
  report.enter(function, iFunction);

  // Dimensions of Spatial Discretization
  int64_t nXs = grid.get_nX();
  int64_t nYs = grid.get_nY();
  int64_t nGCs = grid.get_nGCs();
  int64_t nAlts = grid.get_nAlts();
  int iAlt, iSpec;

  // Time Discretization (TODO: change dt calculation method)
  precision_t dt = time.get_dt();


  // Advance for bulk calculation first, calculate for every altitude
  for (iAlt = nGCs; iAlt < nAlts - nGCs; iAlt++) {
    /** Extract Grid Features **/
    arma_mat x = grid.refx_scgc.slice(iAlt);
    arma_mat xEdges = grid.refx_Left.slice(iAlt);
    arma_mat y = grid.refy_scgc.slice(iAlt);
    arma_mat yEdges = grid.refy_Down.slice(iAlt);

    // Get reference grid dimensions (Assume dx = dy and equidistant)
    arma_vec x_vec = x.col(0);
    precision_t dx = x_vec(1) - x_vec(0);
    precision_t area = dx * dx;
    arma_mat jacobian = grid.sqrt_g_scgc.slice(iAlt);

    /** States preprocessing **/
    /* MASS DENSITY */
    arma_mat rho = rho_scgc.slice(iAlt);

    /* VELOCITY */
    // Get spherical velocity
    arma_mat uVel = velocity_vcgc[0].slice(iAlt);
    arma_mat vVel = velocity_vcgc[1].slice(iAlt);
    // Convert to contravariant (reference) velocity
    arma_mat xVel = uVel % grid.A11_inv_scgc.slice(iAlt) + vVel %
                    grid.A12_inv_scgc.slice(iAlt); // u^1
    arma_mat yVel = uVel % grid.A21_inv_scgc.slice(iAlt) + vVel %
                    grid.A22_inv_scgc.slice(iAlt); // u^2

    // Generate velocity magnitude squared
    arma_mat vel2 = xVel % xVel + yVel % yVel;

    /** Advancing **/
    /* Initialize projection constructs and storages */
    projection_struct rhoP;
    projection_struct xVelP;
    projection_struct yVelP;

    // They are all pure scalar fields without sqrt(g)
    arma_mat rhoL, rhoR, rhoD, rhoU;
    arma_mat xVelL, xVelR, xVelD, xVelU;
    arma_mat yVelL, yVelR, yVelD, yVelU;

    arma_mat velL2, velR2, velD2, velU2;

    /** Initialize Flux and Wave Speed Storages */
    arma_mat eq1FluxLR, eq1FluxDU;
    arma_mat eq1FluxL, eq1FluxR, eq1FluxD, eq1FluxU;

    arma_mat wsL, wsR, wsD, wsU, wsLR, wsDU;

    arma_mat diff; // for Riemann Solver

    /* Projection */
    rhoP = project_to_edges(rho, x, xEdges, y, yEdges, nGCs);
    xVelP = project_to_edges(xVel, x, xEdges, y, yEdges, nGCs);
    yVelP = project_to_edges(yVel, x, xEdges, y, yEdges, nGCs);

    // Resolve Scalar Fields into rho, xVel, yVel, and totalE (without rho)
    rhoL = rhoP.L;
    rhoR = rhoP.R;
    rhoD = rhoP.D;
    rhoU = rhoP.U;

    xVelL = xVelP.L;
    xVelR = xVelP.R;
    xVelD = xVelP.D;
    xVelU = xVelP.U;

    yVelL = yVelP.L;
    yVelR = yVelP.R;
    yVelD = yVelP.D;
    yVelU = yVelP.U;

    velL2 = xVelL % xVelL + yVelL % yVelL;
    velR2 = xVelR % xVelR + yVelR % yVelR;
    velD2 = xVelD % xVelD + yVelD % yVelD;
    velU2 = xVelU % xVelU + yVelU % yVelU;


    /* Calculate Edge Fluxes */
    // Note that dot product between normal vector at edge and flux vector
    // resolves into a pure one component flux or either hat{x} or hat{y}

    // Flux calculated from the left of the edge
    eq1FluxL = rhoL % xVelL % grid.sqrt_g_Left.slice(iAlt);
    // Flux calculated from the right of the edge
    eq1FluxR = rhoR % xVelR % grid.sqrt_g_Left.slice(iAlt);
    // Flux calculated from the down of the edge
    eq1FluxD = rhoD % yVelD % grid.sqrt_g_Down.slice(iAlt);
    // Flux calculated from the up of the edge
    eq1FluxU = rhoU % yVelU % grid.sqrt_g_Down.slice(iAlt);

    /* Wave Speed Calculation */
    wsL = sqrt(velL2);
    wsR = sqrt(velR2);
    wsD = sqrt(velD2);
    wsU = sqrt(velU2);

    wsLR = wsR;

    for (int i = 0; i < nXs + 1; i++) {
      for (int j = 0; j < nYs; j++) {
        if (wsL(i, j) > wsLR(i, j))
          wsLR(i, j) = wsL(i, j);
      }
    }

    wsDU = wsD;

    for (int i = 0; i < nXs; i++) {
      for (int j = 0; j < nYs + 1; j++) {
        if (wsU(i, j) > wsDU(i, j))
          wsDU(i, j) = wsU(i, j);
      }
    }

    /* Calculate average flux at the edges (Rusanov Flux) */
    diff = (rhoR - rhoL) % grid.sqrt_g_Left.slice(
             iAlt); // State difference, need to add sqrt(g)
    eq1FluxLR = (eq1FluxL + eq1FluxR) / 2 + 0.5 * wsLR % diff;
    diff = (rhoU - rhoD) % grid.sqrt_g_Down.slice(iAlt);
    eq1FluxDU = (eq1FluxD + eq1FluxU) / 2 + 0.5 * wsDU % diff;

    /* Update Bulk Scalars and Contravariant velocity */
    // Euler State Update
    for (int j = nGCs; j < nYs - nGCs; j++) {
      for (int i = nGCs; i < nXs - nGCs; i++) {
        precision_t rhoResidual_ij = dx * eq1FluxLR(i + 1, j) -
                                     dx * eq1FluxLR(i, j) +
                                     dx * eq1FluxDU(i, j + 1) -
                                     dx * eq1FluxDU(i, j);
        rho(i, j) = rho(i, j) - dt / area / jacobian(i, j) * rhoResidual_ij;
      }
    }

    /* Re-derive Spherical Velocity and Bulk States */
    // Density
    rho_scgc.slice(iAlt) = rho;

    // Bulk Velocity
    //vel2 = xVel % xVel + yVel % yVel; // Squared Magnitude of Contravariant
    //velocity_vcgc[0].slice(iAlt) = xVel%grid.A11_scgc.slice(iAlt) + yVel%grid.A12_scgc.slice(iAlt);
    //velocity_vcgc[1].slice(iAlt) = xVel%grid.A21_scgc.slice(iAlt) + yVel%grid.A22_scgc.slice(iAlt);

    /* Update specie density */
    for (iSpec = 0; iSpec < nSpecies; iSpec++) {
      //species[iSpec].density_scgc.slice(iAlt) = rho % species[iSpec].concentration_scgc.slice(iAlt);
      species[iSpec].density_scgc.slice(iAlt) = rho / species[iSpec].mass;
    }


    report.exit(function);
    return;
  }
}

void Neutrals::solver_horizontal_RK1_advection(Grid& grid, Times& time) {
  // Function Reporting
  std::string function = "Neutrals::solver_horizontal_RK1_advection";
  static int iFunction = -1;
  report.enter(function, iFunction);

  // Dimensions of Spatial Discretization
  int64_t nXs = grid.get_nX();
  int64_t nYs = grid.get_nY();
  int64_t nGCs = grid.get_nGCs();
  int64_t nAlts = grid.get_nAlts();
  int iAlt, iSpec;

  // Time Discretization (TODO: change dt calculation method)
  precision_t dt = time.get_dt();

  // Advance for bulk calculation first, calculate for every altitude
  for (iAlt = nGCs; iAlt < nAlts - nGCs; iAlt++) {
    /** Extract Grid Features **/
    arma_mat x = grid.refx_scgc.slice(iAlt);
    arma_mat xEdges = grid.refx_Left.slice(iAlt);
    arma_mat y = grid.refy_scgc.slice(iAlt);
    arma_mat yEdges = grid.refy_Down.slice(iAlt);

    // Get reference grid dimensions (Assume dx = dy and equidistant)
    arma_vec x_vec = x.col(0);
    precision_t dx = x_vec(1) - x_vec(0);
    precision_t area = dx * dx;
    arma_mat jacobian = grid.sqrt_g_scgc.slice(iAlt);

    /** States preprocessing **/
    /* MASS DENSITY */
    arma_mat rho = rho_scgc.slice(iAlt);

    /* VELOCITY */
    // Get spherical velocity
    arma_mat uVel = velocity_vcgc[0].slice(iAlt);
    arma_mat vVel = velocity_vcgc[1].slice(iAlt);
    // Convert to contravariant (reference) velocity
    arma_mat xVel = uVel % grid.A11_inv_scgc.slice(iAlt) + vVel %
                    grid.A12_inv_scgc.slice(iAlt); // u^1
    arma_mat yVel = uVel % grid.A21_inv_scgc.slice(iAlt) + vVel %
                    grid.A22_inv_scgc.slice(iAlt); // u^2

    /** Advancing with RK4 **/
    // Setup Containers
    arma_mat rho_0 = rho;

    // FIRST (1) STEP, Compute F_0-> State_1
    // Pass in state vector
    std::vector<arma_mat> state_0;
    state_0.push_back(rho_0);
    state_0.push_back(xVel);
    state_0.push_back(yVel);
    std::vector<arma_mat> f_0_vec = residual_horizontal_rusanov_advection(state_0,
                                    grid, time);

    // Extract Gradients
    arma_mat f_0_eq1 = f_0_vec[0];

    /* Update Bulk Scalars and Contravariant velocity */
    // Euler State Update
    for (int j = nGCs; j < nYs - nGCs; j++) {
      for (int i = nGCs; i < nXs - nGCs; i++)
        rho(i, j) = rho_0(i, j) + dt * f_0_eq1(i, j) / jacobian(i, j);
    }

    /* Re-derive Spherical Velocity and Bulk States */
    // Density
    rho_scgc.slice(iAlt) = rho;

    // Bulk Velocity
    //vel2 = xVel % xVel + yVel % yVel; // Squared Magnitude of Contravariant
    //velocity_vcgc[0].slice(iAlt) = xVel%grid.A11_scgc.slice(iAlt) + yVel%grid.A12_scgc.slice(iAlt);
    //velocity_vcgc[1].slice(iAlt) = xVel%grid.A21_scgc.slice(iAlt) + yVel%grid.A22_scgc.slice(iAlt);

    /* Update specie density */
    for (iSpec = 0; iSpec < nSpecies; iSpec++) {
      //species[iSpec].density_scgc.slice(iAlt) = rho % species[iSpec].concentration_scgc.slice(iAlt);
      species[iSpec].density_scgc.slice(iAlt) = rho / species[iSpec].mass;
    }


    report.exit(function);
    return;
  }
}

void Neutrals::solver_horizontal_RK2_advection(Grid& grid, Times& time) {
  // Function Reporting
  std::string function = "Neutrals::solver_horizontal_RK2_advection";
  static int iFunction = -1;
  report.enter(function, iFunction);

  // Dimensions of Spatial Discretization
  int64_t nXs = grid.get_nX();
  int64_t nYs = grid.get_nY();
  int64_t nGCs = grid.get_nGCs();
  int64_t nAlts = grid.get_nAlts();
  int iAlt, iSpec;

  // Time Discretization (TODO: change dt calculation method)
  precision_t dt = time.get_dt();

  // Advance for bulk calculation first, calculate for every altitude
  for (iAlt = nGCs; iAlt < nAlts - nGCs; iAlt++) {
    /** Extract Grid Features **/
    arma_mat x = grid.refx_scgc.slice(iAlt);
    arma_mat xEdges = grid.refx_Left.slice(iAlt);
    arma_mat y = grid.refy_scgc.slice(iAlt);
    arma_mat yEdges = grid.refy_Down.slice(iAlt);

    // Get reference grid dimensions (Assume dx = dy and equidistant)
    arma_vec x_vec = x.col(0);
    precision_t dx = x_vec(1) - x_vec(0);
    precision_t area = dx * dx;
    arma_mat jacobian = grid.sqrt_g_scgc.slice(iAlt);

    /** States preprocessing **/
    /* MASS DENSITY */
    arma_mat rho = rho_scgc.slice(iAlt);

    /* VELOCITY */
    // Get spherical velocity
    arma_mat uVel = velocity_vcgc[0].slice(iAlt);
    arma_mat vVel = velocity_vcgc[1].slice(iAlt);
    // Convert to contravariant (reference) velocity
    arma_mat xVel = uVel % grid.A11_inv_scgc.slice(iAlt) + vVel %
                    grid.A12_inv_scgc.slice(iAlt); // u^1
    arma_mat yVel = uVel % grid.A21_inv_scgc.slice(iAlt) + vVel %
                    grid.A22_inv_scgc.slice(iAlt); // u^2

    /** Advancing with RK4 **/
    // Setup Containers
    arma_mat rho_0 = rho;
    arma_mat rho_1(nXs, nYs, fill::zeros); // corresponding f_1

    // FIRST (1) STEP, Compute F_0-> State_1
    // Pass in state vector
    std::vector<arma_mat> state_0;
    state_0.push_back(rho_0);
    state_0.push_back(xVel);
    state_0.push_back(yVel);
    std::vector<arma_mat> f_0_vec = residual_horizontal_hlle_advection(state_0,
                                    grid, time);
    // Extract Gradients
    arma_mat f_0_eq1 = f_0_vec[0];

    /* Update Bulk Scalars and Contravariant velocity */
    // Euler State Update
    for (int j = nGCs; j < nYs - nGCs; j++) {
      for (int i = nGCs; i < nXs - nGCs; i++)
        rho_1(i, j) = rho_0(i, j) + dt * f_0_eq1(i, j) / jacobian(i, j);
    }

    // SECOND (2) STEP, Compute F_1-> State_2
    // Pass in state vector
    std::vector<arma_mat> state_1;
    state_1.push_back(rho_1);
    state_1.push_back(xVel);
    state_1.push_back(yVel);
    std::vector<arma_mat> f_1_vec = residual_horizontal_hlle_advection(state_1,
                                    grid, time);
    // Extract Gradients
    arma_mat f_1_eq1 = f_1_vec[0];

    // Summing all steps for final update
    arma_mat f_sum_eq1 = f_0_eq1 + f_1_eq1;

    /* Update Bulk Scalars and Contravariant velocity */
    // Euler State Update
    for (int j = nGCs; j < nYs - nGCs; j++) {
      for (int i = nGCs; i < nXs - nGCs; i++)
        rho(i, j) = rho(i, j) + 0.5 * dt * f_sum_eq1(i, j) / jacobian(i, j);
    }

    /* Re-derive Spherical Velocity and Bulk States */
    // Density
    rho_scgc.slice(iAlt) = rho;

    // Bulk Velocity
    //vel2 = xVel % xVel + yVel % yVel; // Squared Magnitude of Contravariant
    //velocity_vcgc[0].slice(iAlt) = xVel%grid.A11_scgc.slice(iAlt) + yVel%grid.A12_scgc.slice(iAlt);
    //velocity_vcgc[1].slice(iAlt) = xVel%grid.A21_scgc.slice(iAlt) + yVel%grid.A22_scgc.slice(iAlt);

    /* Update specie density */
    for (iSpec = 0; iSpec < nSpecies; iSpec++) {
      //species[iSpec].density_scgc.slice(iAlt) = rho % species[iSpec].concentration_scgc.slice(iAlt);
      species[iSpec].density_scgc.slice(iAlt) = rho / species[iSpec].mass;
    }


    report.exit(function);
    return;
  }
}

void Neutrals::solver_horizontal_RK4_advection(Grid& grid, Times& time) {
  // Function Reporting
  std::string function = "Neutrals::solver_horizontal_RK4_advection";
  static int iFunction = -1;
  report.enter(function, iFunction);

  // Dimensions of Spatial Discretization
  int64_t nXs = grid.get_nX();
  int64_t nYs = grid.get_nY();
  int64_t nGCs = grid.get_nGCs();
  int64_t nAlts = grid.get_nAlts();
  int iAlt, iSpec;

  // Time Discretization (TODO: change dt calculation method)
  precision_t dt = time.get_dt();

  // Advance for bulk calculation first, calculate for every altitude
  for (iAlt = nGCs; iAlt < nAlts - nGCs; iAlt++) {
    /** Extract Grid Features **/
    arma_mat x = grid.refx_scgc.slice(iAlt);
    arma_mat xEdges = grid.refx_Left.slice(iAlt);
    arma_mat y = grid.refy_scgc.slice(iAlt);
    arma_mat yEdges = grid.refy_Down.slice(iAlt);

    // Get reference grid dimensions (Assume dx = dy and equidistant)
    arma_vec x_vec = x.col(0);
    precision_t dx = x_vec(1) - x_vec(0);
    precision_t area = dx * dx;
    arma_mat jacobian = grid.sqrt_g_scgc.slice(iAlt);

    /** States preprocessing **/
    /* MASS DENSITY */
    arma_mat rho = rho_scgc.slice(iAlt);

    /* VELOCITY */
    // Get spherical velocity
    arma_mat uVel = velocity_vcgc[0].slice(iAlt);
    arma_mat vVel = velocity_vcgc[1].slice(iAlt);
    // Convert to contravariant (reference) velocity
    arma_mat xVel = uVel % grid.A11_inv_scgc.slice(iAlt) + vVel %
                    grid.A12_inv_scgc.slice(iAlt); // u^1
    arma_mat yVel = uVel % grid.A21_inv_scgc.slice(iAlt) + vVel %
                    grid.A22_inv_scgc.slice(iAlt); // u^2

    /** Advancing with RK4 **/
    // Setup Containers
    arma_mat rho_0 = rho;
    arma_mat rho_1(nXs, nYs, fill::zeros); // corresponding f_1
    arma_mat rho_2(nXs, nYs, fill::zeros); // corresponding f_2
    arma_mat rho_3(nXs, nYs, fill::zeros); // corresponding f_3

    // FIRST (1) STEP, Compute F_0-> State_1
    // Pass in state vector
    std::vector<arma_mat> state_0;
    state_0.push_back(rho_0);
    state_0.push_back(xVel);
    state_0.push_back(yVel);
    std::vector<arma_mat> f_0_vec = residual_horizontal_rusanov_advection(state_0,
                                    grid, time);
    // Extract Gradients
    arma_mat f_0_eq1 = f_0_vec[0];

    /* Update Bulk Scalars and Contravariant velocity */
    // Euler State Update
    for (int j = nGCs; j < nYs - nGCs; j++) {
      for (int i = nGCs; i < nXs - nGCs; i++)
        rho_1(i, j) = rho_0(i, j) + 0.5 * dt * f_0_eq1(i, j) / jacobian(i, j);
    }

    // SECOND (2) STEP, Compute F_1-> State_2
    // Pass in state vector
    std::vector<arma_mat> state_1;
    state_1.push_back(rho_1);
    state_1.push_back(xVel);
    state_1.push_back(yVel);
    std::vector<arma_mat> f_1_vec = residual_horizontal_rusanov_advection(state_1,
                                    grid, time);
    // Extract Gradients
    arma_mat f_1_eq1 = f_1_vec[0];

    /* Update Bulk Scalars and Contravariant velocity */
    // Euler State Update
    for (int j = nGCs; j < nYs - nGCs; j++) {
      for (int i = nGCs; i < nXs - nGCs; i++)
        rho_2(i, j) = rho_0(i, j) + 0.5 * dt * f_1_eq1(i, j) / jacobian(i, j);
    }

    // THIRD (3) STEP, Compute F_2-> State_3
    // Pass in state vector
    std::vector<arma_mat> state_2;
    state_2.push_back(rho_2);
    state_2.push_back(xVel);
    state_2.push_back(yVel);
    std::vector<arma_mat> f_2_vec = residual_horizontal_rusanov_advection(state_2,
                                    grid, time);
    // Extract Gradients
    arma_mat f_2_eq1 = f_2_vec[0];

    /* Update Bulk Scalars and Contravariant velocity */
    // Euler State Update
    for (int j = nGCs; j < nYs - nGCs; j++) {
      for (int i = nGCs; i < nXs - nGCs; i++)
        rho_3(i, j) = rho_0(i, j) + dt * f_2_eq1(i, j) / jacobian(i, j);
    }

    // FOURTH (4) STEP, Compute F_3
    // Pass in state vector
    std::vector<arma_mat> state_3;
    state_3.push_back(rho_3);
    state_3.push_back(xVel);
    state_3.push_back(yVel);
    std::vector<arma_mat> f_3_vec = residual_horizontal_rusanov_advection(state_3,
                                    grid, time);
    // Extract Gradients
    arma_mat f_3_eq1 = f_3_vec[0];

    // Summing all steps for final update
    arma_mat f_sum_eq1 = f_0_eq1 + 2 * f_1_eq1 + 2 * f_2_eq1 + f_3_eq1;

    /* Update Bulk Scalars and Contravariant velocity */
    // Euler State Update
    for (int j = nGCs; j < nYs - nGCs; j++) {
      for (int i = nGCs; i < nXs - nGCs; i++)
        rho(i, j) = rho(i, j) + dt / 6 * f_sum_eq1(i, j) / jacobian(i, j);
    }

    /* Re-derive Spherical Velocity and Bulk States */
    // Density
    rho_scgc.slice(iAlt) = rho;

    // Bulk Velocity
    //vel2 = xVel % xVel + yVel % yVel; // Squared Magnitude of Contravariant
    //velocity_vcgc[0].slice(iAlt) = xVel%grid.A11_scgc.slice(iAlt) + yVel%grid.A12_scgc.slice(iAlt);
    //velocity_vcgc[1].slice(iAlt) = xVel%grid.A21_scgc.slice(iAlt) + yVel%grid.A22_scgc.slice(iAlt);

    /* Update specie density */
    for (iSpec = 0; iSpec < nSpecies; iSpec++) {
      //species[iSpec].density_scgc.slice(iAlt) = rho % species[iSpec].concentration_scgc.slice(iAlt);
      species[iSpec].density_scgc.slice(iAlt) = rho / species[iSpec].mass;
    }


    report.exit(function);
    return;
  }
}