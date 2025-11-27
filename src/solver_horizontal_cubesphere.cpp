// Copyright 2023, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

// Initial version: F. Cheng, July 2023
// Moved to new solver: August 2025

#include "aether.h"

// ---------------------------------------------------------
// Update States
// ---------------------------------------------------------

void update_states_cubesphere(arma_mat rho,
                              arma_mat &xVel,
                              arma_mat &yVel,
                              arma_mat &temp,
                              arma_mat &drhodt,
                              arma_mat &dlonVeldt,
                              arma_mat &dlatVeldt,
                              arma_mat &dtempdt,
                              cubesphere_chars gridC,
                              cubesphere_chars gridL,
                              cubesphere_chars gridD,
                              precision_t dt,
                              int64_t iZ) {

  arma_mat xMomentum, yMomentum;
  arma_mat rhoE, energy, vel2;

  precision_t cv = 1500.0;

  if (report.test_verbose(2))
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

  if (report.test_verbose(3))
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

  if (report.test_verbose(3))
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

  if (report.test_verbose(3))
    std::cout << "  ---> Normal Velocities\n";

  // Calculate the normal velocity at the boundaries:
  velNormL = xVelP.L % gridL.nXiLon + yVelP.L % gridL.nXiLat;
  velNormR = xVelP.R % gridL.nXiLon + yVelP.R % gridL.nXiLat;
  velNormU = xVelP.U % gridD.nNuLon + yVelP.U % gridD.nNuLat;
  velNormD = xVelP.D % gridD.nNuLon + yVelP.D % gridD.nNuLat;

  if (report.test_verbose(3))
    std::cout << "  ---> Fluxes eq 1\n";

  // Flux calculated from the left of the edge
  eq1FluxL = rhoP.L % velNormL;
  // Flux calculated from the right of the edge
  eq1FluxR = rhoP.R % velNormR;
  // Flux calculated from the down of the edge
  eq1FluxD = rhoP.D % velNormD;
  // Flux calculated from the up of the edge
  eq1FluxU = rhoP.U % velNormU;

  if (report.test_verbose(3))
    std::cout << "  ---> Fluxes eq 2\n";

  eq2FluxL = (xMomentumP.L % velNormL);
  eq2FluxR = (xMomentumP.R % velNormR);
  eq2FluxD = (xMomentumP.D % velNormD);
  eq2FluxU = (xMomentumP.U % velNormU);

  if (report.test_verbose(3))
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
  if (report.test_verbose(3))
    std::cout << "  ---> Diffusive Fluxes\n";

  precision_t cGamma = 5.0 / 3.0;

  wsL.resize(gridC.nXt + 1, gridC.nYt);
  wsR.resize(gridC.nXt + 1, gridC.nYt);
  wsD.resize(gridC.nXt, gridC.nYt + 1);
  wsU.resize(gridC.nXt, gridC.nYt + 1);

  wsL.zeros();
  wsR.zeros();
  wsU.zeros();
  wsD.zeros();

  for (int64_t i = 0; i < gridC.nXt; i++) {
    for (int64_t j = 0; j < gridC.nYt; j++) {
      wsL(i, j) = sqrt(velL2(i, j)) + sqrt(cGamma * (cGamma - 1) * tempP.L(i, j));
      wsR(i, j) = sqrt(velR2(i, j)) + sqrt(cGamma * (cGamma - 1) * tempP.R(i, j));
      wsD(i, j) = sqrt(velD2(i, j)) + sqrt(cGamma * (cGamma - 1) * tempP.D(i, j));
      wsU(i, j) = sqrt(velU2(i, j)) + sqrt(cGamma * (cGamma - 1) * tempP.U(i, j));
    }
  }

  wsLR = wsR;

  for (int64_t i = 0; i < gridC.nXt; i++) {
    for (int64_t j = 0; j < gridC.nYt; j++) {
      if (wsL(i, j) > wsLR(i, j))
        wsLR(i, j) = wsL(i, j);
    }
  }

  wsDU = wsD;

  for (int64_t i = 0; i < gridC.nXt; i++) {
    for (int64_t j = 0; j < gridC.nYt; j++) {
      if (wsU(i, j) > wsDU(i, j))
        wsDU(i, j) = wsU(i, j);
    }
  }

  // ------------------------------------------------
  // Calculate average flux at the edges (Rusanov Flux):

  if (report.test_verbose(3))
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
  if (report.test_verbose(3))
    std::cout << "  ---> Updating equations of state\n";

  precision_t dpdx, dpdn, pp, pm;

  arma_mat ax(gridC.nXt, gridC.nYt), an(gridC.nXt, gridC.nYt);

  ax.zeros();
  an.zeros();
  arma_mat dedt(gridC.nXt, gridC.nYt);
  dedt.zeros();

  arma_mat rhoNew = rho;

  // Only deal with inner cell
  for (int64_t j = gridC.iYfirst_; j < gridC.iYlast_; j++) {
    for (int64_t i = gridC.iXfirst_; i < gridC.iXlast_; i++) {
      precision_t rhoResidual_ij = (gridL.dln(i + 1, j, iZ) * eq1FluxLR(i + 1, j) -
                                    gridL.dln(i, j, iZ) * eq1FluxLR(i, j) +
                                    gridD.dlx(i, j + 1, iZ) * eq1FluxDU(i, j + 1) -
                                    gridD.dlx(i, j, iZ) * eq1FluxDU(i, j));
      drhodt(i, j) = rhoResidual_ij / gridC.dS(i, j, iZ);

      rhoNew(i, j) = rho(i, j) + dt * drhodt(i, j);

      precision_t xMomentumResidual_ij = (gridL.dln(i + 1, j, iZ) * eq2FluxLR(i + 1,
                                          j) -
                                          gridL.dln(i, j, iZ) * eq2FluxLR(i, j) +
                                          gridD.dlx(i, j + 1, iZ) * eq2FluxDU(i, j + 1) -
                                          gridD.dlx(i, j, iZ) * eq2FluxDU(i, j));
      dxVeldt(i, j) = xMomentumResidual_ij / gridC.dS(i, j, iZ) / rhoNew(i, j);

      precision_t yMomentumResidual_ij = (gridL.dln(i + 1, j, iZ) * eq3FluxLR(i + 1,
                                          j) -
                                          gridL.dln(i, j, iZ) * eq3FluxLR(i, j) +
                                          gridD.dlx(i, j + 1, iZ) * eq3FluxDU(i, j + 1) -
                                          gridD.dlx(i, j, iZ) * eq3FluxDU(i, j));
      dyVeldt(i, j) = yMomentumResidual_ij / gridC.dS(i, j, iZ) / rhoNew(i, j);

      // Calculate the gradient in the potential in the cubesphere
      // coordinate system:
      dpdx = 1 / gridC.R(iZ) * gridC.D(i, j) *
             (pressureLR(i + 1, j) - pressureLR(i, j)) / gridC.dxi;
      dpdn = 1 / gridC.R(iZ) * gridC.X(i, j) * gridC.Y(i, j) /
             gridC.D(i, j) *
             (pressureDU(i, j + 1) - pressureDU(i, j)) / gridC.dnu;
      ax(i, j) = (dpdx + dpdn) / rhoNew(i, j);

      dpdx = 1 / gridC.R(iZ) * gridC.X(i, j) * gridC.Y(i, j) /
             gridC.C(i, j) * (pressureLR(i + 1, j) - pressureLR(i, j)) / gridC.dxi;
      dpdn = 1 / gridC.R(iZ) * gridC.C(i, j) *
             (pressureDU(i, j + 1) - pressureDU(i, j)) / gridC.dnu;
      an(i, j) = (dpdx + dpdn) / rhoNew(i, j);

      precision_t energyResidual_ij = (gridL.dln(i + 1, j, iZ) * eq4FluxLR(i + 1, j) -
                                       gridL.dln(i, j, iZ) * eq4FluxLR(i, j) +
                                       gridD.dlx(i, j + 1, iZ) * eq4FluxDU(i, j + 1) -
                                       gridD.dlx(i, j, iZ) * eq4FluxDU(i, j));
      dedt(i, j) = energyResidual_ij / gridC.dS(i, j, iZ);

    }
  }

  // lat is negative because of the Rochi definition of theta:
  dlatVeldt = dyVeldt - (ax % gridC.Atx + an % gridC.Atn);
  dlonVeldt = dxVeldt + ax % gridC.Apx + an % gridC.Apn;
  dtempdt = dedt / rhoNew / cv;

  return;
}


// using namespace Cubesphere_tools;

std::vector<arma_mat> Neutrals::residual_horizontal_rusanov(
  std::vector<arma_mat>& states,
  Grid & grid,
  Times & time,
  int64_t iAlt) {

  // Dimensions of Spatial Discretization
  int64_t nXs = grid.get_nX();
  int64_t nYs = grid.get_nY();
  int64_t nGCs = grid.get_nGCs();
  int64_t nAlts = grid.get_nAlts();

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
  arma_mat rho = states[0];

  /* VELOCITY */
  // Get contravariant velocity
  //arma_mat xVel = states[1]; // u^1
  //arma_mat yVel = states[2]; // u^2

  // Generate contravriant momentum
  arma_mat xMomentum = states[1]; // x1momentum
  arma_mat yMomentum = states[2]; // x2momentum

  // Resolve to contravariant velocity
  arma_mat xVel = xMomentum / rho; // u^1
  arma_mat yVel = yMomentum / rho; // u^2

  // Generate velocity magnitude squared
  arma_mat vel2 = xVel % xVel + yVel % yVel;

  /* TEMP and ENERGY */
  // Generate total energy (rhoE) (TODO: Verify)
  arma_mat rhoE = states[3];

  /** Advancing **/
  /* Initialize projection constructs and storages */
  projection_struct rhoP;
  projection_struct xMomentumP;
  projection_struct yMomentumP;
  projection_struct rhoEP;
  projection_struct gammaP;
  projection_struct tempP;
  projection_struct numberDensityP;

  // They are all pure scalar fields without sqrt(g)
  arma_mat rhoL, rhoR, rhoD, rhoU;
  arma_mat xVelL, xVelR, xVelD, xVelU;
  arma_mat yVelL, yVelR, yVelD, yVelU;
  arma_mat totalEL, totalER, totalED, totalEU;

  arma_mat velL2, velR2, velD2, velU2;
  arma_mat internaleL, internaleR, internaleD, internaleU;
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

  /* Projection */
  rhoP = project_to_edges(rho, x, xEdges, y, yEdges, nGCs);
  xMomentumP = project_to_edges(xMomentum, x, xEdges, y, yEdges, nGCs);
  yMomentumP = project_to_edges(yMomentum, x, xEdges, y, yEdges, nGCs);
  rhoEP = project_to_edges(rhoE, x, xEdges, y, yEdges, nGCs);
  // Also need to project gamma and temp - these should be passed, since
  // they need to be updated for the RK4 scheme:
  gammaP = project_to_edges(gamma_scgc.slice(iAlt), x, xEdges, y, yEdges, nGCs);
  tempP = project_to_edges(temperature_scgc.slice(iAlt), x, xEdges, y, yEdges,
                           nGCs);
  numberDensityP = project_to_edges(density_scgc.slice(iAlt), x, xEdges, y,
                                    yEdges,
                                    nGCs);

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

  totalEL = rhoEP.L / rhoL;
  totalER = rhoEP.R / rhoR;
  totalED = rhoEP.D / rhoD;
  totalEU = rhoEP.U / rhoU;

  velL2 = xVelL % xVelL + yVelL % yVelL;
  velR2 = xVelR % xVelR + yVelR % yVelR;
  velD2 = xVelD % xVelD + yVelD % yVelD;
  velU2 = xVelU % xVelU + yVelU % yVelU;

  internaleL = totalEL - 0.5 * velL2;
  internaleR = totalER - 0.5 * velR2;
  internaleD = totalED - 0.5 * velD2;
  internaleU = totalEU - 0.5 * velU2;

  //pressureL = (gammaP.L - 1) % (rhoP.L % internaleL);
  //pressureR = (gammaP.R - 1) % (rhoP.R % internaleR);
  //pressureD = (gammaP.D - 1) % (rhoP.D % internaleD);
  //pressureU = (gammaP.U - 1) % (rhoP.U % internaleU);

  pressureL = cKB * (numberDensityP.L % tempP.L);
  pressureR = cKB * (numberDensityP.R % tempP.R);
  pressureD = cKB * (numberDensityP.D % tempP.D);
  pressureU = cKB * (numberDensityP.U % tempP.U);

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
  /*
    eq2FluxL = (rhoL % xVelL % xVelL +
                pressureL % grid.g11_upper_Left.slice(iAlt)) %
               grid.sqrt_g_Left.slice(iAlt);
    eq2FluxR = (rhoR % xVelR % xVelR +
                pressureR % grid.g11_upper_Left.slice(iAlt)) %
               grid.sqrt_g_Left.slice(iAlt);
    eq2FluxD = (rhoD % yVelD % xVelD +
                pressureD % grid.g12_upper_Down.slice(iAlt)) %
               grid.sqrt_g_Down.slice(iAlt);
    eq2FluxU = (rhoU % yVelU % xVelU +
                pressureU % grid.g12_upper_Down.slice(iAlt)) %
               grid.sqrt_g_Down.slice(iAlt);
  */
  eq2FluxL = (rhoL % xVelL % xVelL +
              pressureL) %
             grid.sqrt_g_Left.slice(iAlt);
  eq2FluxR = (rhoR % xVelR % xVelR +
              pressureR) %
             grid.sqrt_g_Left.slice(iAlt);
  eq2FluxD = (rhoD % yVelD % xVelD +
              pressureD) %
             grid.sqrt_g_Down.slice(iAlt);
  eq2FluxU = (rhoU % yVelU % xVelU +
              pressureU) %
             grid.sqrt_g_Down.slice(iAlt);
  /*
    eq3FluxL = (rhoL % xVelL % yVelL +
                pressureL % grid.g21_upper_Left.slice(iAlt)) %
               grid.sqrt_g_Left.slice(iAlt);
    eq3FluxR = (rhoR % xVelR % yVelR +
                pressureR % grid.g21_upper_Left.slice(iAlt)) %
               grid.sqrt_g_Left.slice(iAlt);
    eq3FluxD = (rhoD % yVelD % yVelD +
                pressureD % grid.g22_upper_Down.slice(iAlt)) %
               grid.sqrt_g_Down.slice(iAlt);
    eq3FluxU = (rhoU % yVelU % yVelU +
                pressureU % grid.g22_upper_Down.slice(iAlt)) %
               grid.sqrt_g_Down.slice(iAlt);
  */
  eq3FluxL = (rhoL % xVelL % yVelL +
              pressureL) %
             grid.sqrt_g_Left.slice(iAlt);
  eq3FluxR = (rhoR % xVelR % yVelR +
              pressureR) %
             grid.sqrt_g_Left.slice(iAlt);
  eq3FluxD = (rhoD % yVelD % yVelD +
              pressureD) %
             grid.sqrt_g_Down.slice(iAlt);
  eq3FluxU = (rhoU % yVelU % yVelU +
              pressureU) %
             grid.sqrt_g_Down.slice(iAlt);

  eq4FluxL = (rhoEP.L + pressureL) % xVelL % grid.sqrt_g_Left.slice(iAlt);
  eq4FluxR = (rhoEP.R + pressureR) % xVelR % grid.sqrt_g_Left.slice(iAlt);
  eq4FluxD = (rhoEP.D + pressureD) % yVelD % grid.sqrt_g_Down.slice(iAlt);
  eq4FluxU = (rhoEP.U + pressureU) % yVelU % grid.sqrt_g_Down.slice(iAlt);

  /* Wave Speed Calculation */
  wsL = sqrt(velL2) + sqrt(gammaP.L % (gammaP.L - 1.) % tempP.L);
  wsR = sqrt(velR2) + sqrt(gammaP.R % (gammaP.R - 1.) % tempP.R);
  wsD = sqrt(velD2) + sqrt(gammaP.D % (gammaP.D - 1.) % tempP.D);
  wsU = sqrt(velU2) + sqrt(gammaP.U % (gammaP.U - 1.) % tempP.U);

  //wsL = abs(xVelL) + sqrt(gammaP.L % (gammaP.L - 1.) % internaleL);
  //wsR = abs(xVelR) + sqrt(gammaP.R % (gammaP.R - 1.) % internaleR);
  //wsD = abs(yVelD) + sqrt(gammaP.D % (gammaP.D - 1.) % internaleD);
  //wsU = abs(yVelU) + sqrt(gammaP.U % (gammaP.U - 1.) % internaleU);

  // Find the maximum wave speed
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
  /* Why is it + instead of - for the state difference?
   * Because the projection actually works backwards
   * Left states are actually right
   * Right states are actually left
   * Due to the convention in the past codes
   * We keep it this way for consistency
   */

  // State difference, need to add sqrt(g)
  diff = (rhoR - rhoL) % grid.sqrt_g_Left.slice(iAlt);
  eq1FluxLR = (eq1FluxL + eq1FluxR) / 2 + 0.5 * wsLR % diff;
  diff = (rhoU - rhoD) % grid.sqrt_g_Down.slice(iAlt);
  eq1FluxDU = (eq1FluxD + eq1FluxU) / 2 + 0.5 * wsDU % diff;

  if (iAlt == -1) {
    std::cout << "in solver: " << iProc << " " <<
              wsDU(13, 23) << " " <<
              wsU(13, 23) << " " <<
              wsD(13, 23) << " " <<
              gammaP.D(13, 23) << " " <<
              internaleD(13, 23) << " " <<
              gammaP.U(14, 23) << " " <<
              internaleU(13, 23) << " " <<
              diff(13, 22) << "\n";
  }

  diff = (rhoR % xVelR - rhoL % xVelL) % grid.sqrt_g_Left.slice(iAlt);
  eq2FluxLR = (eq2FluxL + eq2FluxR) / 2 + 0.5 * wsLR % diff;
  diff = (rhoU % xVelU - rhoD % xVelD) % grid.sqrt_g_Down.slice(iAlt);
  eq2FluxDU = (eq2FluxD + eq2FluxU) / 2 + 0.5 * wsDU % diff;

  diff = (rhoR % yVelR - rhoL % yVelL) % grid.sqrt_g_Left.slice(iAlt);
  eq3FluxLR = (eq3FluxL + eq3FluxR) / 2 + 0.5 * wsLR % diff;
  diff = (rhoU % yVelU - rhoD % yVelD) % grid.sqrt_g_Down.slice(iAlt);
  eq3FluxDU = (eq3FluxD + eq3FluxU) / 2 + 0.5 * wsDU % diff;

  diff = (rhoR % totalER - rhoL % totalEL) % grid.sqrt_g_Left.slice(iAlt);
  eq4FluxLR = (eq4FluxL + eq4FluxR) / 2 + 0.5 * wsLR % diff;
  diff = (rhoU % totalEU - rhoD % totalED) % grid.sqrt_g_Down.slice(iAlt);
  eq4FluxDU = (eq4FluxD + eq4FluxU) / 2 + 0.5 * wsDU % diff;

  // Setup residual storage for return
  arma_mat eq1_residual(nXs, nYs, fill::zeros);
  arma_mat eq2_residual(nXs, nYs, fill::zeros);
  arma_mat eq3_residual(nXs, nYs, fill::zeros);
  arma_mat eq4_residual(nXs, nYs, fill::zeros);

  // State Update
  // Note the ghost cells WILL NOT BE UPDATED
  for (int j = nGCs; j < nYs - nGCs; j++) {
    for (int i = nGCs; i < nXs - nGCs; i++) {
      precision_t rhoResidual_ij = dx * eq1FluxLR(i + 1, j) -
                                   dx * eq1FluxLR(i, j) +
                                   dx * eq1FluxDU(i, j + 1) -
                                   dx * eq1FluxDU(i, j);
      eq1_residual(i, j) = -1 / area * rhoResidual_ij;
      precision_t xMomentumResidual_ij = dx * eq2FluxLR(i + 1, j) -
                                         dx * eq2FluxLR(i, j) +
                                         dx * eq2FluxDU(i, j + 1) -
                                         dx * eq2FluxDU(i, j);
      eq2_residual(i, j) = -1 / area * xMomentumResidual_ij;
      precision_t yMomentumResidual_ij = dx * eq3FluxLR(i + 1, j) -
                                         dx * eq3FluxLR(i, j) +
                                         dx * eq3FluxDU(i, j + 1) -
                                         dx * eq3FluxDU(i, j);
      eq3_residual(i, j) = -1 / area * yMomentumResidual_ij;
      precision_t rhoEResidual_ij = dx * eq4FluxLR(i + 1, j) -
                                    dx * eq4FluxLR(i, j) +
                                    dx * eq4FluxDU(i, j + 1) -
                                    dx * eq4FluxDU(i, j);
      eq4_residual(i, j) = -1 / area * rhoEResidual_ij;
    }
  }

  if (iAlt == -1) {
    std::cout << "in solver2: " << iProc << " " <<
              eq1_residual(13, 22) << " " <<
              eq2_residual(13, 22) << " " <<
              eq3_residual(13, 22) << " " <<
              eq4_residual(13, 22) << " " <<
              area << "\n";
  }


  // Setup return vector
  std::vector<arma_mat> return_vector;
  return_vector.push_back(eq1_residual);
  return_vector.push_back(eq2_residual);
  return_vector.push_back(eq3_residual);
  return_vector.push_back(eq4_residual);

  return return_vector;
}


//--------------------------------------------------------------------
// New solver using Rochi
//--------------------------------------------------------------------

void Neutrals::solver_horizontal_RK1_rochi(Grid & grid, Times & time) {

  std::string function = "Neutrals::solver_horizontal_RK1_rochi";
  static int iFunction = -1;
  report.enter(function, iFunction);

  precision_t dt = time.get_dt();

  int64_t nXs = grid.get_nX();
  int64_t nYs = grid.get_nY();
  int64_t nGCs = grid.get_nGCs();
  int64_t nAlts = grid.get_nAlts();

  calc_concentration();

  arma_mat temp(nXs, nYs), rho(nXs, nYs), vLon(nXs, nYs), vLat(nXs, nYs);

  arma_mat k1rho(nXs, nYs);
  arma_mat k1vLon(nXs, nYs), k1vLat(nXs, nYs);
  arma_mat k1temp(nXs, nYs);

  int64_t iAlt;

  for (iAlt = nGCs; iAlt < nAlts - nGCs; iAlt++) {

    rho = rho_scgc.slice(iAlt);
    vLon = velocity_vcgc[0].slice(iAlt);
    vLat = velocity_vcgc[1].slice(iAlt);
    temp = temperature_scgc.slice(iAlt);

    // k1 - start at t0, go to t+1/2 to figure out slope at t0 (k1)
    update_states_cubesphere(
      rho, vLon, vLat, temp,
      k1rho, k1vLon, k1vLat, k1temp,
      grid.cubeC, grid.cubeL, grid.cubeD, dt, iAlt);
    // Take full step using k1:
    rho_scgc.slice(iAlt) = rho - k1rho * dt;
    velocity_vcgc[0].slice(iAlt) = vLon - k1vLon * dt;
    velocity_vcgc[1].slice(iAlt) = vLat - k1vLat * dt;
    temperature_scgc.slice(iAlt) = temp - k1temp * dt;
  }

  calc_density_from_mass_concentration();

  report.exit(function);
  return;

}


void Neutrals::solver_horizontal_RK1(Grid & grid, Times & time) {
  // Function Reporting
  std::string function = "Neutrals::solver_horizontal_RK1";
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

  arma_mat x(nXs, nYs), y(nXs, nYs);
  arma_mat jacobian(nXs, nYs), rho(nXs, nYs), rhoE(nXs, nYs), vel2(nXs, nYs);
  arma_mat uVel(nXs, nYs), vVel(nXs, nYs), xVel(nXs, nYs), yVel(nXs, nYs);
  arma_mat xMomentum(nXs, nYs), yMomentum(nXs, nYs);
  arma_mat xMomentum_0(nXs, nYs), yMomentum_0(nXs, nYs);
  arma_mat rho_0(nXs, nYs), rhoE_0(nXs, nYs);
  arma_mat f_0_eq1(nXs, nYs), f_0_eq2(nXs, nYs);
  arma_mat f_0_eq3(nXs, nYs), f_0_eq4(nXs, nYs);

  calc_concentration();

  // Advance for bulk calculation first, calculate for every altitude

  for (iAlt = nGCs; iAlt < nAlts - nGCs; iAlt++) {
    /** Extract Grid Features **/
    x = grid.refx_scgc.slice(iAlt);
    arma_mat xEdges = grid.refx_Left.slice(iAlt);
    y = grid.refy_scgc.slice(iAlt);
    arma_mat yEdges = grid.refy_Down.slice(iAlt);

    // Get reference grid dimensions (Assume dx = dy and equidistant)
    arma_vec x_vec = x.col(0);
    precision_t dx = x_vec(1) - x_vec(0);
    precision_t area = dx * dx;
    jacobian = grid.sqrt_g_scgc.slice(iAlt);

    /** States preprocessing **/
    /* MASS DENSITY */
    rho = rho_scgc.slice(iAlt);

    /* VELOCITY */
    // Get spherical velocity
    uVel = velocity_vcgc[0].slice(iAlt);
    vVel = velocity_vcgc[1].slice(iAlt);
    // Convert to contravariant (reference) velocity
    xVel = uVel % grid.A11_inv_scgc.slice(iAlt) + vVel %
           grid.A12_inv_scgc.slice(iAlt); // u^1
    yVel = uVel % grid.A21_inv_scgc.slice(iAlt) + vVel %
           grid.A22_inv_scgc.slice(iAlt); // u^2
    vel2 = xVel % xVel + yVel % yVel;
    // Generate contravriant momentum (no sqrt(g))
    xMomentum = rho % xVel; // x1momentum
    yMomentum = rho % yVel; // x2momentum

    /* TEMP and ENERGY */
    // Generate total energy (rhoE (no sqrt(g)))
    // (TODO: Verify units)
    rhoE = rho % (temperature_scgc.slice(iAlt) % Cv_scgc.slice(
                    iAlt) + 0.5 * vel2);


    if (iAlt == -1) {
      std::cout << "before solve: " << iProc << " " <<  temperature_scgc(13, 22,
                2) << " " <<
                rhoE(13, 22) << " " << xVel(13, 22) << " " << yVel(13, 22) << " " <<
                rho_scgc(13, 22, 2) << "\n";
    }

    if (iAlt == -2) {
      std::cout << "before solve: " <<
                temperature_scgc(13, 22, 2) << " " <<
                xVel(13, 22) << " " << yVel(13, 22) << " " <<
                rho_scgc(13, 22, 2) << "\n";
    }


    /** Advancing with RK4 **/
    // Setup Containers
    rho_0 = rho;
    xMomentum_0 = xMomentum;
    yMomentum_0 = yMomentum;
    rhoE_0 = rhoE;

    // FIRST (1) STEP, Compute F_0-> State_1
    // Pass in state vector
    std::vector<arma_mat> state_0;
    state_0.push_back(rho_0);
    state_0.push_back(xMomentum_0);
    state_0.push_back(yMomentum_0);
    state_0.push_back(rhoE_0);
    std::vector<arma_mat> f_0_vec = residual_horizontal_rusanov(state_0, grid, time,
                                                                iAlt);
    // Extract Gradients
    f_0_eq1 = f_0_vec[0];
    f_0_eq2 = f_0_vec[1];
    f_0_eq3 = f_0_vec[2];
    f_0_eq4 = f_0_vec[3];

    /* Update Bulk Scalars and Contravariant velocity */
    // Euler State Update
    for (int j = nGCs; j < nYs - nGCs; j++) {
      for (int i = nGCs; i < nXs - nGCs; i++) {
        rho(i, j) = rho_0(i, j) + dt * f_0_eq1(i, j) / jacobian(i, j);
        xMomentum(i, j) = xMomentum_0(i, j) - dt * f_0_eq2(i, j) / jacobian(i, j);
        yMomentum(i, j) = yMomentum_0(i, j) - dt * f_0_eq3(i, j) / jacobian(i, j);
        rhoE(i, j) = rhoE_0(i, j) + dt * f_0_eq4(i, j) / jacobian(i, j);
      }
    }

    if (iAlt == -1) {
      std::cout << "after solve: " << iProc << " " <<
                dt << " " <<
                f_0_eq1(13, 22) << " " <<
                f_0_eq2(13, 22) << " " <<
                f_0_eq3(13, 22) << " " <<
                f_0_eq4(13, 22) << " " <<
                rhoE(13, 22) << " " << xMomentum(13, 22) << " " << yMomentum(13, 22) << " " <<
                rho_scgc(13, 22, 2) << "\n";
    }


    /* Re-derive Spherical Velocity and Bulk States */
    // Density
    rho_scgc.slice(iAlt) = rho;

    // Bulk Velocity
    xVel = xMomentum / rho; // u^1
    yVel = yMomentum / rho; // u^2
    vel2 = xVel % xVel + yVel % yVel; // Squared Magnitude of Contravariant

    if (iProc == 0 && iAlt == -2) {
      std::cout << "a11 : " << grid.A11_scgc(13, 20, 2) << " "
                << grid.A12_scgc(13, 20, 2) << "\n";
      std::cout << "inv : " << grid.A11_inv_scgc(13, 20, 2) << " "
                << grid.A12_inv_scgc(13, 20, 2) << "\n";
    }

    velocity_vcgc[0].slice(iAlt) = xVel % grid.A11_scgc.slice(
                                     iAlt) + yVel % grid.A12_scgc.slice(iAlt);
    velocity_vcgc[1].slice(iAlt) =
      xVel % grid.A21_scgc.slice(iAlt) +
      yVel % grid.A22_scgc.slice(iAlt);

    /* Update temperature */
    temperature_scgc.slice(iAlt) = (rhoE / rho - 0.5 * vel2) / Cv_scgc.slice(iAlt);


    //if (iAlt == 10) {
    std::cout << "after solve: " << iAlt << " " <<
              temperature_scgc(13, 22, 10) << " " <<
              velocity_vcgc[0](13, 22, 10) << " " << velocity_vcgc[1](13, 22, 10) << " " <<
              rho_scgc(13, 22, 10) << "\n";
    //}

    calc_density_from_mass_concentration();
    //assign_bulk_velocity();

  }

  report.exit(function);
  return;
}

/*
void Neutrals::solver_horizontal_RK1(Grid& grid, Times& time) {
  // Function Reporting
  std::string function = "Neutrals::solver_horizontal_RK1";
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

  iAlt = 0;

    arma_mat x(nXs, nYs), y(nXs, nYs);
    arma_mat jacobian(nXs, nYs), rho(nXs, nYs), rhoE(nXs, nYs), vel2(nXs, nYs);
    arma_mat uVel(nXs, nYs), vVel(nXs, nYs), xVel(nXs, nYs), yVel(nXs, nYs);
    arma_mat xMomentum(nXs, nYs), yMomentum(nXs, nYs);
    arma_mat xMomentum_0(nXs, nYs), yMomentum_0(nXs, nYs);
    arma_mat rho_0(nXs, nYs), rhoE_0(nXs, nYs);
    arma_mat f_0_eq1(nXs, nYs), f_0_eq2(nXs, nYs);
    arma_mat f_0_eq3(nXs, nYs), f_0_eq4(nXs, nYs);

// Advance for bulk calculation first, calculate for every altitude
for (iAlt = nGCs; iAlt < nAlts - nGCs; iAlt++) {
  // Extract Grid Features
arma_mat x = grid.refx_scgc.slice(iAlt);
arma_mat xEdges = grid.refx_Left.slice(iAlt);
arma_mat y = grid.refy_scgc.slice(iAlt);
arma_mat yEdges = grid.refy_Down.slice(iAlt);

// Get reference grid dimensions (Assume dx = dy and equidistant)
arma_vec x_vec = x.col(0);
precision_t dx = x_vec(1) - x_vec(0);
precision_t area = dx * dx;
arma_mat jacobian = grid.sqrt_g_scgc.slice(iAlt);

// States preprocessing
// MASS DENSITY
arma_mat rho = rho_scgc.slice(iAlt);

// VELOCITY
// Get spherical velocity
arma_mat uVel = velocity_vcgc[0].slice(iAlt);
arma_mat vVel = velocity_vcgc[1].slice(iAlt);
// Convert to contravariant (reference) velocity
arma_mat xVel = uVel % grid.A11_inv_scgc.slice(iAlt) + vVel %
                grid.A12_inv_scgc.slice(iAlt); // u^1
arma_mat yVel = uVel % grid.A21_inv_scgc.slice(iAlt) + vVel %
                grid.A22_inv_scgc.slice(iAlt); // u^2
arma_mat vel2 = xVel % xVel + yVel % yVel;
// Generate contravriant momentum (no sqrt(g))
arma_mat xMomentum = rho % xVel; // x1momentum
arma_mat yMomentum = rho % yVel; // x2momentum

// TEMP and ENERGY
// Generate total energy (rhoE (no sqrt(g)))
// (TODO: Verify units)
arma_mat rhoE = rho % (temperature_scgc.slice(iAlt) %
                       Cv_scgc.slice(iAlt) + 0.5 * vel2);

// Advancing with RK4
// Setup Containers
arma_mat rho_0 = rho;
arma_mat xMomentum_0 = xMomentum;
arma_mat yMomentum_0 = yMomentum;
arma_mat rhoE_0 = rhoE;

// FIRST (1) STEP, Compute F_0-> State_1
// Pass in state vector
std::vector<arma_mat> state_0;
state_0.push_back(rho_0);
state_0.push_back(xMomentum_0);
state_0.push_back(yMomentum_0);
state_0.push_back(rhoE_0);
std::vector<arma_mat> f_0_vec = residual_horizontal_rusanov(state_0, grid, time,
                                                            iAlt);
// Extract Gradients
arma_mat f_0_eq1 = f_0_vec[0];
arma_mat f_0_eq2 = f_0_vec[1];
arma_mat f_0_eq3 = f_0_vec[2];
arma_mat f_0_eq4 = f_0_vec[3];

// Update Bulk Scalars and Contravariant velocity
// Euler State Update
for (int j = nGCs; j < nYs - nGCs; j++) {
  for (int i = nGCs; i < nXs - nGCs; i++) {
    rho(i, j) = rho_0(i, j) + dt * f_0_eq1(i, j) / jacobian(i, j);
    xMomentum(i, j) = xMomentum_0(i, j) - dt * f_0_eq2(i, j) / jacobian(i, j);
    yMomentum(i, j) = yMomentum_0(i, j) - dt * f_0_eq3(i, j) / jacobian(i, j);
    rhoE(i, j) = rhoE_0(i, j) + dt * f_0_eq4(i, j) / jacobian(i, j);
  }
}

// Re-derive Spherical Velocity and Bulk States
// Density
rho_scgc.slice(iAlt) = rho;

// Bulk Velocity
xVel = xMomentum / rho; // u^1
yVel = yMomentum / rho; // u^2
vel2 = xVel % xVel + yVel % yVel; // Squared Magnitude of Contravariant
velocity_vcgc[0].slice(iAlt) = xVel % grid.A11_scgc.slice(
                                 iAlt) + yVel % grid.A12_scgc.slice(iAlt);
velocity_vcgc[1].slice(iAlt) = xVel % grid.A21_scgc.slice(
                                 iAlt) + yVel % grid.A22_scgc.slice(iAlt);

// Update specie number density and velocity
for (iSpec = 0; iSpec < nSpecies; iSpec++) {
  //species[iSpec].density_scgc.slice(iAlt) = rho % species[iSpec].mass_concentration_scgc.slice(iAlt);
  species[iSpec].density_scgc.slice(iAlt) = rho / species[iSpec].mass;
  species[iSpec].velocity_vcgc[0].slice(iAlt) = velocity_vcgc[0].slice(iAlt);
  species[iSpec].velocity_vcgc[1].slice(iAlt) = velocity_vcgc[1].slice(iAlt);
}

// Update temperature
temperature_scgc.slice(iAlt) = (rhoE / rho - 0.5 * vel2) / Cv_scgc.slice(iAlt);

report.exit(function);
return;
}

*/

void Neutrals::solver_horizontal_RK4(Grid & grid, Times & time) {
  // Function Reporting
  std::string function = "Neutrals::solver_horizontal_RK4";
  static int iFunction = -1;
  report.enter(function, iFunction);

  // Dimensions of Spatial Discretization
  int64_t nXs = grid.get_nX();
  int64_t nYs = grid.get_nY();
  int64_t nGCs = grid.get_nGCs();
  int64_t nAlts = grid.get_nAlts();
  int iAlt, iSpec;

  // Time Discretization (TODO: change dt calculation method)
  precision_t dt = time.get_dt() / 10;

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
    arma_mat vel2 = xVel % xVel + yVel % yVel;
    // Generate contravriant momentum (no sqrt(g))
    arma_mat xMomentum = rho % xVel; // x1momentum
    arma_mat yMomentum = rho % yVel; // x2momentum

    /* TEMP and ENERGY */
    // Generate total energy (rhoE (no sqrt(g)))
    // (TODO: Verify units)
    //arma_mat rhoE = rho % (temperature_scgc.slice(iAlt) % Cv_scgc.slice(
    //                         iAlt) + 0.5 * vel2);
    arma_mat rhoE = rho % (temperature_scgc.slice(iAlt) %
                           Cv_scgc.slice(iAlt) +
                           0.5 * vel2);

    /** Advancing with RK4 **/
    // Setup Containers
    arma_mat rho_0 = rho;
    arma_mat rho_1(nXs, nYs, fill::zeros); // corresponding f_1
    arma_mat rho_2(nXs, nYs, fill::zeros); // corresponding f_2
    arma_mat rho_3(nXs, nYs, fill::zeros); // corresponding f_3

    arma_mat xMomentum_0 = xMomentum;
    arma_mat xMomentum_1(nXs, nYs, fill::zeros); // corresponding f_1
    arma_mat xMomentum_2(nXs, nYs, fill::zeros); // corresponding f_2
    arma_mat xMomentum_3(nXs, nYs, fill::zeros); // corresponding f_3

    arma_mat yMomentum_0 = yMomentum;
    arma_mat yMomentum_1(nXs, nYs, fill::zeros); // corresponding f_1
    arma_mat yMomentum_2(nXs, nYs, fill::zeros); // corresponding f_2
    arma_mat yMomentum_3(nXs, nYs, fill::zeros); // corresponding f_3

    arma_mat rhoE_0 = rhoE;
    arma_mat rhoE_1(nXs, nYs, fill::zeros); // corresponding f_1
    arma_mat rhoE_2(nXs, nYs, fill::zeros); // corresponding f_2
    arma_mat rhoE_3(nXs, nYs, fill::zeros); // corresponding f_3

    // FIRST (1) STEP, Compute F_0-> State_1
    // Pass in state vector
    std::vector<arma_mat> state_0;
    state_0.push_back(rho_0);
    state_0.push_back(xMomentum_0);
    state_0.push_back(yMomentum_0);
    state_0.push_back(rhoE_0);
    std::vector<arma_mat> f_0_vec = residual_horizontal_rusanov(state_0, grid, time,
                                                                iAlt);
    // Extract Gradients
    arma_mat f_0_eq1 = f_0_vec[0];
    arma_mat f_0_eq2 = f_0_vec[1];
    arma_mat f_0_eq3 = f_0_vec[2];
    arma_mat f_0_eq4 = f_0_vec[3];

    /* Update Bulk Scalars and Contravariant velocity */
    // Euler State Update
    for (int j = nGCs; j < nYs - nGCs; j++) {
      for (int i = nGCs; i < nXs - nGCs; i++) {
        rho_1(i, j) = rho_0(i, j) + 0.5 * dt * f_0_eq1(i, j) / jacobian(i, j);
        xMomentum_1(i, j) = xMomentum_0(i, j) + 0.5 * dt * f_0_eq2(i, j) / jacobian(i,
                            j);
        yMomentum_1(i, j) = yMomentum_0(i, j) + 0.5 * dt * f_0_eq3(i, j) / jacobian(i,
                            j);
        rhoE_1(i, j) = rhoE_0(i, j) + 0.5 * dt * f_0_eq4(i, j) / jacobian(i, j);
      }
    }

    // SECOND (2) STEP, Compute F_1-> State_2
    // Pass in state vector
    std::vector<arma_mat> state_1;
    state_1.push_back(rho_1);
    state_1.push_back(xMomentum_1);
    state_1.push_back(yMomentum_1);
    state_1.push_back(rhoE_1);
    std::vector<arma_mat> f_1_vec = residual_horizontal_rusanov(state_1, grid, time,
                                                                iAlt);
    // Extract Gradients
    arma_mat f_1_eq1 = f_1_vec[0];
    arma_mat f_1_eq2 = f_1_vec[1];
    arma_mat f_1_eq3 = f_1_vec[2];
    arma_mat f_1_eq4 = f_1_vec[3];

    /* Update Bulk Scalars and Contravariant velocity */
    // Euler State Update
    for (int j = nGCs; j < nYs - nGCs; j++) {
      for (int i = nGCs; i < nXs - nGCs; i++) {
        rho_2(i, j) = rho_0(i, j) + 0.5 * dt * f_1_eq1(i, j) / jacobian(i, j);
        xMomentum_2(i, j) = xMomentum_0(i, j) + 0.5 * dt * f_1_eq2(i, j) / jacobian(i,
                            j);
        yMomentum_2(i, j) = yMomentum_0(i, j) + 0.5 * dt * f_1_eq3(i, j) / jacobian(i,
                            j);
        rhoE_2(i, j) = rhoE_0(i, j) + 0.5 * dt * f_1_eq4(i, j) / jacobian(i, j);
      }
    }

    // THIRD (3) STEP, Compute F_2-> State_3
    // Pass in state vector
    std::vector<arma_mat> state_2;
    state_2.push_back(rho_2);
    state_2.push_back(xMomentum_2);
    state_2.push_back(yMomentum_2);
    state_2.push_back(rhoE_2);
    std::vector<arma_mat> f_2_vec = residual_horizontal_rusanov(state_2, grid, time,
                                                                iAlt);
    // Extract Gradients
    arma_mat f_2_eq1 = f_2_vec[0];
    arma_mat f_2_eq2 = f_2_vec[1];
    arma_mat f_2_eq3 = f_2_vec[2];
    arma_mat f_2_eq4 = f_2_vec[3];

    /* Update Bulk Scalars and Contravariant velocity */
    // Euler State Update
    for (int j = nGCs; j < nYs - nGCs; j++) {
      for (int i = nGCs; i < nXs - nGCs; i++) {
        rho_3(i, j) = rho_0(i, j) + dt * f_2_eq1(i, j) / jacobian(i, j);
        xMomentum_3(i, j) = xMomentum_0(i, j) + dt * f_2_eq2(i, j) / jacobian(i, j);
        yMomentum_3(i, j) = yMomentum_0(i, j) + dt * f_2_eq3(i, j) / jacobian(i, j);
        rhoE_3(i, j) = rhoE_0(i, j) + dt * f_2_eq4(i, j) / jacobian(i, j);
      }
    }

    // FOURTH (4) STEP, Compute F_3
    // Pass in state vector
    std::vector<arma_mat> state_3;
    state_3.push_back(rho_3);
    state_3.push_back(xMomentum_3);
    state_3.push_back(yMomentum_3);
    state_3.push_back(rhoE_3);
    std::vector<arma_mat> f_3_vec = residual_horizontal_rusanov(state_3, grid, time,
                                                                iAlt);
    // Extract Gradients
    arma_mat f_3_eq1 = f_3_vec[0];
    arma_mat f_3_eq2 = f_3_vec[1];
    arma_mat f_3_eq3 = f_3_vec[2];
    arma_mat f_3_eq4 = f_3_vec[3];

    // Summing all steps for final update
    arma_mat f_sum_eq1 = f_0_eq1 + 2 * f_1_eq1 + 2 * f_2_eq1 + f_3_eq1;
    arma_mat f_sum_eq2 = f_0_eq2 + 2 * f_1_eq2 + 2 * f_2_eq2 + f_3_eq2;
    arma_mat f_sum_eq3 = f_0_eq3 + 2 * f_1_eq3 + 2 * f_2_eq3 + f_3_eq3;
    arma_mat f_sum_eq4 = f_0_eq4 + 2 * f_1_eq4 + 2 * f_2_eq4 + f_3_eq4;

    /* Update Bulk Scalars and Contravariant velocity */
    // Euler State Update
    for (int j = nGCs; j < nYs - nGCs; j++) {
      for (int i = nGCs; i < nXs - nGCs; i++) {
        rho(i, j) = rho(i, j) + dt / 6 * f_sum_eq1(i, j) / jacobian(i, j);
        xMomentum(i, j) = xMomentum(i, j) + dt / 6 * f_sum_eq2(i, j) / jacobian(i, j);
        yMomentum(i, j) = yMomentum(i, j) + dt / 6 * f_sum_eq3(i, j) / jacobian(i, j);
        rhoE(i, j) = rhoE(i, j) + dt / 6 * f_sum_eq4(i, j) / jacobian(i, j);
      }
    }

    /* Re-derive Spherical Velocity and Bulk States */
    // Density
    rho_scgc.slice(iAlt) = rho;

    // Bulk Velocity
    xVel = xMomentum / rho; // u^1
    yVel = yMomentum / rho; // u^2
    vel2 = xVel % xVel + yVel % yVel; // Squared Magnitude of Contravariant
    velocity_vcgc[0].slice(iAlt) = xVel % grid.A11_scgc.slice(
                                     iAlt) + yVel % grid.A12_scgc.slice(iAlt);
    velocity_vcgc[1].slice(iAlt) = xVel % grid.A21_scgc.slice(
                                     iAlt) + yVel % grid.A22_scgc.slice(iAlt);

    /* Update specie number density and velocity */
    for (iSpec = 0; iSpec < nSpecies; iSpec++) {
      //species[iSpec].density_scgc.slice(iAlt) = rho % species[iSpec].concentration_scgc.slice(iAlt);
      species[iSpec].density_scgc.slice(iAlt) = rho / species[iSpec].mass;
      species[iSpec].velocity_vcgc[0].slice(iAlt) = velocity_vcgc[0].slice(iAlt);
      species[iSpec].velocity_vcgc[1].slice(iAlt) = velocity_vcgc[1].slice(iAlt);
    }

    /* Update temperature */
    temperature_scgc.slice(iAlt) = (rhoE / rho - 0.5 * vel2) / Cv_scgc.slice(iAlt);

    report.exit(function);
    return;
  }
}
