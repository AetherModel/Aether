// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <cmath>
#include <iostream>

#include "aether.h"

// ----------------------------------------------------------------------
//  Calculate mass density and number density:
// ----------------------------------------------------------------------

void Neutrals::calc_mass_density(Report &report) {
    
  std::string function = "Neutrals::calc_mass_density";
  static int iFunction = -1;
  report.enter(function, iFunction);
    
  rho_scgc.zeros();
  density_scgc.zeros();
    
  for (int64_t iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    rho_scgc = rho_scgc +
      species[iSpecies].mass * species[iSpecies].density_scgc;
    density_scgc = density_scgc + species[iSpecies].density_scgc;
  }
  report.exit(function);
}

// ----------------------------------------------------------------------
//  Calculate concentration (species number density / total density)
//    Must call calc_mass_density first!!!
// ----------------------------------------------------------------------

void Neutrals::calc_concentration(Report &report) {
    
  std::string function = "Neutrals::calc_concentration";
  static int iFunction = -1;
  report.enter(function, iFunction);
    
  for (int64_t iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    species[iSpecies].concentration_scgc =
      species[iSpecies].density_scgc / density_scgc;
  report.exit(function);  
}

// ----------------------------------------------------------------------
//  Calculate mean major mass
//     Must call calc_mass_density first!!!
// ----------------------------------------------------------------------

void Neutrals::calc_mean_major_mass(Report &report) {
    
  std::string function = "Neutrals::calc_mean_major_mass";
  static int iFunction = -1;
  report.enter(function, iFunction);

  mean_major_mass_scgc = rho_scgc / density_scgc;  
  report.exit(function);
}

// ----------------------------------------------------------------------
//  Calculate pressure
//     Must call calc_mass_density first!!!
// ----------------------------------------------------------------------

void Neutrals::calc_pressure(Report &report) {
    
  std::string function = "Neutrals::calc_pressure";
  static int iFunction = -1;
  report.enter(function, iFunction);

  pressure_scgc = cKB * density_scgc % temperature_scgc;
  report.exit(function);
}

// ----------------------------------------------------------------------
//  Calculate bulk velocity:
// ((sum of species mass * density * velocity) / bulk mass density)
// ----------------------------------------------------------------------

void Neutrals::calc_bulk_velocity(Report &report) {
    
  std::string function = "Neutrals::calc_bulk_velocity";
  static int iFunction = -1;
  report.enter(function, iFunction);
    
  for (int64_t iDir = 0; iDir < 3; iDir++) {
    velocity_vcgc[iDir].zeros();
    for (int64_t iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      velocity_vcgc[iDir] +=
	species[iSpecies].mass *
	species[iSpecies].density_scgc %
	species[iSpecies].velocity_vcgc[iDir];
    velocity_vcgc[iDir] = velocity_vcgc[iDir] / rho_scgc;
  }

  report.exit(function);
}

// ----------------------------------------------------------------------
// Calculate a bunch of derived products:
//   - Specific Heat at Constant Volume (Cv)
//   - Gamma
//   - Kappa
//   - Speed of sound
// ----------------------------------------------------------------------

void Neutrals::calc_specific_heat(Report &report) {

  int64_t iSpecies;

  std::string function = "Neutrals::calc_specific_heat";
  static int iFunction = -1;
  report.enter(function, iFunction);

  Cv_scgc.zeros();
  gamma_scgc.zeros();
  kappa_scgc.zeros();

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    Cv_scgc = Cv_scgc +
              (species[iSpecies].vibe - 2) *
              species[iSpecies].density_scgc *
              cKB / species[iSpecies].mass;
    gamma_scgc = gamma_scgc +
                 species[iSpecies].density_scgc / (species[iSpecies].vibe - 2);
    kappa_scgc = kappa_scgc +
                 species[iSpecies].thermal_cond *
                 species[iSpecies].density_scgc %
                 pow(temperature_scgc, species[iSpecies].thermal_exp);
  }

  Cv_scgc = Cv_scgc / (2 * density_scgc);
  gamma_scgc = gamma_scgc * 2.0 / density_scgc + 1.0;
  kappa_scgc = kappa_scgc / density_scgc;

  sound_scgc = sqrt(cKB *
                    gamma_scgc %
                    temperature_scgc /
                    mean_major_mass_scgc);

  report.exit(function);
  return;
}

//----------------------------------------------------------------------
// Calculate the altitude integral of the different species for EUV
// calculations.  Scale these to be the slant path through the
// atmosphere.  This gets complicated behind the terminator.  All of
// this is taken from Smith and Smith, JGR 1972, vol. 77, page 3592
// ----------------------------------------------------------------------

void Neutrals::calc_chapman(Grid grid, Report &report) {

  int64_t iAlt, iLon, iLat;

  // This is all from Smith and Smith, JGR 1972, vol. 77, page 3592
  // "Numerical evaluation of chapman's grazing incidence integral ch(X,x)"
  // Xp is supposed to be R/H
  // JMB Update: 05/2017.  Corrected a small error in the y-calc for
  // erfc(y)
  //
  // Also Updated the Grazing Integral for SZA > 90.0
  // We now do log-linear interpolation for smoother transitions

  double a = 1.06069630;
  double b = 0.55643831;
  double c = 1.06198960;
  double d = 1.72456090;
  double f = 0.56498823;
  double g = 0.06651874;

  double y, dy;

  precision_t grad_xp, grad_in, Xg, in, int_g, int_p;
  int64_t iiAlt;

  std::string function = "Neutrals::calc_chapman";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int64_t nLons = grid.get_nLons();
  int64_t nLats = grid.get_nLats();
  int64_t nAlts = grid.get_nAlts();
  int64_t nGCs = grid.get_nGCs();

  // New way of doing it with 3D arrays:

  arma_cube integral3d(nLons, nLats, nAlts);
  arma_cube log_int3d(nLons, nLats, nAlts);
  arma_cube xp3d(nLons, nLats, nAlts);
  arma_cube y3d(nLons, nLats, nAlts);
  arma_cube erfcy3d(nLons, nLats, nAlts);

  arma_vec integral1d(nAlts);
  arma_vec log_int1d(nAlts);
  arma_vec xp1d(nAlts);
  arma_vec y1d(nAlts);
  arma_vec erfcy1d(nAlts);
  arma_vec dAlt1d(nAlts);
  arma_vec sza1d(nAlts);
  arma_vec radius1d(nAlts);
  arma_vec H1d(nAlts);

  for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++) {

    species[iSpecies].scale_height_scgc =
      cKB * temperature_scgc /
      (species[iSpecies].mass * grid.gravity_scgc);

    xp3d = grid.radius_scgc / species[iSpecies].scale_height_scgc;
    y3d = sqrt(0.5 * xp3d) % abs(grid.cos_sza_scgc);
    iAlt = nAlts - 1;

    integral3d.fill(0.0);

    integral3d.slice(iAlt) =
      species[iSpecies].density_scgc.slice(iAlt) %
      species[iSpecies].scale_height_scgc.slice(iAlt);

    for (iAlt = nAlts - 1; iAlt >= 0; iAlt--) {
      if (iAlt < nAlts - 1) {
        integral3d.slice(iAlt) = integral3d.slice(iAlt + 1) +
                                 species[iSpecies].density_scgc.slice(iAlt) %
                                 grid.dalt_lower_scgc.slice(iAlt + 1);
      }
    }

    species[iSpecies].rho_alt_int_scgc = integral3d * species[iSpecies].mass;

    erfcy3d = (a + b * y3d) / (c + d * y3d + y3d % y3d);

    for (iLon = 0; iLon < nLons ; iLon++)
      for (iLat = 0; iLat < nLats ; iLat++)
        for (iAlt = 0; iAlt < nAlts ; iAlt++)
          if (y3d(iLon, iLat, iAlt) >= 8.0)
            erfcy3d(iLon, iLat, iAlt) = f / (g + y3d(iLon, iLat, iAlt));

    log_int3d = log(integral3d);

    // Set chapman integrals to max in the lower ghostcells

    species[iSpecies].chapman_scgc.fill(max_chapman);

    for (iLon = 0; iLon < nLons ; iLon++) {
      for (iLat = 0; iLat < nLats ; iLat++) {

        dAlt1d = grid.dalt_lower_scgc.tube(iLon, iLat);
        sza1d = grid.sza_scgc.tube(iLon, iLat);
        integral1d = integral3d.tube(iLon, iLat);
        log_int1d = log_int3d.tube(iLon, iLat);
        xp1d = xp3d.tube(iLon, iLat);
        y1d = y3d.tube(iLon, iLat);
        erfcy1d = erfcy3d.tube(iLon, iLat);
        radius1d = grid.radius_scgc.tube(iLon, iLat);
        H1d = species[iSpecies].scale_height_scgc.tube(iLon, iLat);

        for (iAlt = nGCs; iAlt < nAlts; iAlt++) {
          // This is on the dayside:
          if (sza1d(iAlt) < cPI / 2 || sza1d(iAlt) > 3 * cPI / 2) {
            species[iSpecies].chapman_scgc(iLon, iLat, iAlt) =
              integral1d(iAlt) * sqrt(0.5 * cPI * xp1d(iAlt)) * erfcy1d(iAlt);
          } else {
            // This is on the nghtside of the terminator:

            y = radius1d(iAlt) * abs(cos(sza1d(iAlt) - cPI / 2));

            // This sort of assumes that nGeoGhosts >= 2:
            if (y > radius1d(nGCs)) {

              iiAlt = iAlt;

              while (radius1d(iiAlt - 1) > y)
                iiAlt--;

              iiAlt--;

              // make sure to use the proper cell spacing (iiAlt+1 & lower):
              grad_xp = (xp1d(iiAlt + 1) - xp1d(iiAlt)) / dAlt1d(iiAlt + 1);
              grad_in = (log_int1d(iiAlt + 1) - log_int1d(iiAlt)) /
                        dAlt1d(iiAlt + 1);

              // Linearly interpolate H and X:
              dy = y - radius1d(iiAlt);
              Xg = xp1d(iiAlt) + grad_xp * dy;
              in = log_int1d(iiAlt) + grad_in * dy;

              int_g = exp(in);
              int_p = integral1d(iAlt);
              // Equation (19) Smith & Smith
              species[iSpecies].chapman_scgc(iLon, iLat, iAlt) =
                sqrt(0.5 * cPI * Xg) * (2.0 * int_g - int_p * erfcy1d(iAlt));

              if (species[iSpecies].chapman_scgc(iLon, iLat, iAlt) >
                  max_chapman)
                species[iSpecies].chapman_scgc(iLon, iLat, iAlt) = max_chapman;

            } else {
              // This says that we are in the shadow of the planet:

              species[iSpecies].chapman_scgc(iLon, iLat, iAlt) = max_chapman;
            }
          }
        }  // iAlt
      }  // iLat
    }  // iLon
  }  // iSpecies

  report.exit(function);
  return;
}

// -----------------------------------------------------------------------------
// Calculate thermal conduction
// -----------------------------------------------------------------------------

void Neutrals::calc_conduction(Grid grid, Times time, Report &report) {

  precision_t dt;

  int64_t iLon, iLat;

  std::string function = "Neutrals::calc_conduction";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int64_t nLons = grid.get_nLons();
  int64_t nLats = grid.get_nLats();
  int64_t nAlts = grid.get_nAlts();

  arma_cube rhocvr23d(nLons, nLats, nAlts);
  arma_cube lambda3d(nLons, nLats, nAlts);
  arma_cube prandtl3d(nLons, nLats, nAlts);

  rhocvr23d = rho_scgc % Cv_scgc % grid.radius2_scgc;
  // Need to make this eddy * rho * cv:
  prandtl3d.zeros();
  lambda3d = (kappa_scgc + prandtl3d) % grid.radius2_scgc;

  arma_vec temp1d(nAlts);
  arma_vec lambda1d(nAlts);
  arma_vec rhocvr21d(nAlts);
  arma_vec dalt1d(nAlts);
  arma_vec conduction1d(nAlts);

  for (iLon = 0; iLon < nLons; iLon++) {
    for (iLat = 0; iLat < nLats; iLat++) {

      temp1d = temperature_scgc.tube(iLon, iLat);
      lambda1d = lambda3d.tube(iLon, iLat);
      rhocvr21d = rhocvr23d.tube(iLon, iLat);
      dalt1d = grid.dalt_lower_scgc.tube(iLon, iLat);
      conduction1d.zeros();

      dt = time.get_dt();

      conduction1d = solver_conduction(temp1d, lambda1d, rhocvr21d, dt, dalt1d);

      // We want the sources to be in terms of dT/dt, while the
      // conduction actually solves for Tnew-Told, so divide by dt

      conduction_scgc.tube(iLon, iLat) = conduction1d / dt;
    }  // lat
  }  // lon

  report.exit(function);
}
