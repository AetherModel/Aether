// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "../include/aether.h"

// --------------------------------------------------------------------------
// Calculate the electric field from the potential
// --------------------------------------------------------------------------

void Ions::calc_efield(Grid grid, Report &report) {

  // efield = - grad(potential)
  efield_vcgc = calc_gradient_vector(potential_scgc, grid);

  for (int64_t iComp = 0; iComp < 3; iComp++)
    efield_vcgc[iComp] = -efield_vcgc[iComp];

  // Remove component along b-field (should be zero, anyways!)
  fcube edotb = dot_product(efield_vcgc, grid.bfield_unit_vcgc);

  for (int64_t iComp = 0; iComp < 3; iComp++)
    efield_vcgc[iComp] =
      efield_vcgc[iComp] - edotb % grid.bfield_unit_vcgc[iComp];
}

// --------------------------------------------------------------------------
// Calculate the E x B drift from the electric field and magnetic field
// --------------------------------------------------------------------------

void Ions::calc_exb_drift(Grid grid, Report &report) {
  fcube bmag2 =
    (grid.bfield_mag_scgc) % (grid.bfield_mag_scgc);
  exb_vcgc = cross_product(efield_vcgc, grid.bfield_vcgc);

  for (int64_t iComp = 0; iComp < 3; iComp++)
    exb_vcgc[iComp] = exb_vcgc[iComp] / bmag2;
}

// --------------------------------------------------------------------------
//
// --------------------------------------------------------------------------

std::vector<fcube> Ions::calc_ion_electron_pressure_gradient(int64_t iIon,
    Grid grid,
    Report &report) {
  std::vector<fcube> pressure_gradient_vcgc;
  fcube total_pressure_scgc;

  // Total Pressure =
  //     Ion Pressure + Electron Pressure =
  //     (Ni * Ti + Ne * Te) * k

  total_pressure_scgc =
    (species[iIon].density_scgc %
     species[iIon].temperature_scgc +
     density_scgc %
     electron_temperature_scgc) *
    cKB;

  pressure_gradient_vcgc = calc_gradient_vector(total_pressure_scgc, grid);

  return pressure_gradient_vcgc;
}

// --------------------------------------------------------------------------
// Calculate the ion drift
// --------------------------------------------------------------------------

void Ions::calc_ion_drift(Neutrals neutrals,
                          Grid grid,
                          float dt,
                          Report &report) {

  std::string function = "Ions::calc_ion_drift";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int64_t nX = grid.get_nX();
  int64_t nY = grid.get_nY();
  int64_t nZ = grid.get_nZ();

  report.print(5, "going into calc_efield");
  calc_efield(grid, report);

  // This is for the electron drift motion:
  report.print(5, "going into calc_exb_drift");
  calc_exb_drift(grid, report);

  std::vector<fcube> gravity_vcgc = make_cube_vector(nX, nY, nZ, 3);
  std::vector<fcube> wind_forcing = make_cube_vector(nX, nY, nZ, 3);
  std::vector<fcube> total_forcing = make_cube_vector(nX, nY, nZ, 3);

  report.print(5, "going into collision frequencies");
  calc_ion_neutral_coll_freq(neutrals, report);

  int64_t iIon, iNeutral;

  std::vector<fcube> grad_Pi_plus_Pe;
  fcube rho, rho_nuin, nuin_sum, Nie, sum_rho;

  nuin_sum.set_size(nX, nY, nZ);
  nuin_sum.zeros();

  sum_rho.set_size(nX, nY, nZ);
  sum_rho.zeros();

  fill_electrons(report);

  for (int64_t iComp = 0; iComp < 3; iComp++)
    velocity_vcgc[iComp].zeros();

  for (iIon = 0; iIon < nIons; iIon++) {

    for (int64_t iComp = 0; iComp < 3; iComp++) {
      species[iIon].perp_velocity_vcgc[iComp].zeros();
      species[iIon].par_velocity_vcgc[iComp].zeros();
    }

    if (species[iIon].DoAdvect) {

      // Need mass density for the current ion species:
      rho = species[iIon].mass * species[iIon].density_scgc;
      Nie = cE * species[iIon].density_scgc;

      // Get gradient in pressure:
      report.print(5, "going into pressure gradient");
      grad_Pi_plus_Pe = calc_ion_electron_pressure_gradient(iIon, grid, report);

      // This is assuming that the 3rd dim is radial.
      // Want actual gravity for 3rd dim
      // (negative is due to gravity being positive for some reason):
      gravity_vcgc[2] = -grid.gravity_scgc % rho;

      // Neutral Wind Forcing:
      report.print(5, "neutral winds");

      for (int64_t iComp = 0; iComp < 3; iComp++)
        wind_forcing[iComp].zeros();

      for (iNeutral = 0; iNeutral < nSpecies; iNeutral++) {
        rho_nuin = rho % species[iIon].nu_ion_neutral_vcgc[iNeutral];
        nuin_sum = nuin_sum + species[iIon].nu_ion_neutral_vcgc[iNeutral];

        for (int64_t iComp = 0; iComp < 3; iComp++) {
          wind_forcing[iComp] = wind_forcing[iComp] +
                                rho_nuin % neutrals.velocity_vcgc[iComp];
        }
      }

      // Total Forcing (sum everything - this is A_s):
      for (int64_t iComp = 0; iComp < 3; iComp++) {
        total_forcing[iComp] =
          - grad_Pi_plus_Pe[iComp]
          + gravity_vcgc[iComp]
          + wind_forcing[iComp]
          + Nie % efield_vcgc[iComp];
      }

      std::vector<fcube> a_par = make_cube_vector(nX, nY, nZ, 3);
      std::vector<fcube> a_perp = make_cube_vector(nX, nY, nZ, 3);
      std::vector<fcube> a_x_b;

      if (grid.get_HasBField()) {
        // With a Planetary Magnetic field
        fcube a_dot_b = dot_product(total_forcing, grid.bfield_unit_vcgc);

        for (int64_t iComp = 0; iComp < 3; iComp++) {
          a_par[iComp] = a_dot_b % grid.bfield_unit_vcgc[iComp];
          a_perp[iComp] = total_forcing[iComp] - a_par[iComp];
        }

        a_x_b = cross_product(a_perp, grid.bfield_vcgc);

        for (int64_t iComp = 0; iComp < 3; iComp++)
          species[iIon].perp_velocity_vcgc[iComp] =
            (rho_nuin % a_perp[iComp] + Nie % a_x_b[iComp]) /
            (rho_nuin % rho_nuin +
             Nie % Nie % grid.bfield_mag_scgc % grid.bfield_mag_scgc);
      } else {
        // No Planetary Magnetic field
        for (int64_t iComp = 0; iComp < 3; iComp++) {
          a_par[iComp] = total_forcing[iComp];
          // Steady state:
          species[iIon].par_velocity_vcgc[iComp] =
            a_par[iComp] / rho / nuin_sum;
          //species[iIon].par_velocity_vcgc[iComp] =
          //  (species[iIon].par_velocity_vcgc[iComp] + a_par[iComp] * dt / rho) /
          //  (1 + nuin_sum * dt);
        }
      }

      // Calculate the mass-weighted average total velocity
      sum_rho = sum_rho + rho;

      for (int64_t iComp = 0; iComp < 3; iComp++) {
        velocity_vcgc[iComp] = velocity_vcgc[iComp] +
                               rho % (species[iIon].perp_velocity_vcgc[iComp] +
                                      species[iIon].par_velocity_vcgc[iComp]);
      }

    }  // if DoAdvect

  }  // for iIon

  for (int64_t iComp = 0; iComp < 3; iComp++)
    velocity_vcgc[iComp] = velocity_vcgc[iComp] / sum_rho;

  report.exit(function);
  return;
}


