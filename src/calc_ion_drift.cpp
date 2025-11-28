// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "../include/aether.h"

// --------------------------------------------------------------------------
// Calculate the electric field from the potential
// --------------------------------------------------------------------------

void Ions::calc_efield(Grid &grid) {

  // efield = - grad(potential)
  efield_vcgc = calc_gradient_vector(-1.0 * potential_scgc, grid);

  // Remove component along b-field (should be zero, anyways!)
  arma_cube edotb;
  edotb = dot_product(efield_vcgc, grid.bfield_unit_vcgc);

  for (int64_t iComp = 0; iComp < 3; iComp++)
    efield_vcgc[iComp] =
      efield_vcgc[iComp] - edotb % grid.bfield_unit_vcgc[iComp];
}

// --------------------------------------------------------------------------
// Calculate the E x B drift from the electric field and magnetic field
// --------------------------------------------------------------------------

void Ions::calc_exb_drift(Grid &grid) {
  std::string function = "Ions::calc_exb";
  static int iFunction = -1;
  report.enter(function, iFunction);

  arma_cube bmag2 =
    (grid.bfield_mag_scgc) % (grid.bfield_mag_scgc);
  exb_vcgc = cross_product(efield_vcgc, grid.bfield_vcgc);

  for (int64_t iComp = 0; iComp < 3; iComp++)
    exb_vcgc[iComp] = exb_vcgc[iComp] / bmag2;

  report.exit(function);
}

// --------------------------------------------------------------------------
// Calculate the ion + electron pressure for the specified ion species
// --------------------------------------------------------------------------

std::vector<arma_cube> Ions::calc_ion_electron_pressure_gradient(int64_t iIon,
    Grid grid) {

  std::string function = "Ions::elec_ion_pressure_gradient";
  static int iFunction = -1;
  report.enter(function, iFunction);
  std::vector<arma_cube> pressure_gradient_vcgc;
  arma_cube total_pressure_scgc;

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
  report.exit(function);

  return pressure_gradient_vcgc;
}

// --------------------------------------------------------------------------
// Calculate the ion drift
// --------------------------------------------------------------------------

void Ions::calc_ion_drift(Neutrals &neutrals,
                          Grid &grid,
                          precision_t dt) {

  std::string function = "Ions::calc_ion_drift";
  static int iFunction = -1;
  report.enter(function, iFunction);

  // CHANGE !!!
  dt = dt / 10.0;

  int64_t nX = grid.get_nX();
  int64_t nY = grid.get_nY();
  int64_t nZ = grid.get_nZ();

  set_floor();

  report.print(5, "going into calc_efield");
  calc_efield(grid);

  // This is for the electron drift motion:
  report.print(5, "going into calc_exb_drift");
  calc_exb_drift(grid);

  int64_t iIon, iNeutral, iDim;
  int64_t iComp;

  nuin_sum.zeros();
  sum_rho.zeros();

  fill_electrons();

  for (iComp = 0; iComp < 3; iComp++)
    velocity_vcgc[iComp].zeros();

  for (iIon = 0; iIon < nSpecies; iIon++) {

    for (iComp = 0; iComp < 3; iComp++)
      species[iIon].perp_velocity_vcgc[iComp].zeros();

    if (species[iIon].DoAdvect) {

      // Need mass density for the current ion species:
      rho = species[iIon].mass * species[iIon].density_scgc;
      Nie = cE * species[iIon].density_scgc;

      // Get gradient in pressure:
      report.print(5, "going into pressure gradient");
      grad_Pi_plus_Pe = calc_ion_electron_pressure_gradient(iIon, grid);

      // This is assuming that the 3rd dim is radial.
      // Want actual gravity for 3rd dim
      for (iComp = 0; iComp < 3; iComp ++) {
        gravity_vcgc[iComp] = grid.gravity_vcgc[iComp];
        grad_Pi_plus_Pe[iComp] = grad_Pi_plus_Pe[iComp] / rho;
        efield_acc[iComp] = Nie % efield_vcgc[iComp] / rho;
      }

      // Neutral Wind Forcing:
      report.print(5, "neutral winds");

      for (iComp = 0; iComp < 3; iComp++)
        wind_acc[iComp].zeros();

      nuin_sum.zeros();

      for (iNeutral = 0; iNeutral < neutrals.nSpecies; iNeutral++) {
        nuin = species[iIon].nu_ion_neutral_vcgc[iNeutral];
        nuin_sum = nuin_sum + species[iIon].nu_ion_neutral_vcgc[iNeutral];

        for (iComp = 0; iComp < 3; iComp++) {
          wind_acc[iComp] = wind_acc[iComp] +
                            nuin % neutrals.velocity_vcgc[iComp];
        }
      }

      // Total Forcing (sum everything - this is A_s):
      for (iComp = 0; iComp < 3; iComp++) {
        total_acc[iComp] =
          - grad_Pi_plus_Pe[iComp]
          + gravity_vcgc[iComp]
          + wind_acc[iComp]
          + efield_acc[iComp];
      }

      if (grid.get_HasBField()) {
        // With a Planetary Magnetic field
        arma_cube a_dot_b = dot_product(total_acc, grid.bfield_unit_vcgc);

        for (iComp = 0; iComp < 3; iComp++) {
          a_par[iComp] = a_dot_b % grid.bfield_unit_vcgc[iComp];
          a_perp[iComp] = total_acc[iComp] - a_par[iComp];
        }

        a_x_b = cross_product(a_perp, grid.bfield_vcgc);

        // With floats, this can become 0, which then makes the
        // velocity a nan, so the clamp ensures that the bottom is not 0
        bottom =
          rho % rho % nuin % nuin +
          Nie % Nie % grid.bfield_mag_scgc % grid.bfield_mag_scgc;
        bottom.clamp(1e-32, 1e32);

        for (iComp = 0; iComp < 3; iComp++) {
          // I redefined A to be an acceleration instead of a force, which
          // then changes the definition of top
          top = rho % nuin % a_perp[iComp] + Nie % a_x_b[iComp];
          species[iIon].perp_velocity_vcgc[iComp] = rho % top / bottom;

          // Steady state:
          //species[iIon].par_velocity_vcgc[iComp] =
          //  a_par[iComp] / rho / nuin_sum;
          species[iIon].par_velocity_vcgc[iComp] =
            (species[iIon].par_velocity_vcgc[iComp] + a_par[iComp] * dt) /
            (1 + nuin_sum * dt);

          // These need to change, since they are dependent on the
          // grid. Closed, dipole fieldlines should NOT do this!!!
          species[iIon].par_velocity_vcgc[iComp].slice(nZ - 1).zeros();
          species[iIon].par_velocity_vcgc[iComp].slice(nZ - 2).zeros();
          species[iIon].par_velocity_vcgc[iComp].slice(nZ - 3) =
            species[iIon].par_velocity_vcgc[iComp].slice(nZ - 4);
          species[iIon].par_velocity_vcgc[iComp].clamp(-100, 100);

        }
      } else {
        // No Planetary Magnetic field
        for (iComp = 0; iComp < 3; iComp++) {
          a_par[iComp] = total_acc[iComp];
          // Steady state:
          //species[iIon].par_velocity_vcgc[iComp] =
          //  a_par[iComp] / rho / nuin_sum;
          species[iIon].par_velocity_vcgc[iComp] =
            (species[iIon].par_velocity_vcgc[iComp] + a_par[iComp] * dt / rho) /
            (1 + nuin_sum * dt);
          species[iIon].par_velocity_vcgc[iComp].clamp(-100, 100);

        }
      }

      // Calculate the mass-weighted average total velocity
      sum_rho = sum_rho + rho;

      for (iComp = 0; iComp < 3; iComp++) {
        species[iIon].velocity_vcgc[iComp] =
          species[iIon].perp_velocity_vcgc[iComp] +
          species[iIon].par_velocity_vcgc[iComp];
        velocity_vcgc[iComp] = velocity_vcgc[iComp] +
                               rho % (species[iIon].velocity_vcgc[iComp]);
      }

    }  // if DoAdvect

  }  // for iIon

  // This is the mass weighted total bulk velocity:
  for (iComp = 0; iComp < 3; iComp++)
    velocity_vcgc[iComp] = velocity_vcgc[iComp] / sum_rho;

  report.exit(function);
  return;
}


