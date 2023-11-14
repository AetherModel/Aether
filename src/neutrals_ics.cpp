// Copyright 2023, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md
//
// initial version - A. Ridley - May 27, 2023

#include "aether.h"

// -----------------------------------------------------------------------------
//  Set initial conditions for the neutrals.
//    Two methods implemented so far:
//      - Planet: Use fixed density values in the planet.in file and the
//                temperature profile to set the densities and temperature.
//                Densities are filled with hydrostatic solution.
//      - Msis: Use NRL MSIS to set the densities and temperatures.  If the
//              densities are not found, then set to density in planet.in
//              file and fill with hydrostatic.
// -----------------------------------------------------------------------------

bool Neutrals::initial_conditions(Grid grid,
                                  Times time,
                                  Indices indices) {

  std::string function = "Neutrals::initial_conditions";
  static int iFunction = -1;
  report.enter(function, iFunction);

  // initialize didWork to false, so that we can catch if the initial
  // conditions are not actually set
  bool didWork = false;

  int64_t iLon, iLat, iAlt, iA;
  precision_t alt, r;
  int64_t nAlts = grid.get_nZ();

  report.print(3, "Creating Neutrals initial_condition");

  if (input.get_do_restart()) {
    report.print(1, "Restarting! Reading neutral files!");
    didWork = restart_file(input.get_restartin_dir(), DoRead);

    if (!didWork) {
      report.error("Reading Restart for Neutrals Failed!!!");
      \
    }
  } else {

    json ics = input.get_initial_condition_types();
    std::string icsType = mklower(ics["type"]);

    if (icsType == "msis") {

      report.print(2, "Using MSIS for Initial Conditions");

      Msis msis;

      if (!msis.is_ok()) {
        didWork = false;
        report.error("MSIS Initial Conditions asked for, ");
        report.error("but MSIS is not compiled/working!");
      } else {
        didWork = true;
        msis.set_time(time);
        precision_t f107 = indices.get_f107(time.get_current());
        precision_t f107a = indices.get_f107a(time.get_current());
        msis.set_f107(f107, f107a);
        msis.set_ap(10.0);
        msis.set_locations(grid.geoLon_scgc,
                           grid.geoLat_scgc,
                           grid.geoAlt_scgc);

        // This is just to check if MSIS is actually working:
        if (msis.is_valid_species("Tn"))
          // if it is, fill will temperature:
          temperature_scgc = msis.get_cube("Tn");
        else
          // if it is not, then fill with a value:
          temperature_scgc.fill(initial_temperatures[0]);

        for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
          if (report.test_verbose(3))
            std::cout << "Setting Species : "
                      << species[iSpecies].cName << "\n";

          if (msis.is_valid_species(species[iSpecies].cName)) {
            if (report.test_verbose(3))
              std::cout << "  Found in MSIS!\n";

            species[iSpecies].density_scgc =
              msis.get_cube(species[iSpecies].cName);
          } else {
            if (report.test_verbose(3))
              std::cout << "  NOT Found in MSIS - setting constant\n";

            species[iSpecies].density_scgc.slice(0).
            fill(species[iSpecies].lower_bc_density);
            fill_with_hydrostatic(iSpecies, 1, nAlts, grid);
          }

        } // for species
      } // msis init worked ok
    } // type = msis

    if (icsType == "planet") {

      didWork = true;

      // ---------------------------------------------------------------------
      // This section assumes we want a hydrostatic solution given the
      // temperature profile in the planet.in file.
      // ---------------------------------------------------------------------

      int64_t nLons = grid.get_nLons();
      int64_t nLats = grid.get_nLats();
      int64_t nAlts = grid.get_nAlts();

      // Let's assume that the altitudes are not dependent on lat/lon:

      arma_vec alt1d(nAlts);
      arma_vec temp1d(nAlts);

      arma_mat H2d(nLons, nLats);

      alt1d = grid.geoAlt_scgc.tube(0, 0);

      if (nInitial_temps > 0) {
        for (iAlt = 0; iAlt < nAlts; iAlt++) {
          alt = alt1d(iAlt);

          // Find temperatures:
          if (alt <= initial_altitudes[0])
            temp1d[iAlt] = initial_temperatures[0];

          else {
            if (alt >= initial_altitudes[nInitial_temps - 1])
              temp1d[iAlt] = initial_temperatures[nInitial_temps - 1];

            else {
              // Linear interpolation!
              iA = 0;

              while (alt > initial_altitudes[iA])
                iA++;

              iA--;
              // alt will be between iA and iA+1:
              r = (alt - initial_altitudes[iA]) /
                  (initial_altitudes[iA + 1] - initial_altitudes[iA]);
              temp1d[iAlt] =
                (1.0 - r) * initial_temperatures[iA] +
                (r) * initial_temperatures[iA + 1];
            }
          }
        }
      } else
        temp1d = 200.0;

      // spread the 1D temperature across the globe:
      for (iLon = 0; iLon < nLons; iLon++) {
        for (iLat = 0; iLat < nLats; iLat++)
          temperature_scgc.tube(iLon, iLat) = temp1d;
      }

      // Set the lower boundary condition:
      for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
        species[iSpecies].density_scgc.slice(0).
        fill(species[iSpecies].lower_bc_density);
      }

      calc_scale_height(grid);
      fill_with_hydrostatic(1, nAlts, grid);
    } // type = planet
  }

  if (!didWork)
    report.error("Issue with initial conditions!");

  report.exit(function);
  return didWork;
}



