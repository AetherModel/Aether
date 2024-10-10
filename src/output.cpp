// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "aether.h"

/* ---------------------------------------------------------------------

   Fill output containers for certain output types.  Supported types:
   states - neutral states, ion den., bulk ion vel., temp, elec. temp
   neutrals - neutral states
   ions - Ion densites, temperatures, par & perp velocities, elec temp
   bfield - magnetic coordinates, b-field vector

 -------------------------------------------------------------------- */

std::string get_filename_from_type(std::string type_output) {

  std::string filename = "";

  if (type_output == "neutrals")
    filename = "3DNE";

  if (type_output == "states")
    filename = "3DAL";

  if (type_output == "ions")
    filename = "3DIO";

  if (type_output == "bfield")
    filename = "3DBF";

  if (type_output == "moment")
    filename = "3DMO";

  if (type_output == "gravity")
    filename = "3DGR";

  if (type_output == "corners" || type_output == "grid")
    filename = "3DCO";

  if (type_output == "therm")
    filename = "3DTH";

  return filename;

} 

// -----------------------------------------------------------------------------
//  Fills output containers and outputs them for common output types
// -----------------------------------------------------------------------------

bool output(const Neutrals &neutrals,
            const Ions &ions,
            Grid &grid,
            Times time,
            const Planets &planet) {

  std::string function = "output";
  static int iFunction = -1;
  report.enter(function, iFunction);

  static bool IsFirstTime = true;

  bool didWork = true;

  int nOutputs = input.get_n_outputs();
  static std::vector<OutputContainer> AllOutputContainers;

  if (IsFirstTime) {
    // Initialize all of the output containers for all of the output
    // types requested
    OutputContainer DummyOutputContainer;
    std::string output_dir = "UA/output/";
    DummyOutputContainer.set_directory(output_dir);
    DummyOutputContainer.set_version(aether_version);
    DummyOutputContainer.set_nGhostCells(grid.get_nGCs());

    for (int iOutput = 0; iOutput < nOutputs; iOutput++)
      AllOutputContainers.push_back(DummyOutputContainer);

    IsFirstTime = false;
  }

  report.student_checker_function_name(input.get_is_student(),
                                       input.get_student_name(),
                                       3, "");

  for (int iOutput = 0; iOutput < nOutputs; iOutput++) {

    if (time.check_time_gate(input.get_dt_output(iOutput))) {

      // ------------------------------------------------------------
      // Store time in all of the files:

      AllOutputContainers[iOutput].set_time(time.get_current());

      std::string type_output = input.get_type_output(iOutput);

      // ------------------------------------------------------------
      // Put Lon, Lat, Alt into all output containers:

      if (type_output == "corners" || type_output == "grid") {
        // Cell Corners:
        AllOutputContainers[iOutput].
        store_variable("lon",
                       "longitude",
                       "degrees_east",
                       grid.geoLon_Corner * cRtoD);
        AllOutputContainers[iOutput].
        store_variable("lat",
                       "latitude",
                       "degrees_north",
                       grid.geoLat_Corner * cRtoD);
        AllOutputContainers[iOutput].
        store_variable("z",
                       "height above mean sea level",
                       "m",
                       grid.geoAlt_Corner);
      } else {
        // Cell Centers:
        AllOutputContainers[iOutput].
        store_variable("lon",
                       "longitude",
                       "degrees_east",
                       grid.geoLon_scgc * cRtoD);
        AllOutputContainers[iOutput].
        store_variable("lat",
                       "latitude",
                       "degrees_north",
                       grid.geoLat_scgc * cRtoD);
        AllOutputContainers[iOutput].
        store_variable("z",
                       "height above mean sea level",
                       "m",
                       grid.geoAlt_scgc);
      }

      // ------------------------------------------------------------
      // Put certain variables into each file type

      // Neutral Densities:
      if (type_output == "neutrals" ||
          type_output == "states")
        for (int iSpecies = 0; iSpecies < neutrals.nSpecies; iSpecies++)
          AllOutputContainers[iOutput].
          store_variable("density_" + neutrals.species[iSpecies].cName,
                         neutrals.density_unit,
                         neutrals.species[iSpecies].density_scgc);

      // Neutral Temperature:
      if (type_output == "neutrals" ||
          type_output == "states")
        AllOutputContainers[iOutput].
        store_variable(neutrals.temperature_name + "_neutral",
                       neutrals.temperature_unit,
                       neutrals.temperature_scgc);

      // Bulk Neutral Winds:
      if (type_output == "neutrals" ||
          type_output == "states")
        for (int iDir = 0; iDir < 3; iDir++)
          AllOutputContainers[iOutput].
          store_variable(neutrals.velocity_name[iDir] + "_neutral",
                         neutrals.velocity_unit,
                         neutrals.velocity_vcgc[iDir]);

      // Neutral Species Winds:
      if (type_output == "neutrals")
        for (int iSpecies = 0; iSpecies < neutrals.nSpecies; iSpecies++)
          for (int iDir = 0; iDir < 3; iDir++)
            AllOutputContainers[iOutput].
            store_variable(neutrals.velocity_name[iDir] + "_" +
                           neutrals.species[iSpecies].cName,
                           neutrals.velocity_unit,
                           neutrals.species[iSpecies].velocity_vcgc[iDir]);

      // Ion Densities:
      if (type_output == "ions" ||
          type_output == "states")
        for (int iSpecies = 0; iSpecies <= ions.nSpecies; iSpecies++)
          AllOutputContainers[iOutput].
          store_variable("density_" + ions.species[iSpecies].cName,
                         ions.density_unit,
                         ions.species[iSpecies].density_scgc);

      // Ion Temperatures:
      if (type_output == "ions" ||
          type_output == "states")
        for (int iSpecies = 0; iSpecies <= ions.nSpecies; iSpecies++)
          AllOutputContainers[iOutput].
          store_variable(ions.temperature_name + "_" +
                         ions.species[iSpecies].cName,
                         ions.temperature_unit,
                         ions.species[iSpecies].temperature_scgc);

      // Bulk Ion Temperature:
      if (type_output == "ions" ||
          type_output == "states")
        AllOutputContainers[iOutput].store_variable(ions.temperature_name + "_ion",
                                                    ions.temperature_unit,
                                                    ions.temperature_scgc);

      // Bulk Ion Drifts:
      if (type_output == "states")
        for (int iDir = 0; iDir < 3; iDir++)
          AllOutputContainers[iOutput].store_variable(ions.velocity_name[iDir] + "_ion",
                                                      ions.velocity_unit,
                                                      ions.velocity_vcgc[iDir]);

      // Electric Potential:
      if (type_output == "ions" ||
          type_output == "states")
        AllOutputContainers[iOutput].store_variable(ions.potential_name,
                                                    ions.potential_unit,
                                                    ions.potential_scgc);

      if (type_output == "gravity") {
        AllOutputContainers[iOutput].store_variable("glat",
                                                    "grav Latitude",
                                                    "degrees",
                                                    grid.geoLat_scgc * cRtoD);
        AllOutputContainers[iOutput].store_variable("glon",
                                                    "grav Longitude",
                                                    "degrees",
                                                    grid.geoLon_scgc * cRtoD);
        AllOutputContainers[iOutput].store_variable("lt",
                                                    "Local Time",
                                                    "hours",
                                                    grid.geoLocalTime_scgc);
        AllOutputContainers[iOutput].store_variable("Geast",
                                                    "m/s^2",
                                                    grid.gravity_vcgc[0]);
        AllOutputContainers[iOutput].store_variable("Gnorth",
                                                    "m/s^2",
                                                    grid.gravity_vcgc[1]);
        AllOutputContainers[iOutput].store_variable("Gvertical",
                                                    "m/s^2",
                                                    grid.gravity_vcgc[2]);
        AllOutputContainers[iOutput].store_variable("Gpotential",
                                                    "m^2/s^2",
                                                    grid.gravity_potential_scgc);
        AllOutputContainers[iOutput].store_variable("radius",
                                                    "m",
                                                    grid.radius_scgc);
      }

      if (type_output == "bfield") {
        AllOutputContainers[iOutput].store_variable("mlat",
                                                    "Magnetic Latitude",
                                                    "degrees",
                                                    grid.magLat_scgc * cRtoD);
        AllOutputContainers[iOutput].store_variable("mlon",
                                                    "Magnetic Longitude",
                                                    "degrees",
                                                    grid.magLon_scgc * cRtoD);
        AllOutputContainers[iOutput].store_variable("mlt",
                                                    "Magnetic Local Time",
                                                    "hours",
                                                    grid.magLocalTime_scgc);
        AllOutputContainers[iOutput].store_variable("Beast",
                                                    "nT",
                                                    grid.bfield_vcgc[0]);
        AllOutputContainers[iOutput].store_variable("Bnorth",
                                                    "nT",
                                                    grid.bfield_vcgc[1]);
        AllOutputContainers[iOutput].store_variable("Bvertical",
                                                    "nT",
                                                    grid.bfield_vcgc[2]);
      }

      // Thermal:
      if (type_output == "therm") {
        AllOutputContainers[iOutput].store_variable("Heating_EUV",
                                                    "Heating from EUV",
                                                    "K/s",
                                                    neutrals.heating_euv_scgc);
        AllOutputContainers[iOutput].store_variable("Heating_Chemistry",
                                                    "Heating from Chemistry",
                                                    "K/s",
                                                    neutrals.heating_chemical_scgc);
        AllOutputContainers[iOutput].store_variable("Heating_Transfer",
                                                    "Heating from Ti- Tn Ion Neutral Collisions",
                                                    "K/s",
                                                    neutrals.heating_ion_heat_transfer_scgc);
        AllOutputContainers[iOutput].store_variable("Heating_Ion_Friction",
                                                    "Heating from Friction Ion Neutral Collisions",
                                                    "K/s",
                                                    neutrals.heating_ion_friction_scgc);
        AllOutputContainers[iOutput].store_variable("Conduction",
                                                    "Conduction",
                                                    "K/s",
                                                    neutrals.conduction_scgc);
        AllOutputContainers[iOutput].store_variable("O Rad Cooling",
                                                    "[O] Radiative Cooling",
                                                    "K/s",
                                                    neutrals.O_cool_scgc);
        AllOutputContainers[iOutput].store_variable("NO Rad Cooling",
                                                    "[NO] Radiative Cooling",
                                                    "K/s",
                                                    neutrals.NO_cool_scgc);
      }

      if (type_output == "moment") {
        AllOutputContainers[iOutput].store_variable("accel_cent_east",
                                                    "Logitudinal Centripetal Acceleration",
                                                    "m/s^2",
                                                    grid.cent_acc_vcgc[0]);
        AllOutputContainers[iOutput].store_variable("accel_cent_north",
                                                    "Latitudinal Centripetal Acceleration",
                                                    "m/s^2",
                                                    grid.cent_acc_vcgc[1]);
        AllOutputContainers[iOutput].store_variable("accel_cent_up",
                                                    "Radial Centripetal Acceleration",
                                                    "m/s^2",
                                                    grid.cent_acc_vcgc[2]);
      }

      // ------------------------------------------------------------
      // Set output file names

      std::string filename = get_filename_from_type(type_output);

      if (filename.length() < 4) {
        report.print(0, "type_output : " + type_output);
        report.error("File output type not found!");
        didWork = false;
      } else {
        if (grid.get_IsGeoGrid())
          filename = filename + "G_";
        else
          filename = filename + "M_";

        if ((int64_t(input.get_dt_output(iOutput)) % 60) == 0)
          filename = filename + time.get_YMD_HM0();
        else
          filename = filename + time.get_YMD_HMS();

        if (nMembers > 1)
          filename = filename + "_" + cMember;

        filename = filename + "_" + cGrid;

        report.print(0, "Writing file : " + filename);
        AllOutputContainers[iOutput].set_filename(filename);

        // ------------------------------------------------------------
        // write output container

        if (!AllOutputContainers[iOutput].write()) {
          report.error("Error in writing output container!");
          didWork = false;
        }
      }

      // ------------------------------------------------------------
      // Clear variables for next time

      AllOutputContainers[iOutput].clear_variables();

    }

    if (!didWork)
      break;
  }

  report.exit(function);
  return didWork;
}

