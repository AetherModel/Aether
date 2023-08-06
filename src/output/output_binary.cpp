// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "aether.h"

#include <fstream>

//----------------------------------------------------------------------
// Figure out which variables to the binary file.
// ----------------------------------------------------------------------

int write_binary_all_3d(std::string file_name,
                        std::string type_output,
                        Neutrals neutrals,
                        Ions ions,
                        Grid grid,
                        Times time,
                        Planets planet) {

  int iErr = 0;
  std::ofstream binary;
  binary.open(file_name, ios::binary | ios::out);

  // Everything gets (geo) lon/lat/alt:
  output_variable_3d(binary, grid.geoLon_scgc);
  output_variable_3d(binary, grid.geoLat_scgc);
  output_variable_3d(binary, grid.geoAlt_scgc);

  // Neutral States
  if (type_output == "neutrals" ||
      type_output == "states") {
    for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      output_variable_3d(binary, neutrals.species[iSpecies].density_scgc);

    output_variable_3d(binary, neutrals.temperature_scgc);
  }

  // Ion States:
  if (type_output == "ions" ||
      type_output == "states") {
    for (int iSpecies = 0; iSpecies < nIons + 1; iSpecies++)
      output_variable_3d(binary, ions.species[iSpecies].density_scgc);

    for (int iComp = 0; iComp < 3; iComp++)
      output_variable_3d(binary, ions.velocity_vcgc[iComp]);

    output_variable_3d(binary, ions.potential_scgc);
  }

  // Magnetic Field Variables:
  if (type_output == "bfield") {
    output_variable_3d(binary, grid.magLat_scgc);
    output_variable_3d(binary, grid.magLon_scgc);
    output_variable_3d(binary, grid.magLocalTime_scgc);

    for (int iComp = 0; iComp < 3; iComp++)
      output_variable_3d(binary, grid.bfield_vcgc[iComp]);
  }

  binary.close();
  return iErr;

}

//----------------------------------------------------------------------
// Write header files. Supported types:
// 1. neutrals: Lon/Lat/Alt + nSpecies + temp
// 2. ions: Lon/Lat/Alt + nIons + [e-] + Vi + potential
// 3. states: Lon/Lat/Alt + nSpecies + temp + nIons + [e-] + Vi + potential
//----------------------------------------------------------------------

int write_header(std::string file_name,
                 std::string type_output,
                 Neutrals neutrals,
                 Ions ions,
                 Grid grid,
                 Times time,
                 Planets planet) {

  int iErr = 0;
  std::ofstream header;
  header.open(file_name);

  // Everything gets lon/lat/alt:
  int nVars = 3;

  if (type_output == "neutrals" ||
      type_output == "states")
    // All neutrals, temperature
    nVars = nVars + nSpecies + 1;

  if (type_output == "ions" ||
      type_output == "states")
    // All ions, electrons, Vi, potential:
    nVars = nVars + nIons + 1 + 3 + 1;

  if (type_output == "bfield")
    nVars = nVars + 6;

  header << "\n";
  header << "BLOCKS\n";
  header << "1\n";
  header << "1\n";
  header << "1\n";

  header << "\n";
  header << "TIME\n";
  std::vector<int> iCurrent = time.get_iCurrent();

  for (int i = 0; i < 7; i++)
    header << iCurrent[i] << "\n";

  header << "\n";
  header << "VERSION\n";
  header << "0.10\n";

  header << "\n";
  header << "NUMERICAL VALUES\n";
  header << nVars << " nVars\n";
  header << grid.get_nX() << " nLons\n";
  header << grid.get_nY() << " nLats\n";
  header << grid.get_nZ() << " nAlts\n";

  header << "\n";
  header << "VARIABLE LIST\n";
  header << "1 Longitude (radians)\n";
  header << "2 Latitude (radians)\n";
  header << "3 Altitude (m)\n";

  int iVar = 4;

  if (type_output == "neutrals" ||
      type_output == "states") {

    for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      header << iVar << " "
             << neutrals.species[iSpecies].cName << " "
             << neutrals.density_unit << "\n";
      iVar++;
    }

    header << iVar << " "
           << neutrals.temperature_name << " "
           << neutrals.temperature_unit << "\n";
    iVar++;
  }

  if (type_output == "ions" ||
      type_output == "states") {

    for (int iSpecies = 0; iSpecies < nIons + 1; iSpecies++) {
      header << iVar << " "
             << ions.species[iSpecies].cName << " "
             << neutrals.density_unit << "\n";
      iVar++;
    }

    header << iVar << " Ion Velocity (East) (m/s)\n";
    header << iVar << " Ion Velocity (North) (m/s)\n";
    header << iVar << " Ion Velocity (Vertical) (m/s)\n";
    header << iVar << " Potential (V)\n";
    iVar++;

  }

  if (type_output == "bfield") {
    header << iVar << " Magnetic Latitude (radians)\n";
    iVar++;
    header << iVar << " Magnetic Longitude (radians)\n";
    iVar++;
    header << iVar << " Magnetic Local Time (hours)\n";
    iVar++;
    header << iVar << " Beast (nT)\n";
    iVar++;
    header << iVar << " Bnorth (nT)\n";
    iVar++;
    header << iVar << " Bvertical (nT)\n";
    iVar++;
  }

  header << "\n";
  header << "END\n";
  header << "\n";

  header.close();
  return iErr;
}

//----------------------------------------------------------------------
// Output the different file types to binary files.
//----------------------------------------------------------------------

int output(Neutrals neutrals,
           Ions ions,
           Grid grid,
           Times time,
           Planets planet) {

  int iErr = 0;

  int nOutputs = input.get_n_outputs();

  int64_t nLons = grid.get_nLons();
  int64_t nLats = grid.get_nLats();
  int64_t nAlts = grid.get_nAlts();
  double time_array[1];

  std::string function = "output";
  static int iFunction = -1;
  report.enter(function, iFunction);

  for (int iOutput = 0; iOutput < nOutputs; iOutput++) {

    if (time.check_time_gate(input.get_dt_output(iOutput))) {

      grid.calc_sza(planet, time);
      grid.calc_gse(planet, time);
      grid.calc_mlt();

      std::string time_string;
      std::string file_name;
      std::string file_pre;

      std::string type_output = input.get_type_output(iOutput);
      std::string output_dir = "UA/output/";

      if (type_output == "neutrals")
        file_pre = "3DNEU";

      if (type_output == "states")
        file_pre = "3DALL";

      if (type_output == "bfield")
        file_pre = "3DBFI";

      time_string = time.get_YMD_HMS();
      std::string file_ext = ".header";
      file_name = output_dir + "/" + file_pre + "_" + time_string + file_ext;

      // Create the file:
      report.print(0, "Writing file : " + file_name);

      iErr = write_header(file_name,
                          type_output,
                          neutrals,
                          ions,
                          grid,
                          time,
                          planet);

      file_ext = ".bin";
      file_name = output_dir + "/" + file_pre + "_" + time_string + file_ext;

      // Create the file:
      report.print(0, "Writing file : " + file_name);

      iErr = write_binary_all_3d(file_name,
                                 type_output,
                                 neutrals,
                                 ions,
                                 grid,
                                 time,
                                 planet);

    }  // if time check
  }  // for iOutput

  report.exit(function);
  return iErr;
}
