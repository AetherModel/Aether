// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "../include/aether.h"

#include <fstream>

//----------------------------------------------------------------------
// Output a given variable to the binary file. 
// ----------------------------------------------------------------------


void output_variable_3d(std::ofstream &binary,
                        fcube value) {

  // Get the size of the cube:
  
  int64_t nX = value.n_rows;
  int64_t nY = value.n_cols;
  int64_t nZ = value.n_slices;
  int64_t iX, iY, iZ, index;

  int64_t nPts = nX * nY * nZ;
  int64_t iTotalSize = nPts * sizeof(float);

  // Create a temporary c-array to use to output the variable
  float *tmp_s3gc = static_cast<float*>(malloc(iTotalSize));

  // Move the data from the cube to the c-array
  
  for (iX = 0; iX < nX; iX++) {
    for (iY = 0; iY < nY; iY++) {
      for (iZ = 0; iZ < nZ; iZ++) {
        index = iX*nY*nZ + iY*nZ + iZ;
        tmp_s3gc[index] = value(iX, iY, iZ);
      }
    }
  }

  // Output the data to the binary file
  binary.write((char *) tmp_s3gc, iTotalSize);  

  // delete the c-array
  free(tmp_s3gc);
}

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
    for (int iSpecies=0; iSpecies < nSpecies; iSpecies++)
      output_variable_3d(binary, neutrals.species[iSpecies].density_scgc);
    output_variable_3d(binary, neutrals.temperature_scgc);
  }

  // Ion States:
  if (type_output == "ions" ||
      type_output == "states") {
    for (int iSpecies=0; iSpecies < nIons+1; iSpecies++) 
      output_variable_3d(binary, ions.species[iSpecies].density_scgc);
  }

  // Magnetic Field Variables:
  if (type_output == "bfield") {
    output_variable_3d(binary, grid.magLat_scgc);
    output_variable_3d(binary, grid.magLon_scgc);
    output_variable_3d(binary, grid.magLocalTime_scgc);
    for (int iComp=0; iComp < 3; iComp++) 
      output_variable_3d(binary, grid.bfield_vcgc[iComp]);
  }

  binary.close();
  return iErr;
  
}

//----------------------------------------------------------------------
// Write header files. Supported types:
// 1. neutrals: Lon/Lat/Alt + nSpecies + temp
// 2. ions: Lon/Lat/Alt + nIons + [e-]
// 3. states: Lon/Lat/Alt + nSpecies + temp + nIons + [e-]
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
    nVars = nVars + nSpecies + 1;

  if (type_output == "ions" ||
      type_output == "states")
    nVars = nVars + nIons + 1;
  
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
  for (int i = 0; i < 7; i++) header << iCurrent[i] << "\n";

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
    
    for (int iSpecies=0; iSpecies < nSpecies; iSpecies++) {
      header << iVar << " "
	     << neutrals.species[iSpecies].cName << " "
	     << "(" << neutrals.density_unit << ")\n";
      iVar++;
    }

    header << iVar << " "
	   << neutrals.temperature_name << " "
	   << neutrals.temperature_unit << "\n";
    iVar++;
  }    

  if (type_output == "ions" ||
      type_output == "states") {
    
    for (int iSpecies=0; iSpecies < nIons+1; iSpecies++) {
      header << iVar << " "
	     << ions.species[iSpecies].cName << " "
	     << neutrals.density_unit << "\n";
      iVar++;
    }

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
           Planets planet,
           Inputs args,
           Report &report) {

  int iErr = 0;

  int nOutputs = args.get_n_outputs();

  int64_t nLons = grid.get_nLons();
  int64_t nLats = grid.get_nLats();
  int64_t nAlts = grid.get_nAlts();
  double time_array[1];

  std::string function = "output";
  static int iFunction = -1;
  report.enter(function, iFunction);

  for (int iOutput = 0; iOutput < nOutputs; iOutput++) {

    if (time.check_time_gate(args.get_dt_output(iOutput))) {

      grid.calc_sza(planet, time, report);
      grid.calc_gse(planet, time, report);
      grid.calc_mlt(report);

      std::string time_string;
      std::string file_name;
      std::string file_pre;

      std::string type_output = args.get_type_output(iOutput);
      std::string output_dir = args.get_output_directory();

      if (type_output == "neutrals") file_pre = "3DNEU";
      if (type_output == "states") file_pre = "3DALL";
      if (type_output == "bfield") file_pre = "3DBFI";

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
      
      
      /*
      NcFile ncdf_file(file_name, NcFile::replace);

      // Add dimensions:
      NcDim lonDim = ncdf_file.addDim("Longitude", nLons);
      NcDim latDim = ncdf_file.addDim("Latitude", nLats);
      NcDim altDim = ncdf_file.addDim("Altitude", nAlts);

      NcDim timeDim = ncdf_file.addDim("Time", 1);

      // Define the Coordinate Variables

      // Define the netCDF variables for the 3D data.
      // First create a vector of dimensions:

      std::vector<NcDim> dimVector;
      dimVector.push_back(lonDim);
      dimVector.push_back(latDim);
      dimVector.push_back(altDim);

      NcVar timeVar = ncdf_file.addVar("Time", ncDouble, timeDim);
      NcVar lonVar = ncdf_file.addVar("Longitude", ncFloat, dimVector);
      NcVar latVar = ncdf_file.addVar("Latitude", ncFloat, dimVector);
      NcVar altVar = ncdf_file.addVar("Altitude", ncFloat, dimVector);

      timeVar.putAtt(UNITS, "seconds");
      lonVar.putAtt(UNITS, "radians");
      latVar.putAtt(UNITS, "radians");
      altVar.putAtt(UNITS, "meters");

      std::vector<size_t> startp, countp;
      startp.push_back(0);
      startp.push_back(0);
      startp.push_back(0);

      countp.push_back(nLons);
      countp.push_back(nLats);
      countp.push_back(nAlts);

      // Output time:

      time_array[0] = time.get_current();
      timeVar.putVar(time_array);

      // Output longitude, latitude, altitude 3D arrays:

      output_variable_3d(startp, countp, grid.geoLon_scgc, lonVar);
      output_variable_3d(startp, countp, grid.geoLat_scgc, latVar);
      output_variable_3d(startp, countp, grid.geoAlt_scgc, altVar);

      // ----------------------------------------------
      // Neutral Densities and Temperature
      // ----------------------------------------------

      if (type_output == "neutrals" ||
          type_output == "states") {

        // Output all species densities:
        std::vector<NcVar> denVar;
        for (int iSpecies=0; iSpecies < nSpecies; iSpecies++) {
          if (report.test_verbose(3))
            std::cout << "Outputting Var : "
                      << neutrals.species[iSpecies].cName << "\n";
          denVar.push_back(ncdf_file.addVar(neutrals.species[iSpecies].cName,
                                            ncFloat, dimVector));
          denVar[iSpecies].putAtt(UNITS, neutrals.density_unit);
          output_variable_3d(startp, countp,
                             neutrals.species[iSpecies].density_scgc,
                             denVar[iSpecies]);
        }

        // Output bulk temperature:
        NcVar tempVar = ncdf_file.addVar(neutrals.temperature_name,
                                         ncFloat, dimVector);
        tempVar.putAtt(UNITS, neutrals.temperature_unit);
        output_variable_3d(startp, countp, neutrals.temperature_scgc, tempVar);

        // Output SZA
        NcVar szaVar = ncdf_file.addVar("Solar Zenith Angle",
          ncFloat, dimVector);
        szaVar.putAtt(UNITS, "degrees");
        output_variable_3d(startp, countp, grid.sza_scgc * cRtoD, szaVar);

      }

      // ----------------------------------------------
      // Ion Densities and Ion Temperature and Electron Temperature
      // ----------------------------------------------

      if (type_output == "ions" ||
          type_output == "states") {

        // Output all species densities:
        std::vector<NcVar> ionVar;
        for (int iSpecies=0; iSpecies < nIons; iSpecies++) {
          if (report.test_verbose(3))
            std::cout << "Outputting Var : "
                      << ions.species[iSpecies].cName << "\n";
          ionVar.push_back(ncdf_file.addVar(ions.species[iSpecies].cName,
                                            ncFloat, dimVector));
          ionVar[iSpecies].putAtt(UNITS, neutrals.density_unit);
          output_variable_3d(startp, countp,
                             ions.species[iSpecies].density_scgc,
                             ionVar[iSpecies]);
        }

        ionVar.push_back(ncdf_file.addVar("e-", ncFloat, dimVector));
        ionVar[nIons].putAtt(UNITS, neutrals.density_unit);
        output_variable_3d(startp, countp, ions.density_scgc, ionVar[nIons]);
      }

      // ----------------------------------------------
      // Magnetic field
      // ----------------------------------------------

      if (type_output == "bfield") {
        NcVar mLatVar = ncdf_file.addVar("Magnetic Latitude",
                                         ncFloat, dimVector);
        mLatVar.putAtt(UNITS, "radians");
        output_variable_3d(startp, countp, grid.magLat_scgc, mLatVar);

        NcVar mLonVar = ncdf_file.addVar("Magnetic Longitude",
                                         ncFloat, dimVector);
        mLonVar.putAtt(UNITS, "radians");
        output_variable_3d(startp, countp, grid.magLat_scgc, mLonVar);

        NcVar mLTVar = ncdf_file.addVar("Magnetic Local Time",
                                         ncFloat, dimVector);
        mLTVar.putAtt(UNITS, "hours");
        output_variable_3d(startp, countp, grid.magLocalTime_scgc, mLTVar);

        // Output magnetic field components:

        NcVar bxVar = ncdf_file.addVar("Bx", ncFloat, dimVector);
        bxVar.putAtt(UNITS, "nT");
        output_variable_3d(startp, countp, grid.bfield_vcgc[0], bxVar);

        NcVar byVar = ncdf_file.addVar("By", ncFloat, dimVector);
        byVar.putAtt(UNITS, "nT");
        output_variable_3d(startp, countp, grid.bfield_vcgc[1], bxVar);

        NcVar bzVar = ncdf_file.addVar("Bz", ncFloat, dimVector);
        bzVar.putAtt(UNITS, "nT");
        output_variable_3d(startp, countp, grid.bfield_vcgc[2], bxVar);
      }  // if befield

      ncdf_file.close();

      */
    }  // if time check
  }  // for iOutput

  report.exit(function);
  return iErr;
}
