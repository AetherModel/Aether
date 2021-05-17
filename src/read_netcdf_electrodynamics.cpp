
#include <iostream>
#include <string>
#include <stdlib.h>
#include <netcdf>
#include "read_amie_bin.cpp"
#include "time_conversion.cpp"
#include "../include/electrodynamics.h"
#include <vector>
using namespace netCDF;
//g++ read_netcdf_electrodynamics.cpp -o -std=c++11

//done with binary currently, not netcdf
void Electrodynamics::read_netcdf_electrodynamics_file(std::string filename, Report &report){
    int reclen;
    int nLats, nMlts, nTimes, nVars;
    std::vector<double> real_times;
    fvec mlats_struct, mlts_struct;
    std::vector<fmat> potential_struct, energy_flux_struct, average_energy_struct, 
    ion_energy_flux_struct, ion_average_energy_struct;

    std::cout << filename << "\n";

    FILE *infile;

    char* char_arr;
    char_arr = &filename[0];
    infile = fopen(char_arr, "rb");
    if (infile == NULL) {
        std::cout << "can't find the file!\n";
    } else {

        // read in number of lats, mlts, times:
        fread(&reclen, sizeof(reclen), 1, infile);
        fread(&nLats, sizeof(nLats), 1, infile);
        fread(&nMlts, sizeof(nMlts), 1, infile);
        fread(&nTimes, sizeof(nTimes), 1, infile);
        std::cout << nLats << " "
            << nMlts << " "
            << nTimes << "\n";
        fread(&reclen, sizeof(reclen), 1, infile);

        std::vector<float> lats, mlts;
        lats = read_in_float_array(infile, nLats);
        mlts = read_in_float_array(infile, nMlts);
        mlats_struct = conv_to<fmat>::from(lats);
        mlts_struct = conv_to<fmat>::from(mlts);


        // read in number of variables
        fread(&reclen, sizeof(reclen), 1, infile);
        fread(&nVars, sizeof(nVars), 1, infile);
        std::cout << nVars << "\n";
        fread(&reclen, sizeof(reclen), 1, infile);

        std::vector<std::string> Vars;
        for (int i = 0; i < nVars; i++) {
        Vars.push_back(read_in_string(infile));
        }
        std::vector<int> itime;
        std::vector<float> indices;

        fmat values;
        std::vector<fmat> values_one_time;
        std::vector<std::vector<fmat>> all_values;

        for (int it = 0; it < nTimes; it++) {

            itime = read_in_int_array(infile, 6);
            itime.erase(itime.begin());
            itime.push_back(0);
            itime.push_back(0);
            double real_time = time_int_to_real(itime);
            real_times.push_back(real_time);

            indices = read_in_float_array(infile, 13);
            
            for (int iv = 0; iv < nVars; iv++) { 
                values = read_in_fmat_array(infile, nMlts, nLats);
                if (it == 0) {
                values_one_time.push_back(values);
                } else {
                values_one_time[iv] = values;
                }
            }
            potential_struct.push_back(values_one_time[0]);
            energy_flux_struct.push_back(values_one_time[1]);
            average_energy_struct.push_back(values_one_time[2]);
            if (nVars == 5){
              ion_energy_flux_struct.push_back(values_one_time[3]);
              ion_average_energy_struct.push_back(values_one_time[4]);
            }
            all_values.push_back(values_one_time);
          }      
        }
    fclose(infile);
   
    input_electrodynamics_struct obj = input_electrodynamics_struct();
    obj.nLats = nLats;
    obj.nMlts = nMlts;
    obj.times = real_times;
    obj.mlats = mlats_struct;
    obj.mlts = mlts_struct;

    obj.potential = potential_struct;
    obj.energy_flux = energy_flux_struct;
    obj.average_energy = average_energy_struct;
    obj.ion_energy_flux = ion_energy_flux_struct;
    obj.ion_average_energy = ion_average_energy_struct;

    if (nVars == 5){
      obj.DoesIncludeIonPrecip = 1;
    }

    input_electrodynamics.push_back(obj);
}
