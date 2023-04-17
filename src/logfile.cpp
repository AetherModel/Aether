#include "aether.h"
#include "logfile.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

Logfile::Logfile(Indices indices,
		 Inputs inputs,
		 Report &report){
  // should add a switch that alters this between trunc and app
  // restart should be app, while starting from scratch should be trunc
  logfilestream.open(inputs.get_logfile(), ofstream::trunc);
  cout << inputs.get_logfile();
  logfilestream << "year month day hour minute second ";

  if (!header){
    vector<string> species_array = inputs.get_species_vector();
    for (int i = 0; i < indices.all_indices_array_size(); i++){
      if (indices.get_nValues(i)>0){
	logfilestream << indices.get_name(i) << " ";
      }
    } 
    logfilestream << "min_temp max_temp mean_temp ";
    string name;
    for (int i = 0; i < species_array.size(); i++){
      name = species_array.at(i);
      logfilestream << name << "_min " << name << "_max " << name << "_mean ";
    }
    logfilestream << "\n";
    header = true;
  }

  
}

//-------------------------------------------------------------
// Get min, mean, and max of an arma_cube
//-------------------------------------------------------------

vector<precision_t> get_min_mean_max(arma_cube value) {
  vector<precision_t> mmm(3);
  mmm[0] = value.min();
  mmm[1] = accu(value)/value.n_elem;
  mmm[2] = value.max();
  return mmm;
}

//-------------------------------------------------------------
// Get min, mean, and max of either a neutral or ion species
//-------------------------------------------------------------

vector<precision_t> get_min_mean_max_density(std::string name,
					     Neutrals neutrals,
					     Ions ions,
					     Report report) {
  vector<precision_t> mmm(3);
  
  int id = neutrals.get_species_id(name, report);
  if (id > -1){
    mmm = get_min_mean_max(neutrals.species[id].density_scgc);
  } else {
    id = ions.get_species_id(name, report);
    if (id > -1) {
      mmm = get_min_mean_max(ions.species[id].density_scgc);
    } else {
      mmm[0] = 0.0;
      mmm[1] = 0.0;
      mmm[2] = 0.0;
    }
  }
  return mmm;
}

//-------------------------------------------------------------
// Add a new line to the log file
//-------------------------------------------------------------

void Logfile::write_logfile(Times time,
			    Neutrals neutrals,
			    Ions ions,
			    Inputs inputs,
			    Indices indices,
			    Report report) {

  // output time first:
  vector<int> itime = time.get_iCurrent();
  vector<string> species_array = inputs.get_species_vector();
  for (int i = 0; i <= 6; ++i){
    logfilestream << std::setw(2) << itime.at(i) << " ";
  }

  // output indices next:
  logfilestream.precision(4);
  for (int i = 0; i < indices.all_indices_array_size(); i++){
    if (indices.get_nValues(i) > 0)
      logfilestream << indices.get_index(time.get_current(), i) << " ";
  }

  // Output Neutral Temperatures:
  vector<precision_t> temp_stats = get_min_mean_max(neutrals.temperature_scgc);
  for (int i = 0; i < 3; i++)
    logfilestream << temp_stats[i] << " ";

  // Output densities of requested species:
  for(int i = 0; i<species_array.size(); i++){
    vector<precision_t> density_stats =
      get_min_mean_max_density(species_array[i], neutrals, ions, report);
    for (int i = 0; i < 3; i++)
      logfilestream << density_stats[i] << " ";
  }
  logfilestream << "\n";
}

void Logfile::close_logfile(){
    logfilestream.close();
}

