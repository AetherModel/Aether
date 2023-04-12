#include "aether.h"
#include "logfile.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

Logfile::Logfile(Inputs input, Report &report){
    logfilestream.open(input.get_logfile(), ofstream::app);
    cout << input.get_logfile();
    logfilestream << "year month day hour minute second";

}

void Logfile::write_logfile(Times time, Neutrals neutrals, Ions ions, Inputs inputs, Indices indices, Report report){
    vector<int> itime = time.get_iCurrent();
    vector<string> species_array = inputs.get_species_vector();
    if (!header){
        for (int i = 0; i < indices.all_indices_array_size(); i++){
            if (indices.get_nValues(i)>0){
                logfilestream << indices.get_name(i) << " ";
            }
        } 
        logfilestream << "min_temp max_temp mean_temp ";
        /*string name;
        for (int i = 0; i<species_array.size(); i++){
            name = species_array.at(i);
            logfilestream << name << "_min " << name << "_max " << name << "_mean ";
        }*/
        logfilestream << "\n";
        header = true;
    }
    for (int i = 0; i <= 6; ++i){
        logfilestream << itime.at(i) << " ";
    }
    for (int i = 0; i < indices.all_indices_array_size(); i++){
        logfilestream << indices.get_index(time.get_current(), i) << " ";
    }
    /*for(int i = 0; i<species_array.size(); i++){
        vector<double> temps = Neutral_Ion_Data(species_array.at(i), neutrals, ions, report);
        for (int i=0; i<temps.size(); i++){
        logfilestream << temps.at(i) << " ";
        }
    }*/
    logfilestream << "\n";
}

void Logfile::close_logfile(){
    logfilestream.close();
}

/*vector<double> Logfile::Neutral_Ion_Data(string name, Neutrals neutrals, Ions ions, Report report){
    int id = neutrals.get_species_id(name, report);
    bool IsNeutral = false;
    bool IsIon =  false;
    vector<double> idk (3,0.0);
    if (id > -1){
        IsNeutral = true;
    }
    else{
        id = ions.get_species_id(name, report);
    }
    if (id > -1){
        IsIon = true;
    }
    if (IsNeutral){
        return neutrals.temp_data(neutrals.temperature_scgc);
    }
    if (IsIon){
        return ions.temp_data(ions.temperature_scgc);
    }
    return idk;
    //idk for if it is not neutral or ion
}*/