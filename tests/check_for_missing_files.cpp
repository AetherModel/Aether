#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <map>

#include "../include/aether.h"
#include "mpi.h"

Times time;
Report report;

bool check_settings(Inputs input) {
    std::map<std::string, std::vector<std::string>> settings_2keys;
    std::vector<std::string> settings_1key;
    std::vector<std::string> missing_settings;

    settings_2keys = {{"Restart", {"OutDir", "InDir", "do", "dt"}}, 
                      {"Logfile", {"name", "dt", "append", "species"}}, 
                      {"GeoGrid", {"AltFile", "IsUniformAlt", "MinAlt", "dAlt", "MinLat", "MaxLat", 
                       "MinLon", "MaxLon"}}, {"Student", {"is", "name"}}, {"CubeSphere", {"is"}},
                      {"Euv", {"Model", "HeatingEfficiency", "dt", "File"}}, 
                      {"Debug", {"dt"}}, {"Outputs", {"type", "dt"}},
                      {"GeoBlockSize", {"nLons", "nLats", "nAlts", "nBlocksLon", "nBlocksLat"}}, 
                      {"Ensembles", {"nMembers"}},
                      {"Debug", {"iVerbose", "iProc"}}};

    settings_1key = {"BField", "Seed", "AuroraFile", "ChemistryFile", "CollisionsFile", "IndicesLookupFile", "OmniwebFiles",
                    "F107File", "Planet", "PlanetCharacteristicsFile", "PlanetSpeciesFile", "ElectrodynamicsFile", 
                    "DoCalcBulkIonTemp", "Perturb"};

    bool is_ok = true;
    for(std::pair<std::string, std::vector<std::string>> keys : settings_2keys)
        for(std::string key2 : keys.second)
            if(!input.check_settings(keys.first, key2))
                missing_settings.push_back("[" + keys.first + ", " + key2 + "]");
    
    for(std::string setting : settings_1key)
        if(!input.check_settings(setting))
            missing_settings.push_back("[" + setting + "]");

    if(!is_ok){
        std::string error = "Missing settings: " + 
        std::accumulate(missing_settings.begin(), missing_settings.end(), std::string(", "));
        throw std::invalid_argument(error);
        return false;
    }
    return true;
}

int main() {
    Inputs input(time, report);
    if (!input.is_ok())
      throw std::string("input initialization failed!");
    
    assert(check_settings(input));
}