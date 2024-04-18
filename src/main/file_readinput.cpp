
// (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

/*
  This is an example code to demonstrate how to use the json library.

  Here is how I compiled the code:
  g++ newMain.cpp -o json_test.exe -std=c++11 -I /mnt/c/Users/rutvi/Documents/json-develop/include
 */

#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

int main() {
  std::ifstream infile_ptr("../inputs/aether_in.json");
  json commands;
  infile_ptr >> commands;

  //this is how to access the values attached to the main keys (commands)
  std::cout << "Read the ather input json file!\n\n";
  std::cout << "Debug: " << commands.at("Debug") << "\n";
  std::cout << "Planet: " << commands.at("Planet") << "\n";
  std::cout << "Starttime: " << commands.at("Starttime") << "\n";
  std::cout << "Endtime: " << commands.at("Endtime") << "\n";
  std::cout << "F107File: " << commands.at("F107File") << "\n";
  std::cout << "BField: " << commands.at("BField") << "\n";
  std::cout << "Output: " << commands.at("Output") << "\n" << "\n";

  //individual elements in an array can also be accessed
  //for example, the individual values in the starttime array
  std::cout << "Starttime Year: " << commands.at("Starttime").at(0) << "\n";
  infile_ptr.close();

}