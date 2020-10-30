// (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <string>
#include <fstream>
#include <vector>
#include <sstream> 
#include <iostream>

// -------------------------------------------------------------------
// Read the file until it gets to a # as the first character
// -------------------------------------------------------------------

std::string find_next_hash(std::ifstream &file_ptr) {

  std::string line;

  if (!file_ptr.is_open()) {
    std::cout << "Could not open file!\n";
  } else {

    while (getline(file_ptr,line)) {
      if (line[0] == '#') return line;
      std::cout << line << "\n";
    }
    
  }

  line = "";
  return line;
}

// -------------------------------------------------------------------
// This can read a generic comma separated value file
// This is somewhat generic in that it can read from the current
// position to the first line that is blank, so it can be used
// multiple times per file.
// Returns a 2D array of strings.
// -------------------------------------------------------------------

std::vector<std::vector<std::string>> read_csv(std::ifstream &file_ptr) {

  std::vector<std::vector<std::string>> data;
  std::vector<std::string> row;

  std::string line, col;

  line = "  ";

  if (!file_ptr.is_open()) {
    std::cout << "File is not open!\n";
  } else {

    int IsFirstTime = 1;
    while (getline(file_ptr,line) && line.length() > 1) {
      std::stringstream ss(line);
      int j = 0;
      while (getline(ss, col, ',')) {
	std::cout << col << "\n";
	if (IsFirstTime) row.push_back(col);
	else row[j]=col;
	j++;
      }
      data.push_back(row);
    }
    
  }
  
  return data;

}

