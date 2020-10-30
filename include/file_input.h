// (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef AETHER_INCLUDE_FILE_INPUT_H_
#define AETHER_INCLUDE_FILE_INPUT_H_

#include <string>
#include <fstream>
#include <vector>

std::string find_next_hash(std::ifstream &file_ptr);
std::vector<std::vector<std::string>> read_csv(std::ifstream &file_ptr);

#endif // AETHER_INCLUDE_FILE_INPUT_H_
