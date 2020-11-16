// (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef AETHER_INCLUDE_READ_F107_FILE_H_
#define AETHER_INCLUDE_READ_F107_FILE_H_

#include <vector>
#include <string>

int read_f107_file(std::string f107_file,
		   std::vector<double> &time,
		   std::vector<float> &f107);

#endif // AETHER_INCLUDE_READ_F107_FILE_H_
