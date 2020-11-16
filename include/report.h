// (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef AETHER_INCLUDE_REPORT_H_
#define AETHER_INCLUDE_REPORT_H_

#include <string>

void report(int iLevel, std::string output_string, int iVerbose);
int test_verbose(int iLevel, int iVerbose);

#endif // AETHER_INCLUDE_REPORT_H_
