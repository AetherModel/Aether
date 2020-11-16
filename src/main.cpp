// (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "../include/times.h"
#include "../include/inputs.h"
#include "../include/euv.h"

int main() {

  int iErr = 0;

  Times time;
  Inputs input(time);
  Euv euv(input);
  Planets planet(input);
  Indices indices(input);
  
  return iErr;

}
