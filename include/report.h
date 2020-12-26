// (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef AETHER_INCLUDE_REPORT_H_
#define AETHER_INCLUDE_REPORT_H_

#include <string>
#include <vector>

class Report {

public:

  // Functions:

  Report();
  void set_verbose(int input);
  void print(int iLevel, std::string output_string);
  int test_verbose(int iLevel);
  int get_verbose();
  void enter(std::string, int &iFunction);
  void exit(std::string);
  void times();
  
private:

  int iVerbose;  
  
  struct item_struct {

    std::string entry;
    int nTimes;
    float timing_total;
    unsigned long long timing_start;
    int iLevel;
    int iStringPosBefore;
    int iLastEntry;

  };
  
  std::vector<item_struct> entries;
  int nEntries;
  std::string current_entry;
  int iCurrentFunction = -1;

  std::string divider;
  int divider_length;
  
  int iLevel;

};

#endif // AETHER_INCLUDE_REPORT_H_
