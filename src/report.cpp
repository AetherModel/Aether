
#include <string>
#include <iostream>

void report(int iLevel, std::string output_string, int iVerbose) {

  if (iLevel <= iVerbose) {

    for (int iL=0;iL<iLevel;iL++) std::cout << "=";
    std::cout << "> " << output_string << "\n";
    
  }

}

int test_verbose(int iLevel, int iVerbose) {

  int iPass = 0;
  if (iLevel <= iVerbose) {
    iPass = 1;
    for (int iL=0;iL<iLevel;iL++) std::cout << "=";
    std::cout << "> ";
  }

  return iPass;

}
