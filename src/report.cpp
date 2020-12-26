
#include <string>
#include <iostream>
#include <chrono>

#include "../include/report.h"


// -----------------------------------------------------------------------
// 
// -----------------------------------------------------------------------

Report::Report() {

  current_entry = "";
  nEntries = 0;
  iVerbose = 0;
  divider = ">";
  divider_length = divider.length();
  iLevel = 0;  
  
}

// -----------------------------------------------------------------------
// 
// -----------------------------------------------------------------------

void Report::enter(std::string input, int &iFunction) {

  int iOldStrLen = current_entry.length();

  current_entry = current_entry + divider + input;

  int iEntry = -1;

  if (iFunction > -1)
    if (current_entry == entries[iFunction].entry) 
      iEntry = iFunction;
  if (iEntry == -1) {
    for (int i=0; i < nEntries; i++) 
      if (current_entry == entries[i].entry) iEntry = i;
  }
  if (iEntry == -1) {
    item_struct tmp;
    tmp.entry = current_entry;
    tmp.nTimes = 0;
    tmp.timing_total = 0.0;
    tmp.iStringPosBefore = iOldStrLen;
    tmp.iLastEntry = iCurrentFunction;
    entries.push_back(tmp);
    nEntries++;
    iEntry = nEntries-1;
    iFunction = iEntry;
  }

  // This was taken from
  // https://stackoverflow.com/questions/19555121/how-to-get-current-timestamp-in-milliseconds-since-1970-just-the-way-java-gets
  unsigned long long now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();

  entries[iEntry].timing_start = now;
  iLevel++;
  entries[iEntry].iLevel = iLevel;
  iCurrentFunction = iEntry;
  
  print(iLevel,"Entering function : "+current_entry);
  
}

// -----------------------------------------------------------------------
// 
// -----------------------------------------------------------------------

void Report::exit(std::string input) {

  int iEntry = -1;
  iEntry = iCurrentFunction;
  //for (int i=0; i < nEntries; i++) 
  //  if (current_entry == entries[i].entry) iEntry = i;
  //std::cout << "iEntry : " << iEntry << " " << iCurrentFunction << " " << nEntries << "\n";
  if (iEntry == -1) {
    std::cout << "Report::exit Error!!! Could not find valid entry!\n";
    std::cout << "current_entry : " << current_entry << "\n";
  } else {
    unsigned long long now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    entries[iEntry].timing_total = entries[iEntry].timing_total +
      float(now - entries[iEntry].timing_start)/1000.0;
    entries[iEntry].nTimes++;

    // std::size_t pos = current_entry.find(input);
    // if (pos == std::string::npos) {
    //   std::cout << "Report::exit Error!!! Could not find input : " << input << "!\n";
    //   std::cout << "   current_entry : " << current_entry << "\n";
    // } else {
    // current_entry = current_entry.substr(0,pos-divider_length);
    current_entry = current_entry.substr(0,entries[iEntry].iStringPosBefore);
    iCurrentFunction = entries[iEntry].iLastEntry;
    iLevel--;
    // }
  }
  
}

// -----------------------------------------------------------------------
// 
// -----------------------------------------------------------------------

void Report::times() {

  std::cout << "Timing Summary :\n";
  for (int i=0; i < nEntries; i++) {
    std::cout << entries[i].entry << "\n";
    for (int j=0; j < entries[i].iLevel; j++) std::cout << "  ";
    std::cout << "nTimes called : " << entries[i].nTimes << "\n";
    for (int j=0; j < entries[i].iLevel; j++) std::cout << "  ";
    std::cout << "timing_total (s) : " << entries[i].timing_total << "\n";
  }

}

// -----------------------------------------------------------------------
// 
// -----------------------------------------------------------------------

void Report::print(int iLevel, std::string output_string) {

  if (iLevel <= iVerbose) {

    for (int iL=0;iL<iLevel;iL++) std::cout << "=";
    std::cout << "> " << output_string << "\n";
    
  }

}

// -----------------------------------------------------------------------
// 
// -----------------------------------------------------------------------

int Report::test_verbose(int iLevel) {

  int iPass = 0;
  if (iLevel <= iVerbose) {
    iPass = 1;
    for (int iL=0;iL<iLevel;iL++) std::cout << "=";
    std::cout << "> ";
  }

  return iPass;

}

// -----------------------------------------------------------------------
// 
// -----------------------------------------------------------------------

void Report::set_verbose(int input) {
  iVerbose = input;
}

// -----------------------------------------------------------------------
// 
// -----------------------------------------------------------------------

int Report::get_verbose() {
  return iVerbose;
}
