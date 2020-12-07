// (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include "../include/inputs.h"
#include "../include/indices.h"
#include "../include/read_f107_file.h"

Indices::Indices(Inputs args) {
  int iErr;
  std::string file;

  iErr = 0;

  // Read F10.7 file (if set):

  file = args.get_f107_file();
  if (file.length() > 0) {
    std::vector<double> time;
    std::vector<float> f107array;
    iErr = read_f107_file(file, time, f107array);
    if (iErr == 0)
      iErr = set_f107(time, f107array);
    else std::cout << "ERROR in reading f107 file!!!\n";
  }

}

int Indices::set_f107(std::vector<double> time,
		      std::vector<float> f107array) {

  int iErr = set_index(time, f107array, f107);

  // We want to set the 81-day average.  This is somewhat complicated,
  // since it seems like the f107 file does have exactly 24 hour
  // spaced data.  It seems like there are 3 points per day.  Let's
  // just ignore this fact.

  // Let's simply start at the start time and then progress forward 24
  // hours at a time

  double currenttime = time[0];
  long nTimes = time.size(), itime = 0;
  double eightone = 81.0 * 86400.0;

  ind_time_pair tmp;

  while (currenttime < time[nTimes-1]-eightone) {

    long isub = itime, nSubs = 0;
    double sumf107 = 0, sumtime = 0;

    while (time[isub] < currenttime + eightone) {
      sumf107 += f107array[isub];
      sumtime += time[isub];
      isub++;
      nSubs++;
    }

    tmp.time  = sumtime/nSubs;
    tmp.index = sumf107/nSubs;

    f107a.push_back(tmp);

    itime++;
    currenttime = time[itime];

  }

  // Let's hold the last 81 days constant, which means we just put one
  // single value at the end which is equal to the last average value:

  tmp.time = time[nTimes-1];
  // tmp.index is already set!!!
  f107a.push_back(tmp);

  return iErr;

}

float Indices:: get_f107(double time) {

  return get_index(time, f107);

}

float Indices:: get_f107a(double time) {

  return get_index(time, f107a);

}

float Indices::get_index(double time, std::vector<ind_time_pair> index) {

  long iLow, iMid, iHigh;

  iLow = 0;
  iHigh = index.size()-1;
  iMid = (iHigh+iLow)/2;

  while (iHigh-iLow > 1) {

    // cout << "top: " << iLow << " " << iMid << " " << iHigh << " " << index[iMid].time-time <<"\n";

    // Break if iMid <= time < iMid+1
    if (index[iMid].time == time) break;
    if (index[iMid].time <= time && index[iMid+1].time > time) break;
    // Upper Half:
    if (index[iMid].time < time) {
      iLow = iMid;
      iMid = (iHigh+iLow)/2;
    } else {
      iHigh = iMid;
      iMid = (iHigh+iLow)/2;
    }

    // cout << "bot: " << iLow << " " << iMid << " " << iHigh << " " << index[iMid].time-time <<"\n";

  }

  // At this point, time should be between iMid and iMid+1:

  float x = (time-index[iMid].time)/(index[iMid+1].time-index[iMid].time);
  float value = (1.0-x) * index[iMid].index + x * index[iMid+1].index;

  return value;

}

int Indices::set_index(std::vector<double> time,
		       std::vector<float> indexarray,
		       std::vector<ind_time_pair> &index) {

  int iErr;
  ind_time_pair tmp;

  if (time.size() != indexarray.size()) {
    std::cout << "In set_index. Size of time and index arrays don't match!\n";
    iErr = 1;
  } else {

    for (int i=0; i<time.size(); i++) {
      tmp.time = time[i];
      tmp.index = indexarray[i];
      index.push_back(tmp);
    }

  }

  return iErr;

}
