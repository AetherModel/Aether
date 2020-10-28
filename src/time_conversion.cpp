
#include <vector>
#include <iostream>
#include "../include/constants.h"
#include "../include/time_conversion.h"

// -----------------------------------------------------------------------------
// day_of_year: Convert from year, month, day to day of year
// -----------------------------------------------------------------------------

int day_of_year(int year, int month, int day) {

  int doy;

  doy = 0;

  for (int i=1; i<month; i++) {
    doy += days_of_month[i-1];
    if (year % 4 == 0 && i == 2) doy++;
  }
  doy += day;

  return doy;

}

// -----------------------------------------------------------------------------
// time_int_to_real: convert from year, month, day, hour, minute, second
//                   to seconds since reference date
// -----------------------------------------------------------------------------

double time_int_to_real(std::vector<int> itime) {

  int nYears = itime[0] - reference_year;
  int nLeaps = nYears/4;
  int nDays = day_of_year(itime[0], itime[1], itime[2])-1;

  double timereal =
    double(itime[6])/1000.0d + // milliseconds
    double(itime[5]) + // seconds
    double(itime[4]) * seconds_per_minute +
    double(itime[3]) * seconds_per_hour + 
    double(nDays+nLeaps) * seconds_per_day + 
    double(nYears) * seconds_per_year; 

  return timereal;
  
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------

void time_real_to_int(double timereal, std::vector<int> &itime) {

  int nYears = int(timereal/seconds_per_year);
  int nLeaps = nYears/4;
  int nDays = int((timereal - (double(nYears)*seconds_per_year))/seconds_per_day);

  // This can happen in some circumstances, like the first few days of the year:
  if (nDays < nLeaps) {
    nYears = (timereal - (double(nLeaps)*seconds_per_day))/seconds_per_year;
    nLeaps = nYears/4;
    nDays = (timereal - (double(nYears)*seconds_per_year))/seconds_per_day;
    // This should be very rare:
    if (nDays < nLeaps) {
      nYears = (timereal - (double(nLeaps)*seconds_per_day))/seconds_per_year;
      nLeaps = nYears/4;
      nDays = (timereal - (double(nYears)*seconds_per_year))/seconds_per_day;
    }
  }

  // Subtract off the leap days:
  nDays -= nLeaps;

  // Calculate how much time is left, after subtracting off years and days:
  double timeleft = timereal
    - double(nYears) * seconds_per_year
    - double(nDays + nLeaps) * seconds_per_day;

  // Calculate hours and subtract them:
  int nHours = timeleft/seconds_per_hour;
  timeleft -= double(nHours) * seconds_per_hour;

  // Calculate minutes and subtract them:
  int nMinutes = timeleft/seconds_per_minute;
  timeleft -= double(nMinutes) * seconds_per_minute;

  // Calculate milliseconds:
  int nSeconds = timeleft;
  int nMillis = int((timeleft-double(nSeconds))*1000.0);

  itime[0] = nYears + reference_year;

  // This is to get the month right.  We need to see if we are in a leap year.
  // If we are, then we need to make sure that Feb has one extra day. This is
  // sort of a hack so we can keep days_of_month as a constant.
  int nMonths = 1;
  int iLeap = 0;
  if (itime[0] % 4 == 0) iLeap=1;
  std::vector<int> add_on_days {0,iLeap,0,0,0,0,0,0,0,0,0,0};
  while (nDays > (days_of_month[nMonths-1]+add_on_days[nMonths-1])) {
    nDays = nDays - (days_of_month[nMonths-1]+add_on_days[nMonths-1]);
    nMonths++;
  }
  itime[1] = nMonths;
  itime[2] = nDays+1;
  itime[3] = nHours;
  itime[4] = nMinutes;
  itime[5] = nSeconds;
  itime[6] = nMillis;
  
}

// -----------------------------------------------------------------------------
// Convert from integer time to actual Julian Day
// -----------------------------------------------------------------------------

double time_int_to_jday(std::vector<int> itime) {

  // For this, we are going to use a relationship between our
  // time reference system and the Julian Day reference system.
  
  double our_time = time_int_to_real(itime);
  double our_time_in_days = our_time/seconds_per_day;
  double jd = julian_day_of_reference + our_time_in_days;
  
  return jd;

}

// -----------------------------------------------------------------------------
// testing
// -----------------------------------------------------------------------------

int test_time_routines() {

  std::vector<int> itime {1970, 1, 1, 0, 0, 0, 0};
  double timeout, timecheck;
  int i;
  int iErr;

  iErr = 0;
  
  timeout = time_int_to_real(itime);
  // This assumes a reference time of Jan. 1, 1965 0000 UT:
  timecheck = 1.5776640e+08;
  
  for (i=0; i<6; i++) std::cout << itime[i] << " ";
  std::cout << " --> " << timeout << " compares to : " << timecheck << "\n";
  if (abs(timecheck-timeout) > 1.0) {
    iErr = 1;
    std::cout << "Fails!!!\n";
  } else std::cout << "Passes!!!\n";

  time_real_to_int(timecheck, itime);
  for (i=0; i<6; i++) std::cout << itime[i] << " ";
  std::cout << "\n";

  itime[0]=2000;
  itime[1]=1;
  itime[2]=1;

  double jd_test = time_int_to_jday(itime);
  std::cout << "Test Julian Day = " << jd_test << "\n";
  std::cout << "Julian Day 2000 = " << j2000 << "\n";
  
  return iErr;
  
}



// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------

