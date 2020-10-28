// (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef AETHER_INCLUDE_CONSTANTS_H_
#define AETHER_INCLUDE_CONSTANTS_H_

// Header file to define a bunch of physical constants

const float avogadros_number = 6.02214086e23; //
const float boltzmanns_constant = 1.38064852e-23; // m2 kg /s2 /K
const float mass_proton = 1.6726219e-27; // kg
const float amu = mass_proton;
const float mass_electron = 9.1094e-31; // Kg
const float planck_constant = 6.6261e-34; // Js
const float element_charge = 1.6022e-19; // C (J/eV)
const float speed_light = 2.9979e8; // m/s
const float univ_gas_constant = avogadros_number*boltzmanns_constant; // J/(moleK)
const float r_gas = univ_gas_constant*1.0E+07; // erg/(moleK)
const float pi = 3.141592653589793; 
const float twopi = 2*pi;
const float dtor = pi/180.0;
const float rtod = 180.0/pi;
const float p00 = 1.0e5; // Pa
const float T00 = 273.0; // K
const float SpecificHeatVolume = 3./2.; 
const float SpecificHeatPressure = SpecificHeatVolume + 1.0;
const float gamma_const = SpecificHeatPressure/SpecificHeatVolume;

// A bunch of constants for converting time between seconds and arrays:
const std::vector<int> days_of_month {31,28,31,30,31,30,31,31,30,31,30,31};
const int reference_year = 1965; // I think that this has to be the year after a leap year.
const double julian_day_of_reference = 2438762.0;
const double j2000 = 2451545.0;
const double seconds_per_year = 31536000.0;
const double seconds_per_day = 86400.0;
const double seconds_per_hour = 3600.0;
const double seconds_per_minute = 60.0;

// Conversion
const float mtokm = 1.0/1000.0;
const float kmtom = 1000.0;
const float pcm3topm3 = 1e6; // /cm3 to /m3
const float pcm2topm2 = 1e4; // /cm2 to /m2
const float pcmtopm = 100.0; // /cm to /m
const float atom = 1.0e-10; // angstrom to m

// Stellar constants:
const float au_to_m = 1.495978707e11; // m
const float gravitational_constant = 6.67408e-11; // m3/kg/s2
// This is for our sun:
const float solar_mass = 1.989e30; // kg

#endif // AETHER_INCLUDE_CONSTANTS_H_
