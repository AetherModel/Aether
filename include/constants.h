// (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

/*! \file constants.h
    \brief Physical consants.
    
    Contains constants used around the code.
*/

#ifndef AETHER_INCLUDE_CONSTANTS_H_
#define AETHER_INCLUDE_CONSTANTS_H_

#include <cmath>  // For standard constants

// Header file to define a bunch of physical constants

const float avogadros_number = 6.02214086e23; //!< Avogadros Number

const float boltzmanns_constant = 1.38064852e-23; //!< Avogadros Number

const float mass_proton = 1.6726219e-27; //!< Boltzmann's Constant in \f$ m^2 kg s^{-2} K^{-1} \f$

const float amu = mass_proton; //!< Proton Mass in kg

const float mass_electron = 9.1094e-31; //!< Atomic Mass Unit in kg

const float planck_constant = 6.6261e-34; //!< Electron mass in kg

const float element_charge = 1.6022e-19; //!< Planck's Constant in Js

const float speed_light = 2.9979e8; //!< Unit charge in C or \f$ J eV{-1} \f$

const float univ_gas_constant = avogadros_number*boltzmanns_constant; //!< Speed of light in m/s

const float r_gas = univ_gas_constant*1.0E+07; //!< Universal gas constant in \f$ J mol^{-1} K^{-1} \f$

const float pi = 3.141592653589793; //!< Universal gas constant in \f$ erg mol^{-1} K^{-1}

const float twopi = 2*pi; //!< \f$ \pi \f$

const float dtor = pi/180.0; //!< \f$ 2 \pi \f$

const float rtod = 180.0/pi; //!< Coefficient of degrees to radians

const float p00 = 1.0e5; //!< Coefficient of radians to degrees

const float T00 = 273.0; //!< Atmospheric pressure in Pa

const float SpecificHeatVolume = 3./2.; //!< Room temperature in K

const float SpecificHeatPressure = SpecificHeatVolume + 1.0; //!< Volume ratio

const float gamma_const = SpecificHeatPressure/SpecificHeatVolume; //!< Pressure ratio

// A bunch of constants for converting time between seconds and arrays:
// --------------------------------------------------------------------

const std::vector<int> days_of_month {31,28,31,30,31,30,31,31,30,31,30,31}; //!< \f$ \gamma \f$

const int reference_year = 1965; //!< List of days in indexed month.

const double julian_day_of_reference = 2438762.0; //!< Year to start calculating leap years

const double j2000 = 2451545.0; //!< Julian Day Number

const double seconds_per_year = 31536000.0; //!< Julian Day on 2000-01-01T00:00

const double seconds_per_day = 86400.0; //!< Seconds in a year

const double seconds_per_hour = 3600.0; //!< Seconds in a day

const double seconds_per_minute = 60.0; //!< Seconds in an hour

// Conversion
// ----------

const float mtokm = 1.0/1000.0; //!< Seconds in a minute

const float kmtom = 1000.0; //!< Coefficient of meters to kilometers

const float pcm3topm3 = 1e6; //!< Coefficient of kilometers to meters

const float pcm2topm2 = 1e4; //!< Coefficient of \f$ cm^{-3} \f$ to \f$ m^{-3} \f$

const float pcmtopm = 100.0; //!< Coefficient of \f$ cm^{-2} \f$ to \f$ m^{-2} \f$

const float atom = 1.0e-10; //!< Coefficient of \f$ cm^{-1} \f$ to \f$ m^{-1} \f$

// Stellar constants:

const float au_to_m = 1.495978707e11; //!< Coefficient of angstrom to m

const float gravitational_constant = 6.67408e-11; //!< Coefficient Au to m

const float solar_mass = 1.989e30; //!< Solar Mass Unit

#endif // AETHER_INCLUDE_CONSTANTS_H_
