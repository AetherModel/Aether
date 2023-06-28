// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_EUV_H_
#define INCLUDE_EUV_H_

/**************************************************************
 * \class Euv
 *
 * \brief Defines the Extreme Ultraviolet radiation above the atmosphere
 * 
 * The Euv class defines the EUV environment above the atmosphere. It
 * does this through the use of a CSV file that contains a bunch of
 * information. Namely:
 *  - Upper and lower wavelengths of bins in the EUV spectrum
 *  - How to relate the wavelengths to brightness (using a model like EUVAC)
 *  - Absorption, ionization, and dissociation cross sections for species
 *
 * \author Aaron Ridley
 *
 * \date 2021/03/28 
 *
 **************************************************************/

#include <vector>
#include <string>

class Euv {

public:

  /// number of wavelengths in spectrum: 
  int nWavelengths;

  // number of lines in the EUV CSV file:
  int nLines;

  /// struct to describe a single line in the EUV CSV file:
  struct waveinfotype {
    /// Key name of the line (describes wavelength or species acting upon):
    std::string name;
    /// If a cross-section, what is species becoming:
    std::string to;
    /// Type of cross-section (abs, ion, diss):
    std::string type;
    /// Unit of cross-section, wavelength, etc (in the row):
    std::string units;
    /// Any notes for the particular row:
    std::string note;
    /// The actual numerical values of the cross-section/wavelength/whatever
    std::vector<float> values;
  };

  /// All of the information in the EUV CSV file
  std::vector<waveinfotype> waveinfo;

  /// EUV Spectrum, lower wavelength of the bins:
  std::vector<float> wavelengths_short;
  
  /// EUV Spectrum, upper wavelength of the bins:
  std::vector<float> wavelengths_long;
  
  /// EUV Spectrum, energy of bin:
  std::vector<float> wavelengths_energy;

  /// EUV Spectrum, intensity flux of each spectral bin at 1 AU:
  std::vector<float> wavelengths_intensity_1au;

  /// EUV Spectrum, intensity flux of each spectral bin scaled to Earth:
  std::vector<float> wavelengths_intensity_top;

  /// EUVAC model linear coefficients (1):
  std::vector<float> euvac_f74113;
  /// EUVAC model linear coefficients (2):
  std::vector<float> euvac_afac;

  /// NEUVAC model linear coefficients (1-3):
  std::vector<float> neuvac_s1;
  std::vector<float> neuvac_s2;
  std::vector<float> neuvac_s3;

  /// NEUVAC model powers (1-2):
  std::vector<float> neuvac_p1;
  std::vector<float> neuvac_p2;

  /// NEUVAC model intercept:
  std::vector<float> neuvac_int;
  
  // --------------------------------------------------------------------
  // Functions:

  /**********************************************************************
     \brief Initialize the Euv class
     \param input info about how user has configured things
     \param report allow reporting to occur
   **/
  Euv(Inputs args, Report report);

  /**********************************************************************
     \brief Compute the EUV spectrum given F107 and F107a
     \param time The times within the model (dt is needed)
     \param indices Need the F107 and F107a
     \param report allow reporting to occur
   **/
  int euvac(Times time, Indices indices, Report &report);

  /**********************************************************************
     \brief Compute the EUV spectrum given F107 and F107a (new version)
     \param time The times within the model (dt is needed)
     \param indices Need the F107 and F107a
     \param report allow reporting to occur
   **/
  int neuvac(Times time, Indices indices, Report &report);

  /**********************************************************************
     \brief Scale the EUV spectrum given the star - planet distance
     \param planet needed to compute the star - planet distance
     \param time Needed to compute orbital position around star
     \param report allow reporting to occur
   **/
  int scale_from_1au(Planets planet, Times time, Report report);

  /**********************************************************************
     \brief Pairs rows in the EUV CSV file with neutral and ions

     Reads through each row in the EUV CSV file and figures out whether
     the row is abs, ion, diss, and then figures out which neutral it is
     acting on and which neutral or ion results from the action 
     (e.g. O + photon -> O+, identifies O as ionization "loss" and
     O+ as an ionization "source")

     \param neutrals Needs names of the neutrals, stores lines in Neutrals
     \param ions Needs names of the ions
     \param input info about how user has configured things
     \param report allow reporting to occur
   **/
  bool pair_euv(Neutrals &neutrals,
		Ions ions,
		Inputs input,
		Report report);

  /**********************************************************************
     \brief Check to see if internal state of class is ok
   **/
  
  bool is_ok();
  
private:

  /**********************************************************************
     \brief Read in the EUV CSV file

     Read in the EUV CSV file that describes all of the wavelengths and
     cross sections (and any other EUV - related things that are a 
     function of wavelength)

     \param input info about how user has configured things
     \param report allow reporting to occur
   **/
  bool read_file(Inputs args, Report report);

  /**********************************************************************
     \brief Interprets the EUV CSV rows and returns the relevant row

     Find the correct row in the EUV CSV file information, and return
     the values in that row.

     \param item The string value to search for in the first column
     \param item2 If not blank, the string value to search for in the 2nd col.
     \return values The values in the CSV row that matches the item (and item2) 
     \param report Allow reporting to occur
   **/
  bool slot_euv(std::string item,
		std::string item2,
		std::vector<float> &values,
		Report report);

  /// An internal variable to hold the state of the class
  bool IsOk;
};

#endif  // INCLUDE_EUV_H_
