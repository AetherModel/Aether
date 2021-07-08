// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_GRID_H_
#define INCLUDE_GRID_H_

// ----------------------------------------------------------------------------
// Grid class
// ----------------------------------------------------------------------------

class Grid {

public:

  int get_IsGeoGrid();
  void set_IsGeoGrid(int value);

  int64_t get_nPointsInGrid();

  int64_t get_nX();
  int64_t get_nY();
  int64_t get_nZ();

  int64_t get_nLons();
  int64_t get_nLats();
  int64_t get_nAlts();

  int64_t get_nGCs();

  // Armidillo Cube Versions:
  fcube geoLon_scgc, geoX_scgc;
  fcube geoLat_scgc, geoY_scgc;
  fcube geoAlt_scgc, geoZ_scgc;
  fcube geoLocalTime_scgc;
  
  // These define the magnetic grid:
  // Armidillo Cube Versions:
  fcube magLon_scgc, magX_scgc;
  fcube magLat_scgc, magY_scgc;
  fcube magAlt_scgc, magZ_scgc;
  fcube magLocalTime_scgc;

  // These are the locations of the magnetic poles:
  //  ll -> lat, lon, radius independent
  fvec mag_pole_north_ll;
  fvec mag_pole_south_ll;

  // pole gse -> needs to be for each altitude, so we can compute
  // magnetic local time. We want to use some GSE conversion function,
  // so this type has to a vector of fcubes:
  std::vector<fcube> mag_pole_north_gse;
  std::vector<fcube> mag_pole_south_gse;
  
  std::vector<fcube> GSE_XYZ_vcgc;

  std::string altitude_name = "Altitude";
  std::string altitude_unit = "meters";

  std::string longitude_name = "Longitude";
  std::string longitude_unit = "radians";

  std::string latitude_name = "Latitude";
  std::string latitude_unit = "radians";

  // These are derived variables from the grid:

  // Switch to armadillo variables (float cubes):
  fcube radius_scgc;
  fcube radius2_scgc;
  fcube radius2i_scgc;
  fcube gravity_scgc;

  fcube sza_scgc;
  fcube cos_sza_scgc;

  fcube dalt_center_scgc;
  fcube dalt_lower_scgc;

  std::vector<fcube> bfield_vcgc;
  fcube bfield_mag_scgc;

  Grid(int nX_in, int nY_in, int nZ_in, int nGCs_in);

  void calc_sza(Planets planet, Times time, Report &report);
  void calc_gse(Planets planet, Times time, Report &report);
  void calc_mlt(Report &report);
  void fill_grid(Planets planet, Report &report);
  void fill_grid_radius(Planets planet, Report &report);
  void init_geo_grid(Planets planet, Inputs input, Report &report);
  void fill_grid_bfield(Planets planet, Inputs input, Report &report);

  void init_mag_grid(Planets planet, Inputs input, Report &report);
 private:

  int IsGeoGrid;

  int64_t nX, nLons;
  int64_t nY, nLats;
  int64_t nZ, nAlts;

  int nGCs; // number of ghostcells

};

#endif  // INCLUDE_GRID_H_
