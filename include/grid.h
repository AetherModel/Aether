// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_GRID_H_
#define INCLUDE_GRID_H_

#include "mpi.h"

// ----------------------------------------------------------------------------
// Grid class
// ----------------------------------------------------------------------------

class Grid {

public:

  // Armidillo Cube Versions:
  arma_cube geoLon_scgc, geoX_scgc;
  arma_cube geoLat_scgc, geoY_scgc;
  arma_cube geoAlt_scgc, geoZ_scgc;
  arma_cube geoLocalTime_scgc;

  arma_cube geoLon_Left;
  arma_cube geoLon_Down;
  arma_cube geoLon_Corner;

  arma_cube geoLat_Left;
  arma_cube geoLat_Down;
  arma_cube geoLat_Corner;

  arma_cube geoAlt_Below;
  arma_cube geoAlt_Corner;

  // These define the magnetic grid:
  // Armidillo Cube Versions:
  arma_cube magLon_scgc, magX_scgc;
  arma_cube magLat_scgc, magY_scgc;
  arma_cube magAlt_scgc, magZ_scgc;
  arma_cube magLocalTime_scgc;

  // These are the locations of the magnetic poles:
  //  ll -> lat, lon, radius independent
  arma_vec mag_pole_north_ll;
  arma_vec mag_pole_south_ll;

  // pole gse -> needs to be for each altitude, so we can compute
  // magnetic local time. We want to use some GSE conversion function,
  // so this type has to a vector of arma_cubes:
  std::vector<arma_cube> mag_pole_north_gse;
  std::vector<arma_cube> mag_pole_south_gse;

  std::vector<arma_cube> GSE_XYZ_vcgc;

  std::string altitude_name = "Altitude";
  std::string altitude_unit = "meters";

  std::string longitude_name = "Longitude";
  std::string longitude_unit = "radians";

  std::string latitude_name = "Latitude";
  std::string latitude_unit = "radians";

  // These are derived variables from the grid:

  // Switch to armadillo variables (precision_t cubes):
  arma_cube radius_scgc;
  arma_cube radius2_scgc;
  arma_cube radius2i_scgc;
  arma_cube gravity_scgc;

  arma_cube sza_scgc;
  arma_cube cos_sza_scgc;

  arma_cube dalt_center_scgc;
  arma_cube dalt_lower_scgc;
  arma_cube dalt_ratio_scgc;
  arma_cube dalt_ratio_sq_scgc;

  arma_cube dlon_center_scgc;
  arma_cube dlon_center_dist_scgc;

  arma_cube dlat_center_scgc;
  arma_cube dlat_center_dist_scgc;

  std::vector<arma_cube> bfield_vcgc;
  arma_cube bfield_mag_scgc;
  std::vector<arma_cube> bfield_unit_vcgc;

  Grid(int nX_in, int nY_in, int nZ_in, int nGCs_in);

  int get_IsGeoGrid();
  bool get_HasBField();
  void set_IsGeoGrid(int value);

  int64_t get_nPointsInGrid();

  int64_t get_nX();
  int64_t get_nY();
  int64_t get_nZ();

  int64_t get_nLons();
  int64_t get_nLats();
  int64_t get_nAlts();

  int64_t get_nGCs();

  void calc_sza(Planets planet, Times time, Report &report);
  void calc_gse(Planets planet, Times time, Report &report);
  void calc_mlt(Report &report);
  void fill_grid(Planets planet, Report &report);
  void fill_grid_radius(Planets planet, Report &report);
  bool init_geo_grid(Quadtree quadtree,
		     Planets planet,
		     Inputs input,
		     Report &report);
  void create_sphere_connection(Quadtree quadtree,
                                Inputs input,
                                Report &report);
  void create_sphere_grid(Quadtree quadtree, Inputs input, Report &report);
  void create_cubesphere_connection(Quadtree quadtree,
                                    Inputs input,
                                    Report &report);
  void create_cubesphere_grid(Quadtree quadtree, Inputs input, Report &report);
  void create_altitudes(Planets planet, Inputs input, Report &report);
  void fill_grid_bfield(Planets planet, Inputs input, Report &report);
  bool read_restart(std::string dir);
  bool write_restart(std::string dir);
  void report_grid_boundaries();

  // Need to move these to private at some point:

  bool IsLatLonGrid;
  bool IsCubeSphereGrid;
  bool DoesTouchNorthPole;
  bool DoesTouchSouthPole;
  /// The processor to the East/Right/X+:
  int iProcXp;
  /// The processor to the West/Left/X-:
  int iProcXm;
  /// The processor to the North/Up/Y+:
  int iProcYp;
  /// The processor to the South/Down/Y-:
  int iProcYm;

  arma_vec edge_Xp;
  arma_vec edge_Yp;
  arma_vec edge_Xm;
  arma_vec edge_Ym;

  int64_t iRoot;
  int64_t iRootXp;
  int64_t iRootXm;
  int64_t iRootYp;
  int64_t iRootYm;

  struct messages_struct {
    int64_t iFace;
    int64_t iProc_to;
    int64_t iSizeTotal;
    int64_t iTag;
    bool IsPole;
    bool DoReverseX;
    bool DoReverseY;
    bool XbecomesY;

    /// Variables needed for asynchronous message passing
    MPI_Request requests;
    precision_t* buffer;
    precision_t* rbuffer;
  };

  std::vector<messages_struct> interchanges;

  messages_struct make_new_interconnection(int64_t iDir,
					   int64_t nVars,
					   int64_t iProc_to,
					   arma_vec edge_center,
					   bool IsPole,
					   bool DoReverseX,
					   bool DoReverseY,
					   bool XbecomesY);

  bool send_one_face(int64_t iFace);
  bool receive_one_face(int64_t iFace);

  // interpolation functions
  // Estimate the value of the point at (lon_in, lat_in, alt_in)
  precision_t interp_linear(const arma_cube &data,
                            const precision_t lon,
                            const precision_t lat,
                            const precision_t alt);
  // the position of a vector of points is specified by
  // three vectors of lon, lat and alt
  std::vector<precision_t> interp_linear(const arma_cube &data,
                                         const std::vector<precision_t> &Lons,
                                         const std::vector<precision_t> &Lats,
                                         const std::vector<precision_t> &Alts);

 private:

  int IsGeoGrid;
  bool HasBField;

  int64_t nX, nLons;
  int64_t nY, nLats;
  int64_t nZ, nAlts;

  int nGCs; // number of ghostcells

  // interpolation members
  // The struct representing the range of a spherical grid
  struct sphere_range {
    precision_t lon_min;
    precision_t lon_max;
    precision_t dLon;
    precision_t lat_min;
    precision_t lat_max;
    precision_t dLat;
    precision_t alt_min;
    precision_t alt_max;
  };
  // The struct representing the range of a cubesphere grid
  struct cubesphere_range {
    // The minimum value and delta change of row and col
    // We don't use row_max and col_max because they are not promised to be
    // greater than min, for example the right norm of suface 2 expands along
    // the -x axis. drow and dcol can be negative, and boundary checking will
    // compare the theoretical index with 0 and nLon or nLat
    precision_t row_min;
    precision_t drow;
    precision_t col_min;
    precision_t dcol;
    // Range of altitude
    precision_t alt_min;
    precision_t alt_max;
    // The surface number of the grid
    int64_t surface_number;
    // The axis that row and col expands along
    // 0 means x-axis, 1 means y-axix, and 2 means z-axis
    int64_t row_direction;
    int64_t col_direction;
    // Used to promise that one and only one processor
    // returns the interpolation value and all others return cNinf
    bool row_min_exclusive;
    bool row_max_exclusive;
    bool col_min_exclusive;
    bool col_max_exclusive;
  };

  // Return the index of the last element that has altitude smaller than or euqal to the input
  uint64_t search_altitude(const precision_t alt_in);

  // Calculate the range of a spherical grid
  void get_sphere_grid_range(struct sphere_range &sr);
  // Calculate the range of a cubesphere grid
  void get_cubesphere_grid_range(struct cubesphere_range &cr);
  
  // Helper function for interpolation, so that grid range is only
  // calculated once no matter how many points
  precision_t interp_sphere_linear_helper(const arma_cube &data,
                                          const sphere_range &sr,
                                          const precision_t lon_in,
                                          const precision_t lat_in,
                                          const precision_t alt_in);
  precision_t interp_cubesphere_linear_helper(const arma_cube &data,
                                              const cubesphere_range &cr,
                                              const precision_t lon_in,
                                              const precision_t lat_in,
                                              const precision_t alt_in);
};

#endif  // INCLUDE_GRID_H_
