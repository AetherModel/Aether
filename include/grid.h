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

  std::vector<arma_cube> rad_unit_vcgc;
  arma_cube gravity_potential_scgc;
  std::vector<arma_cube> gravity_vcgc;

  std::vector<arma_cube> cent_acc_vcgc;

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

  void calc_sza(Planets planet, Times time);
  void calc_gse(Planets planet, Times time);
  void calc_mlt();

  void calc_grid_spacing(Planets planet);
  void calc_alt_grid_spacing();
  void calc_lat_grid_spacing();
  void calc_long_grid_spacing();
  void fill_grid_radius(Planets planet);
  void calc_rad_unit(Planets planet);
  void calc_gravity(Planets planet);
  bool init_geo_grid(Quadtree quadtree,
		     Planets planet);
  void create_sphere_connection(Quadtree quadtree);
  void create_sphere_grid(Quadtree quadtree);
  void create_cubesphere_connection(Quadtree quadtree);
  void create_cubesphere_grid(Quadtree quadtree);
  void create_altitudes(Planets planet);
  void fill_grid_bfield(Planets planet);
  bool read_restart(std::string dir);
  bool write_restart(std::string dir);
  void report_grid_boundaries();
  void calc_cent_acc(Planets planet);

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

  /**
   * \brief Set the interpolation coefficients
   * \param Lons The longitude of points
   * \param Lats The latitude of points
   * \param Alts The altitude of points
   * \pre This instance is an geo grid
   * \pre Lons, Lats and Alts have the same size
   * \return true if the function succeeds, false if the instance is not a
   *         geo grid or the size of Lons, Lats and Alts are not the same.
   */
  bool set_interpolation_coefs(const std::vector<precision_t> &Lons,
                               const std::vector<precision_t> &Lats,
                               const std::vector<precision_t> &Alts);
  /**
   * \brief Create a map of geographic locations to data and do the interpolation
   * \param data The value at the positions of geoLon, geoLat, and geoAlt
   * \pre The size of the data should be the same as the geoLat/Lon/Alt_scgc
   * \return A vector of estimated value at the points set by the last
   *         set_interpolation_coefs function call if the function succeeds,
   *         an empty vector if the data is not the same size as the geo grid.
   */
  std::vector<precision_t> get_interpolation_values(const arma_cube &data) const;

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

  // The index and coefficient used for interpolation
  // Each point is processed by the function set_interpolation_coefs and stored
  // in the form of this structure.
  // If the point is out of the grid, in_grid = false and all other members are undefined
  struct interp_coef_t {
    // The point is inside the cube of [iRow, iRow+1], [iCol, iCol+1], [iAlt, iAlt+1]
    uint64_t iRow;
    uint64_t iCol;
    uint64_t iAlt;
    // The coefficients along row, column and altitude
    precision_t rRow;
    precision_t rCol;
    precision_t rAlt;
    // Whether the point is within this grid or not
    bool in_grid;
  };

  // Return the index of the last element that has altitude smaller than or euqal to the input
  uint64_t search_altitude(const precision_t alt_in) const;

  // Calculate the range of a spherical grid
  void get_sphere_grid_range(struct sphere_range &sr) const;
  // Calculate the range of a cubesphere grid
  void get_cubesphere_grid_range(struct cubesphere_range &cr) const;

  // Helper function for set_interpolation_coefs
  void set_interp_coef_sphere(const sphere_range &sr,
                              const precision_t lon_in,
                              const precision_t lat_in,
                              const precision_t alt_in);
  void set_interp_coef_cubesphere(const cubesphere_range &cr,
                                  const precision_t lon_in,
                                  const precision_t lat_in,
                                  const precision_t alt_in);

  // Processed interpolation coefficients
  std::vector<struct interp_coef_t> interp_coefs;
};

#endif  // INCLUDE_GRID_H_
