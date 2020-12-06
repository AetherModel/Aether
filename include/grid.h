// (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef AETHER_INCLUDE_GRID_H_
#define AETHER_INCLUDE_GRID_H_

#include "inputs.h"
#include "sizes.h"
#include "planets.h"
#include "times.h"

// We need a naming convention for the variables that are defined on
// the grid.  These could then match the formulas that are used to find
// the given points on the grid.
/*

  _ - delimiter between main variable name and descriptor
  1 - indication of whether variable is scalar (s) or vector (v)
  2 - physical dimensions:
      3 - 3d (e.g., lon, lat, alt or x, y, z)
      


For example:
  _s3gc : scalar variable, 3d, include ghost cells, cell centers
  _31ne : 


 */

// These are mapping functions from 3d to 1d arrays. They have to be
// pretty precise, which is a bit scary to me.
//
// Assume [Lon][Lat][Alt] layout

// Scalars, 3D, Include Ghostcells, Cell Centers:
#define ijk_geo_s3gc(i,j,k) \
  ((i)*long(nGeoLatsG)*long(nGeoAltsG) + \
   (j)*long(nGeoAltsG) + \
   (k))
#define ijk_mag_s3gc(i,j,k) \
  ((i)*long(nMagLatsG)*long(nGeoAltsG) + \
   (j)*long(nMagAltsG) + \
   (k))

// Scalars, 3D, Include Ghostcells, Cell Edges (Altitude):
#define ijk_geo_s3ge3(i,j,k) \
  ((i)*long(nGeoLatsG)*long(nGeoAltsG+1) + \
   (j)*long(nGeoAltsG+1) + \
   (k));

// Vectors, 3D, Include Ghostcells, Cell Centers:
#define ijk_geo_v3gc(i,j,k,l) \
  ((i)*long(nGeoLatsG)*long(nGeoAltsG)*long(3) + \
   (j)*long(nGeoAltsG)*long(3) + \
   (k)*long(3) + \
   (l))


class Grid {

public:

  int get_IsGeoGrid();
  void set_IsGeoGrid(int value);

  // These define the geographic grid:
  float *geoLon_s3gc, *geoX_s3gc;
  float *geoLat_s3gc, *geoY_s3gc;
  float *geoAlt_s3gc, *geoZ_s3gc;
  
  // These define the magnetic grid:
  float *magLon_s3gc, *magX_s3gc;
  float *magLat_s3gc, *magY_s3gc;
  float *magAlt_s3gc, *magZ_s3gc;
  float *magLocalTime_s3gc;

  std::string altitude_name = "Altitude";
  std::string altitude_unit = "meters";

  std::string longitude_name = "Longitude";
  std::string longitude_unit = "radians";

  std::string latitude_name = "Latitude";
  std::string latitude_unit = "radians";
  
  // These are derived variables from the grid:
  float *radius_s3gc;
  float *radius_sq_s3gc;
  float *radius_inv_sq_s3gc;
  float *gravity_s3gc;
  float *sza_s3gc, *cos_sza_s3gc;

  float *geoX_edges_s3ge, *geoX_cell_width_s3gc;

  float *alt_edges;
  float *alt_cell_width;
  float *dalt_center_s3gc;
  float *dalt_lower_s3gc;

  Grid(int nX, int nY, int nZ);

  void calc_sza(Planets planet, Times time, Report &report);
  void fill_grid(Planets planet, Report &report);
  void fill_grid_radius(Planets planet, Report &report);
  void init_geo_grid(Planets planet, Inputs input, Report &report);
  
 private:

  int IsGeoGrid;

};

#endif // AETHER_INCLUDE_GRID_H_
