// Copyright 2025, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "aether.h"

// -----------------------------------------------------------------------------
//  Fill in XYZ in geo and mag coordinates
// -----------------------------------------------------------------------------

void Grid::calc_grid_spacing(Planets planet) {
  int64_t iLon, iLat, iAlt;

  report.print(3, "starting calc_grid_spacing");

  calc_alt_grid_spacing();
  calc_lat_grid_spacing();
  calc_long_grid_spacing();

  std::vector<arma_cube> lon_lat_radius;
  lon_lat_radius.push_back(geoLon_scgc);
  lon_lat_radius.push_back(geoLat_scgc);
  lon_lat_radius.push_back(radius_scgc);

  std::vector<arma_cube> xyz;
  xyz = transform_llr_to_xyz_3d(lon_lat_radius);
  geoX_scgc = xyz[0];
  geoY_scgc = xyz[0];
  geoZ_scgc = xyz[0];

  report.print(3, "ending calc_grid_spacing");
}

// ---------------------------------------
// Grid spacing for altitude:
// ---------------------------------------

void Grid::calc_alt_grid_spacing() {

  int64_t iAlt;

  report.print(4, "starting calc_alt_grid_spacing");

  for (iAlt = 1; iAlt < nAlts - 1; iAlt++) {
    dalt_center_scgc.slice(iAlt) =
      (geoAlt_scgc.slice(iAlt + 1) - geoAlt_scgc.slice(iAlt - 1)) / 2.0;
    dalt_lower_scgc.slice(iAlt) =
      geoAlt_scgc.slice(iAlt) - geoAlt_scgc.slice(iAlt - 1);
    dr_lower_scgc.slice(iAlt) =
      radius_scgc.slice(iAlt) - radius_scgc.slice(iAlt - 1);
  }

  dalt_center_scgc.slice(0) = dalt_center_scgc.slice(1);
  dalt_center_scgc.slice(nAlts - 1) = dalt_center_scgc.slice(nAlts - 2);

  dalt_lower_scgc.slice(0) = dalt_lower_scgc.slice(1);
  dr_lower_scgc.slice(0) = dr_lower_scgc.slice(1);
  iAlt = nAlts - 1;
  dalt_lower_scgc.slice(iAlt) =
    geoAlt_scgc.slice(iAlt) - geoAlt_scgc.slice(iAlt - 1);
  dr_lower_scgc.slice(iAlt) =
    radius_scgc.slice(iAlt) - radius_scgc.slice(iAlt - 1);

  // For a stretched grid, calculate some useful quantities:
  // lower is defined for the current cell, which
  // means that upper(iAlt) is lower(iAlt+1)
  // ratio = upper / lower
  for (iAlt = 0; iAlt < nAlts - 1; iAlt++)
    dalt_ratio_scgc.slice(iAlt) =
      dalt_lower_scgc.slice(iAlt + 1) / dalt_lower_scgc.slice(iAlt);

  iAlt = nAlts - 1;
  dalt_ratio_scgc.slice(iAlt) = dalt_ratio_scgc.slice(iAlt - 1);

  // Need the square of the ratio:
  dalt_ratio_sq_scgc = dalt_ratio_scgc % dalt_ratio_scgc;

  // Calculate the one-sided 3rd order gradient coefficients for the lower BC:

  arma_mat h1, h2, h3, h4, MeshH1, MeshH2, MeshH3, MeshH4;

  if (HasZdim & nAlts > nGCs + 5) {
    for (iAlt = 0; iAlt < nGCs; iAlt++) {
      h1 = dalt_lower_scgc.slice(iAlt + 2);
      h2 = dalt_lower_scgc.slice(iAlt + 3);
      h3 = dalt_lower_scgc.slice(iAlt + 4);
      h4 = dalt_lower_scgc.slice(iAlt + 5);
      MeshH1 = h1;
      MeshH2 = h1 + h2;
      MeshH3 = h1 + h2 + h3;
      MeshH4 = h1 + h2 + h3 + h4;
      MeshCoef1s3rdp1.slice(iAlt) =
        -1.0 * ( MeshH2 % MeshH3 % MeshH4 + MeshH1 % MeshH3 % MeshH4 +
                 MeshH1 % MeshH2 % MeshH4 + MeshH1 % MeshH2 % MeshH3) /
        (MeshH1 % MeshH2 % MeshH3 % MeshH4);
      MeshCoef1s3rdp2.slice(iAlt) =
        1.0 * ( MeshH2 % MeshH3 % MeshH4) / (h1 % h2 % (h2 + h3) % (h2 + h3 + h4));
      MeshCoef1s3rdp3.slice(iAlt) =
        -1.0 * ( MeshH1 % MeshH3 % MeshH4) / (MeshH2 % h2 % h3 % (h3 + h4));
      MeshCoef1s3rdp4.slice(iAlt) =
        1.0 * ( MeshH1 % MeshH2 % MeshH4) / (MeshH3 % (h3 + h2) % h3 % h4);
      MeshCoef1s3rdp5.slice(iAlt) =
        -1.0 * ( MeshH1 % MeshH2 % MeshH3) / (MeshH4 % (h2 + h3 + h4) % (h3 + h4) % h4);
    }
  }

  report.print(4, "ending calc_alt_grid_spacing");
  return;
}

// ---------------------------------------
// Grid spacing for latitude:
// ---------------------------------------

void Grid::calc_lat_grid_spacing() {

  int64_t iLat;

  report.print(4, "starting calc_lat_grid_spacing");

  for (iLat = 1; iLat < nLats - 1; iLat++) {
    dlat_center_scgc.col(iLat) =
      (geoLat_scgc.col(iLat + 1) - geoLat_scgc.col(iLat - 1)) / 2.0;
  }

  // Bottom (one sided):
  iLat = 0;
  dlat_center_scgc.col(iLat) =
    geoLat_scgc.col(iLat + 1) - geoLat_scgc.col(iLat);
  // Top (one sided):
  iLat = nLats - 1;
  dlat_center_scgc.col(iLat) =
    geoLat_scgc.col(iLat) - geoLat_scgc.col(iLat - 1);

  // Make this into a distance:
  dlat_center_dist_scgc = dlat_center_scgc % radius_scgc;
  report.print(4, "ending calc_lat_grid_spacing");
}

// ---------------------------------------
// Grid spacing for longitude:
// ---------------------------------------

void Grid::calc_long_grid_spacing() {

  int64_t iLon;

  report.print(4, "starting calc_long_grid_spacing");

  for (iLon = 1; iLon < nLons - 1; iLon++)
    dlon_center_scgc.row(iLon) =
      (geoLon_scgc.row(iLon + 1) - geoLon_scgc.row(iLon - 1)) / 2.0;

  // Bottom (one sided):
  iLon = 0;
  dlon_center_scgc.row(iLon) =
    geoLon_scgc.row(iLon + 1) - geoLon_scgc.row(iLon);
  // Top (one sided):
  iLon = nLons - 1;
  dlon_center_scgc.row(iLon) =
    geoLon_scgc.row(iLon) - geoLon_scgc.row(iLon - 1);

  // Make this into a distance:
  dlon_center_dist_scgc =
    dlon_center_scgc % radius_scgc % abs(cos(geoLat_scgc));

  report.print(4, "ending calc_long_grid_spacing");

}

// ---------------------------------------
// Grid spacing for magnetic longitude:
// ---------------------------------------

void Grid::calc_maglong_grid_spacing() {

  int64_t iLon;

  report.print(4, "starting calc_maglong_grid_spacing");

  for (iLon = 1; iLon < nLons - 1; iLon++)
    dmlon_center_scgc.row(iLon) =
      (magLon_scgc.row(iLon + 1) - magLon_scgc.row(iLon - 1)) / 2.0;

  // Bottom (one sided):
  iLon = 0;
  dmlon_center_scgc.row(iLon) =
    magLon_scgc.row(iLon + 1) - magLon_scgc.row(iLon);
  // Top (one sided):
  iLon = nLons - 1;
  dmlon_center_scgc.row(iLon) =
    magLon_scgc.row(iLon) - magLon_scgc.row(iLon - 1);

  // Make this into a distance:
  dmlon_center_dist_scgc =
    dmlon_center_scgc % radius_scgc % abs(cos(magLat_scgc));

  report.print(4, "ending calc_maglong_grid_spacing");
}

// -----------------------------------------------------------------------------
//  Calaculate Grid Spacing for Dipole Grid
// -----------------------------------------------------------------------------

void Grid::calc_dipole_grid_spacing(Planets planet) {

  int64_t iLon, iLat, iAlt;

  report.print(3, "starting calc_grid_spacing");

  // This is close, but may need to be adjusted later.
  // These quantities are obtained from integrating the scale factor (h)
  // The along-field-line distance (alt) should be right, but the lat distance
  // is the shortest distance from a point to the adjacent field line, not the adjacent cell.

  report.print(3, "starting alt");
  calc_alt_dipole_grid_spacing();
  report.print(3, "starting lat");
  calc_lat_dipole_grid_spacing();
  report.print(3, "starting long");
  calc_long_dipole_grid_spacing();

  std::vector<arma_cube> lon_lat_radius;
  lon_lat_radius.push_back(geoLon_scgc);
  lon_lat_radius.push_back(geoLat_scgc);
  lon_lat_radius.push_back(radius_scgc);
  std::vector<arma_cube> xyz;

  xyz = transform_llr_to_xyz_3d(lon_lat_radius);
  geoX_scgc = xyz[0];
  geoY_scgc = xyz[0];
  geoZ_scgc = xyz[0];

  report.print(3, "ending calc_grid_spacing");
}

// for sanity (only marginally helpful):
inline arma_mat delTm(arma_mat theta) {
  return (sqrt(3 * cos(theta) % cos(theta) + 1));
}
inline arma_cube delTc(arma_cube theta) {
  return (sqrt(3 * cos(theta) % cos(theta) + 1));
}

// -----------------------------------------------------------------------------
// Grid spacing for altitude:
//   - Dipole grid needs to be handled differently!
// -----------------------------------------------------------------------------

void Grid::calc_alt_dipole_grid_spacing() {

  int64_t iAlt;
  precision_t planetRadius;

  for (iAlt = 1; iAlt < nAlts - 1; iAlt++) {

    dalt_center_scgc.slice(iAlt) =
      abs(magAlt_scgc.slice(iAlt + 1) % sin(magLat_scgc.slice(iAlt + 1))
          % (1 / delTm(magLat_scgc.slice(iAlt + 1)))
          - magAlt_scgc.slice(iAlt - 1) % sin(magLat_scgc.slice(iAlt - 1))
          % (1 / delTm(magLat_scgc.slice(iAlt - 1)))) * 2;

    dalt_lower_scgc.slice(iAlt) =
      abs(magAlt_scgc.slice(iAlt) % sin(magLat_scgc.slice(iAlt))
          % (1 / delTm(magLat_scgc.slice(iAlt)))
          - magAlt_scgc.slice(iAlt - 1) % sin(magLat_scgc.slice(iAlt - 1))
          % (1 / delTm(magLat_scgc.slice(iAlt - 1)))) * 2;
    dr_lower_scgc.slice(iAlt) =
      radius_scgc.slice(iAlt) - radius_scgc.slice(iAlt - 1);
  }

  dalt_center_scgc.slice(0) = dalt_center_scgc.slice(1);
  dalt_center_scgc.slice(nAlts - 1) = dalt_center_scgc.slice(nAlts - 2);

  dalt_lower_scgc.slice(0) = dalt_lower_scgc.slice(1);
  dr_lower_scgc.slice(0) = dr_lower_scgc.slice(1);
  iAlt = nAlts - 1;
  dalt_lower_scgc.slice(iAlt) =
    geoAlt_scgc.slice(iAlt) - geoAlt_scgc.slice(iAlt - 1);
  dr_lower_scgc.slice(iAlt) =
    radius_scgc.slice(iAlt) - radius_scgc.slice(iAlt - 1);

  // For a stretched grid, calculate some useful quantities:
  // lower is defined for the current cell, which
  // means that upper(iAlt) is lower(iAlt+1)
  // ratio = upper / lower
  for (iAlt = 0; iAlt < nAlts - 1; iAlt++)
    dalt_ratio_scgc.slice(iAlt) =
      dalt_lower_scgc.slice(iAlt + 1) / dalt_lower_scgc.slice(iAlt);

  iAlt = nAlts - 1;
  dalt_ratio_scgc.slice(iAlt) = dalt_ratio_scgc.slice(iAlt - 1);

  // Need the square of the ratio:
  dalt_ratio_sq_scgc = dalt_ratio_scgc % dalt_ratio_scgc;
}

// ---------------------------------------
// Grid spacing for latitude:
// Again, different for the dipole...
//  - uhoh, might not be right. not actually perpendicular to q-p, but no way around that, i think.
// ---------------------------------------

void Grid::calc_lat_dipole_grid_spacing() {

  int64_t iLat;

  for (iLat = 1; iLat < nLats - 1; iLat++) {
    dlat_center_scgc.col(iLat) =
      abs(magAlt_scgc.col(iLat + 1) % sin(magLat_scgc.col(iLat + 1))
          % (1 / delTc(magLat_scgc.col(iLat + 1)))
          - magAlt_scgc.col(iLat - 1) % sin(magLat_scgc.col(iLat - 1))
          % (1 / delTc(magLat_scgc.col(iLat - 1)))) * 2;
  }

  // Bottom (one sided):
  iLat = 0;
  dlat_center_scgc.col(iLat) =
    geoLat_scgc.col(iLat + 1) - geoLat_scgc.col(iLat);
  // Top (one sided):
  iLat = nLats - 1;
  dlat_center_scgc.col(iLat) =
    geoLat_scgc.col(iLat) - geoLat_scgc.col(iLat - 1);

  // Make this into a distance:
  dlat_center_dist_scgc = dlat_center_scgc % radius_scgc;
}

// ---------------------------------------
// Grid spacing for longitude:
// ---------------------------------------

void Grid::calc_long_dipole_grid_spacing() {

  int64_t iLon;

  for (iLon = 1; iLon < nLons - 1; iLon++)
    dlon_center_scgc.row(iLon) =
      (magLon_scgc.row(iLon + 1) - magLon_scgc.row(iLon - 1)) / 2.0;

  // this might be fine for the dipole, if it works for the geo grid...

  // Bottom (one sided):
  iLon = 0;
  dlon_center_scgc.row(iLon) =
    magLon_scgc.row(iLon + 1) - magLon_scgc.row(iLon);
  // Top (one sided):
  iLon = nLons - 1;
  dlon_center_scgc.row(iLon) =
    magLon_scgc.row(iLon) - magLon_scgc.row(iLon - 1);

  // Make this into a distance:
  dlon_center_dist_scgc =
    // dlon_center_scgc % radius_scgc % abs(cos(geoLat_scgc));
    dlon_center_scgc % magAlt_scgc % cos(magLat_scgc);
}
