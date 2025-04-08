// Copyright 2025, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "aether.h"

bool grid_match(Grid gGrid,
                Grid mGrid,
                Quadtree gQuadtree,
                Quadtree mQuadtree) {

  // Let's do magnetic to geographic first:

  int64_t iX, mnX = mGrid.get_nX();
  int64_t iY, mnY = mGrid.get_nY();
  int64_t iZ, mnZ = mGrid.get_nZ();
  int64_t mGCs = mGrid.get_nGCs();
  precision_t lon, lat;
  precision_t normX, normY, normZ;
  arma_vec norms(3);
  int64_t iNode;

  for (iX = mGCs; iX < mnX - mGCs; iX++) {
    for (iY = mGCs; iY < mnY - mGCs; iY++) {
      for (iZ = mGCs; iZ < mnZ - mGCs; iZ++) {
        lon = mGrid.geoLon_scgc(iX, iY, iZ);
        lat = mGrid.geoLat_scgc(iX, iY, iZ);
        if (gGrid.iGridShape_ == gGrid.iSphere_) {
          norms(0) = lon / cPI;
          norms(1) = lat / cPI;
          norms(2) = 0.0;
          iNode = gQuadtree.find_point(norms);
        } else {
          norms = sphere_to_cube(lon, lat);
          iNode = gQuadtree.find_point(norms);
        }
        if (report.test_verbose(6))
          std::cout << "lon, lat, node: " << lon*cRtoD << " "
                    << lat*cRtoD << " "
                    << norms(0) << " "
                    << norms(1) << " "
                    << norms(2) << " "
                    << iNode << "\n";
      }
    }
  }

  return true;
}
