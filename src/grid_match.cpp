// Copyright 2025, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "aether.h"

// -----------------------------------------------------------------------------
// Send arrays of variables to other processors on the given grid.
// -----------------------------------------------------------------------------

bool exchange_information(int64_t *nPointsToPass,
                          std::vector<precision_t *> varToSend,
                          int64_t *nPointsToReceive,
                          std::vector<precision_t *> varToReceive) {

  int64_t jNode, iPt, iTag, iProcTo, iProcFrom;
  std::vector<MPI_Request> requests(nGrids);

  // Here we send the message into the wind:
  //   - if it is the same processor, just copy the information
  //   - if it is a different processor, send the data
  for (jNode = 0; jNode < nGrids ; jNode++) {
    if (jNode == iGrid) {
      for (iPt = 0; iPt < nPointsToPass[jNode]; iPt ++) {
        varToReceive[jNode][iPt] = varToSend[jNode][iPt];
      }
    } else {
      iProcTo = iMember * nGrids + jNode;
      // iTag is a unique id allowing all processors to
      // communicate asynchronously
      iTag = iProc * 10000 + iProcTo;
      MPI_Isend(varToSend[jNode],
                nPointsToPass[jNode] * sizeof(precision_t),
                MPI_BYTE,
                iProcTo,
                iTag,
                aether_comm,
                &requests[jNode]);
    }
  }

  // Wait for everyone to get the information that was sent:
  for (jNode = 0; jNode < nGrids ; jNode++)
    if (jNode != iGrid)
      MPI_Wait(&requests[jNode], MPI_STATUS_IGNORE);

  // Receive it into the receiving array:
  for (jNode = 0; jNode < nGrids ; jNode++)
    if (jNode != iGrid) {
      iProcFrom = iMember * nGrids + jNode;
      // Rebuid the unique id:
      iTag = iProcFrom * 10000 + iProc;
      MPI_Recv(varToReceive[jNode],
               nPointsToReceive[jNode] * sizeof(precision_t),
               MPI_BYTE,
               jNode,
               iTag,
               aether_comm,
               MPI_STATUS_IGNORE);
    }

  MPI_Barrier(aether_comm);
  return true;
}

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
  int64_t jNode, kNode;
  int64_t *nPointsToPass = static_cast<int64_t*>(malloc(nGrids * sizeof(int64_t)));
  int64_t *nPointsToReceive = static_cast<int64_t*>(malloc(nGrids * sizeof(int64_t)));
  int64_t *nPointsDummy = static_cast<int64_t*>(malloc(nGrids * sizeof(int64_t)));

  for (jNode = 0; jNode < nGrids ; jNode++)
    nPointsToPass[jNode] = 0;

  // This is not the most efficient way to do this, but the first pass, let's
  // just count how many points we need to send to the other processors:
  for (iX = mGCs; iX < mnX - mGCs; iX++) {
    for (iY = mGCs; iY < mnY - mGCs; iY++) {
      for (iZ = mGCs; iZ < mnZ - mGCs; iZ++) {
        lon = mGrid.geoLon_scgc(iX, iY, iZ);
        lat = mGrid.geoLat_scgc(iX, iY, iZ);
        if (gGrid.iGridShape_ == gGrid.iSphere_) {
          norms(0) = lon / cPI;
          norms(1) = lat / cPI;
          norms(2) = 0.0;
          jNode = gQuadtree.find_point(norms);
        } else {
          norms = sphere_to_cube(lon, lat);
          jNode = gQuadtree.find_point(norms);
        }
        if (jNode < 0 || jNode >= nGrids) {
            std::cout << "out of bounds!!! " << jNode << "\n";
        }
        nPointsToPass[jNode] = nPointsToPass[jNode]+1;
        /* std::cout << "lon, lat, node: " << lon*cRtoD << " "
            << lat*cRtoD << " "
            << norms(0) << " "
            << norms(1) << " "
            << norms(2) << " "
            << jNode << " "
            << iProc << " "
            << nPoints[jNode] << "\n"; */
      }
    }
  }
  std::cout << "made it here: " << iProc << "\n";
  MPI_Barrier(aether_comm);

  for (jNode = 0; jNode < nGrids ; jNode++)
    std::cout << "nPtsToPass : " << iProc << " " << nPointsToPass[jNode] << "\n";

  std::cout << "sending number of points :\n";

  // This section sends the number of points that need to be transfered to each processor.
  // Then the processor saves the number of points, so it can be remembered, and both the
  // sender and receiver will have the information.
  for (jNode = 0; jNode < nGrids ; jNode++) {
    if (jNode == iGrid) {
      for (kNode = 0; kNode < nGrids ; kNode++)
        nPointsDummy[kNode] = nPointsToPass[kNode];
    }
    MPI_Bcast(nPointsDummy, nGrids, MPI_INT64_T, jNode, aether_comm);
    nPointsToReceive[jNode] = nPointsDummy[iGrid];
  }

  MPI_Barrier(aether_comm);

  for (jNode = 0; jNode < nGrids ; jNode++) {
    std::cout << "nPtsToReceive : " << iProc << " " << jNode << " " << nPointsToReceive[jNode] << "\n";
    MPI_Barrier(aether_comm);
  }

  //  Now we need to create an array of send points and an array of receive points.
  std::vector<precision_t *> latsToPass(nGrids);
  std::vector<precision_t *> lonsToPass(nGrids);
  std::vector<precision_t *> altsToPass(nGrids);
  for (jNode = 0; jNode < nGrids ; jNode++) {
    latsToPass[jNode] = static_cast<precision_t*>(malloc(nPointsToPass[jNode] * sizeof(precision_t)));
    lonsToPass[jNode] = static_cast<precision_t*>(malloc(nPointsToPass[jNode] * sizeof(precision_t)));
    altsToPass[jNode] = static_cast<precision_t*>(malloc(nPointsToPass[jNode] * sizeof(precision_t)));
  }

  std::vector<precision_t *> latsToInterTo(nGrids);
  std::vector<precision_t *> lonsToInterTo(nGrids);
  std::vector<precision_t *> altsToInterTo(nGrids);
  for (jNode = 0; jNode < nGrids ; jNode++) {
    latsToInterTo[jNode] = static_cast<precision_t*>(malloc(nPointsToReceive[jNode] * sizeof(precision_t)));
    lonsToInterTo[jNode] = static_cast<precision_t*>(malloc(nPointsToReceive[jNode] * sizeof(precision_t)));
    altsToInterTo[jNode] = static_cast<precision_t*>(malloc(nPointsToReceive[jNode] * sizeof(precision_t)));
  }

  // now, the second pass, let's store the information so we can pass it:
  for (jNode = 0; jNode < nGrids ; jNode++)
    nPointsToPass[jNode] = 0;
  for (iX = mGCs; iX < mnX - mGCs; iX++) {
    for (iY = mGCs; iY < mnY - mGCs; iY++) {
      for (iZ = mGCs; iZ < mnZ - mGCs; iZ++) {
        lon = mGrid.geoLon_scgc(iX, iY, iZ);
        lat = mGrid.geoLat_scgc(iX, iY, iZ);
        if (gGrid.iGridShape_ == gGrid.iSphere_) {
          norms(0) = lon / cPI;
          norms(1) = lat / cPI;
          norms(2) = 0.0;
          jNode = gQuadtree.find_point(norms);
        } else {
          norms = sphere_to_cube(lon, lat);
          jNode = gQuadtree.find_point(norms);
        }
        latsToPass[jNode][nPointsToPass[jNode]] = lat;
        lonsToPass[jNode][nPointsToPass[jNode]] = lon;
        altsToPass[jNode][nPointsToPass[jNode]] = mGrid.geoAlt_scgc(iX, iY, iZ);
        nPointsToPass[jNode] = nPointsToPass[jNode]+1;
      }
    }
  }
  bool didWork;
  didWork = exchange_information(nPointsToPass,
                                 latsToPass,
                                 nPointsToReceive,
                                 latsToInterTo);


  return true;
}
