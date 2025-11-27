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
      for (iPt = 0; iPt < nPointsToPass[jNode]; iPt ++)
        varToReceive[jNode][iPt] = varToSend[jNode][iPt];
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

// -----------------------------------------------------------------------------
// This function:
//   on the requesting information side:
//     - figures out which processor each point of the other grid is on
//     - counts the points for each processor
//     - exchanges how many points to pass for each processor
//     - makes lists of coordinates to send to each processor
//     - sends those lists
//   on the interpolator side:
//     - builds interpolators for the requested information
// -----------------------------------------------------------------------------

bool grid_match(Grid &gGrid,
                Grid &mGrid,
                Quadtree gQuadtree,
                Quadtree mQuadtree) {

  std::string function = "grid_match";
  static int iFunction = -1;
  report.enter(function, iFunction);

  // Let's do magnetic to geographic first:

  int64_t iX, mnX = mGrid.get_nX();
  int64_t iY, mnY = mGrid.get_nY();
  int64_t iZ, mnZ = mGrid.get_nZ();
  int64_t mGCs = mGrid.get_nGCs();
  precision_t lon, lat;
  precision_t normX, normY, normZ;
  arma_vec norms(3);
  int64_t jNode, kNode;
  int64_t *nPointsToPass = static_cast<int64_t*>(malloc(nGrids * sizeof(
                                                          int64_t)));
  int64_t *nPointsToReceive = static_cast<int64_t*>(malloc(nGrids * sizeof(
                                                             int64_t)));
  int64_t *nPointsDummy = static_cast<int64_t*>(malloc(nGrids * sizeof(int64_t)));

  for (jNode = 0; jNode < nGrids ; jNode++)
    nPointsToPass[jNode] = 0;

  // This is not the most efficient way to do this, but the first pass, let's
  // just count how many points we need to send to the other processors:
  mGrid.gridToGridMap.set_size(mnX, mnY, mnZ);

  for (iX = mGCs; iX < mnX - mGCs; iX++) {
    for (iY = mGCs; iY < mnY - mGCs; iY++) {
      for (iZ = mGCs; iZ < mnZ - mGCs; iZ++) {
        lon = mGrid.geoLon_scgc(iX, iY, iZ);
        lat = mGrid.geoLat_scgc(iX, iY, iZ);

        if (gGrid.iGridShape_ == iSphere_) {
          norms(0) = lon / cPI;
          norms(1) = lat / cPI;
          norms(2) = 0.0;
          jNode = gQuadtree.find_point(norms);
        } else {
          norms = sphere_to_cube(lon, lat);
          jNode = gQuadtree.find_point(norms);
        }

        if (jNode < 0 || jNode >= nGrids)
          std::cout << "out of bounds!!! " << jNode << "\n";

        mGrid.gridToGridMap(iX, iY, iZ) = jNode;
        nPointsToPass[jNode] = nPointsToPass[jNode] + 1;
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

  MPI_Barrier(aether_comm);

  if (report.test_verbose(3)) {
    for (jNode = 0; jNode < nGrids ; jNode++)
      std::cout << "nPtsToPass : " << iProc << " " << nPointsToPass[jNode] << "\n";

    std::cout << "sending number of points :\n";
  }

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

  if (report.test_verbose(3)) {
    for (jNode = 0; jNode < nGrids ; jNode++)
      std::cout << "nPtsToReceive : " << iProc << " " << jNode << " " <<
                nPointsToReceive[jNode] << "\n";
  }

  //  Now we need to create an array of send points and an array of receive points.
  std::vector<precision_t *> latsToPass(nGrids);
  std::vector<precision_t *> lonsToPass(nGrids);
  std::vector<precision_t *> altsToPass(nGrids);
  std::vector<precision_t *> latsToInterTo(nGrids);
  std::vector<precision_t *> lonsToInterTo(nGrids);
  std::vector<precision_t *> altsToInterTo(nGrids);

  for (jNode = 0; jNode < nGrids ; jNode++) {
    latsToPass[jNode] = static_cast<precision_t*>(malloc(nPointsToPass[jNode] *
                                                         sizeof(precision_t)));
    lonsToPass[jNode] = static_cast<precision_t*>(malloc(nPointsToPass[jNode] *
                                                         sizeof(precision_t)));
    altsToPass[jNode] = static_cast<precision_t*>(malloc(nPointsToPass[jNode] *
                                                         sizeof(precision_t)));
    latsToInterTo[jNode] = static_cast<precision_t*>(malloc(
                                                       nPointsToReceive[jNode] * sizeof(precision_t)));
    lonsToInterTo[jNode] = static_cast<precision_t*>(malloc(
                                                       nPointsToReceive[jNode] * sizeof(precision_t)));
    altsToInterTo[jNode] = static_cast<precision_t*>(malloc(
                                                       nPointsToReceive[jNode] * sizeof(precision_t)));
  }

  // now, the second pass, let's store the information so we can pass it:
  for (jNode = 0; jNode < nGrids ; jNode++)
    nPointsToPass[jNode] = 0;

  for (iX = mGCs; iX < mnX - mGCs; iX++) {
    for (iY = mGCs; iY < mnY - mGCs; iY++) {
      for (iZ = mGCs; iZ < mnZ - mGCs; iZ++) {
        lon = mGrid.geoLon_scgc(iX, iY, iZ);
        lat = mGrid.geoLat_scgc(iX, iY, iZ);

        if (gGrid.iGridShape_ == iSphere_) {
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
        nPointsToPass[jNode] = nPointsToPass[jNode] + 1;
      }
    }
  }

  bool didWork;
  // Pass first coordinate (lons)
  didWork = exchange_information(nPointsToPass,
                                 lonsToPass,
                                 nPointsToReceive,
                                 lonsToInterTo);
  // Pass second coordinate (lats)
  didWork = exchange_information(nPointsToPass,
                                 latsToPass,
                                 nPointsToReceive,
                                 latsToInterTo);
  // Pass third coordinate (alts):
  didWork = exchange_information(nPointsToPass,
                                 altsToPass,
                                 nPointsToReceive,
                                 altsToInterTo);

  if (report.test_verbose(2)) {
    for (jNode = 0; jNode < nGrids ; jNode++) {
      std::cout << "Received the following points from iGrid = " << jNode << "\n";
      std::cout << " -> points received : " << nPointsToReceive[jNode] << "\n";

      for (int64_t iPt = 0; iPt < nPointsToReceive[jNode]; iPt++)
        std::cout << "  -> " << iPt << " "
                  << lonsToInterTo[jNode][iPt] << " "
                  << latsToInterTo[jNode][iPt] << " "
                  << altsToInterTo[jNode][iPt] << "\n";
    }
  }

  struct grid_to_grid_t oneGrid;

  int64_t nPts;

  for (jNode = 0; jNode < nGrids ; jNode++) {
    // These are backwards now, since we will switch sender and reciever:
    oneGrid.nPts = nPointsToReceive[jNode];
    oneGrid.nPtsReceive = nPointsToPass[jNode];
    oneGrid.iProcTo = iMember * nGrids + jNode;

    if (report.test_verbose(2))
      std::cout << "Making interpolation coefficients for : " << jNode
                << "; points : " << oneGrid.nPts << "\n";

    if (oneGrid.nPts > 0) {
      // Interpolation function takes vectors,
      // so transfer these arrays to vectors:
      std::vector<precision_t> Lons(oneGrid.nPts);
      std::vector<precision_t> Lats(oneGrid.nPts);
      std::vector<precision_t> Alts(oneGrid.nPts);

      for (int64_t iPt = 0; iPt < oneGrid.nPts; iPt++) {
        Lons[iPt] = lonsToInterTo[jNode][iPt];
        Lats[iPt] = latsToInterTo[jNode][iPt];
        Alts[iPt] = altsToInterTo[jNode][iPt];
      }

      oneGrid.interpCoefs = gGrid.get_interpolation_coefs(Lons, Lats, Alts);
    }

    gGrid.gridToGridCoefs.push_back(oneGrid);
  }

  report.exit(function);
  return didWork;
}

bool get_data_from_other_grid(Grid &gGrid,
                              Grid &mGrid,
                              arma_cube &gData,
                              arma_cube &mData) {

  std::string function = "get_data_from_other_grid";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int64_t jNode, iPt;
  std::vector<precision_t *> dataToSend(nGrids);
  std::vector<precision_t *> dataToReceive(nGrids);
  int64_t *nPointsToSend = static_cast<int64_t*>(malloc(nGrids * sizeof(
                                                          int64_t)));
  int64_t *nPointsToReceive = static_cast<int64_t*>(malloc(nGrids * sizeof(
                                                             int64_t)));

  for (jNode = 0; jNode < nGrids ; jNode++) {
    if (report.test_verbose(2))
      std::cout << "nPts : " << jNode << " " << gGrid.gridToGridCoefs[jNode].nPts <<
                "\n";

    nPointsToSend[jNode] = gGrid.gridToGridCoefs[jNode].nPts;
    nPointsToReceive[jNode] = gGrid.gridToGridCoefs[jNode].nPtsReceive;
    dataToSend[jNode] = static_cast<precision_t*>(malloc(
                                                    gGrid.gridToGridCoefs[jNode].nPts * sizeof(precision_t)));
    dataToReceive[jNode] = static_cast<precision_t*>(malloc(
                                                       gGrid.gridToGridCoefs[jNode].nPtsReceive * sizeof(precision_t)));
    std::vector<precision_t> values = gGrid.get_interpolation_values(gData,
                                      gGrid.gridToGridCoefs[jNode].interpCoefs);

    for (iPt = 0; iPt < gGrid.gridToGridCoefs[jNode].nPts; iPt++) {
      dataToSend[jNode][iPt] = values[iPt];

      if (report.test_verbose(2))
        std::cout << "datatosend : " << iPt << " " << dataToSend[jNode][iPt] << "\n";
    }
  }

  bool didWork = exchange_information(nPointsToSend,
                                      dataToSend,
                                      nPointsToReceive,
                                      dataToReceive);
  int64_t iX, mnX = mGrid.get_nX();
  int64_t iY, mnY = mGrid.get_nY();
  int64_t iZ, mnZ = mGrid.get_nZ();
  int64_t mGCs = mGrid.get_nGCs();
  std::vector<int64_t> iCounter(nGrids);

  for (jNode = 0; jNode < nGrids ; jNode++)
    iCounter[jNode] = 0;

  for (iX = mGCs; iX < mnX - mGCs; iX++) {
    for (iY = mGCs; iY < mnY - mGCs; iY++) {
      for (iZ = mGCs; iZ < mnZ - mGCs; iZ++) {
        jNode = mGrid.gridToGridMap(iX, iY, iZ);

        if (report.test_verbose(2)) {
          std::cout << "unpacking point : " << iX << " " << iY << " " << iZ << " " <<
                    jNode << " "
                    << iCounter[jNode] << " " << dataToReceive[jNode][iCounter[jNode]] << "\n";
        }

        mData(iX, iY, iZ) = dataToReceive[jNode][iCounter[jNode]];
        iCounter[jNode] = iCounter[jNode] + 1;
      }
    }
  }

  report.exit(function);
  return true;

}