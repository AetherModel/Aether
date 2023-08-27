// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "aether.h"

// -----------------------------------------------------------------------------
// This is where all of the exchange messages routines will sit.
//
// Notes:
//   - We are going to try to do asynchronous communications, which means we
//     need to do the following:
//      - pack all variables for all four faces
//      - send all messages
//      - receive all messages
//      - unpack all four faces
//   - To do this, we need to make send and receive buffers which can't be
//     touched until everything is complete, so we will build a structure
//     that contains both the send and receive buffers.
//
//   - Direction standard:
//     iDir == 0 => face 0 => right
//     iDir == 1 => face 1 = up
//     iDir == 2 => face 2 = left
//     iDir == 3 => face 3 = down
//     This is the side we are dealing with for the process.  For example,
//       iDir == 0 in sending could be iDir == 2 in receiving for blocks
//       near the equator.
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// Pack variables for message passing
//   value is variable to pack
//   packed is the buffer to pack it in
//   iCounter is the location in the buffer to start
//   nG is the number of ghost cells
//   iDir is the direction of the message pass:
//     0 - right
//     1 - top
//     2 - left
//     3 - bottom
//
// cells (assume gc = 2):
// 0  1 | 2 3 4 ... n-gc*2 n-gc-1 | n-gc n-1
//                       0      1 |    2   3 4  ... n-gc*2 n-gc-1 | n-gc n-1
// -----------------------------------------------------------------------------

bool pack_border(const arma_cube &value,
                 precision_t *packed,
                 int64_t *iCounter,
                 int64_t nG,
                 int iDir) {

  bool DidWork = true;
  static int64_t nX = value.n_rows;
  static int64_t nY = value.n_cols;
  static int64_t nZ = value.n_slices;

  int64_t iXstart, iXend;
  int64_t iYstart, iYend;

  // ----------------------------
  // left / right message passing
  if (iDir == 0 || iDir == 2) {
    iYstart = nG;
    iYend = nY - nG;

    if (iDir == 0) {
      // left -> right
      iXstart = nX - 2 * nG;
      iXend = nX - nG;
    } else {
      // right -> left
      iXstart = nG;
      iXend = 2 * nG;
    }
  }

  // ----------------------------
  // top / bottom message passing
  if (iDir == 1 || iDir == 3) {
    iXstart = nG;
    iXend = nX - nG;

    if (iDir == 1) {
      // bottom -> top
      iYstart = nY - 2 * nG;
      iYend = nY - nG;
    } else {
      // top -> bottom
      iYstart = nG;
      iYend = 2 * nG;
    }
  }

  try {
    for (int64_t iZ = nG; iZ < nZ - nG; iZ++) {
      for (int64_t iY = iYstart; iY < iYend; iY++) {
        for (int64_t iX = iXstart; iX < iXend; iX++) {
          packed[*iCounter] = value(iX, iY, iZ);
          *iCounter = *iCounter + 1;
        }
      }
    }
  } catch (...) {
    DidWork = false;
  }

  return DidWork;
}

// -----------------------------------------------------------------------------
// Unpack variables after message passing
//   value is variable to unpack into
//   packed is the buffer that contains the values
//   iCounter is the location in the buffer to start
//   nG is the number of ghost cells
//   iDir is the direction of the message pass:
//     0 - right
//     1 - top
//     2 - left
//     3 - bottom
//   DoReverseX and DoReverseY are because packing always happens from
//     lower left to upper right, while face we are unpacking too may
//     have a different (left - right and up - down) geometry
//   XbecomesY as above, but in this case, the L-R and U-D could change.
//     This is really for the CubeSphere dealing with the top/bottom to
//     sides geometry.
// -----------------------------------------------------------------------------

bool unpack_border(arma_cube &value,
                   precision_t *packed,
                   int64_t *iCounter,
                   int64_t nG,
                   int iDir,
                   bool DoReverseX,
                   bool DoReverseY,
                   bool XbecomesY) {

  bool DidWork = true;
  static int64_t nX = value.n_rows;
  static int64_t nY = value.n_cols;
  static int64_t nZ = value.n_slices;

  int64_t iXstart, iXend;
  int64_t iYstart, iYend;
  int64_t xInc = 1, yInc = 1;

  int64_t iXOff = 0;
  int64_t nCx = nX - 2 * nG;

  // This is for over the pole message passing on one processor:
  if (nProcs == 1 && (iDir == 1 || iDir == 3)) {
    if (nCx % 2 > 0) {
      std::cout << "If you are running on one processor, and you want to\n";
      std::cout << "model the whole Earth, it is highly recommended that\n";
      std::cout << "nLons is EVEN, so that the message passing across the\n";
      std::cout << "pole works as it should.\n";
    }

    iXOff = nCx / 2;
  }

  // ----------------------------
  // left / right message passing
  if (iDir == 0 || iDir == 2) {
    iYstart = nG;
    iYend = nY - nG;

    if (iDir == 2) {
      // left -> right
      iXstart = 0;
      iXend = nG;
    } else {
      // right -> left
      iXstart = nX - nG;
      iXend = nX;
    }
  }

  // ----------------------------
  // top / bottom message passing
  if (iDir == 1 || iDir == 3) {
    iXstart = nG;
    iXend = nX - nG;

    if (iDir == 3) {
      // bottom -> top
      iYstart = 0;
      iYend = nG;
    } else {
      // top -> bottom
      iYstart = nY - nG;
      iYend = nY;
    }
  }

  try {
    int64_t iXp, iYp;

    for (int64_t iZ = nG; iZ < nZ - nG; iZ++) {
      if (XbecomesY) {
        for (int64_t iX = iXstart; iX < iXend; iX += xInc) {
          iXp = iX;

          if (DoReverseX)
            iXp = iXend - 1 - (iX - iXstart);

          for (int64_t iY = iYstart; iY < iYend; iY++) {
            iYp = iY;

            if (DoReverseY)
              iYp = iYend - 1 - (iY - iYstart);

            value(iXp, iYp, iZ) = packed[*iCounter];
            *iCounter = *iCounter + 1;
          }
        }
      } else {
        for (int64_t iY = iYstart; iY < iYend; iY++) {
          iYp = iY;

          if (DoReverseY)
            iYp = iYend - 1 - (iY - iYstart);

          for (int64_t iX = iXstart; iX < iXend; iX += xInc) {
            iXp = iX;

            if (DoReverseX)
              iXp = iXend - 1 - (iX - iXstart);

            if (iXOff > 0) {
              iXp = (iXp + iXOff) % nCx;

              if (iXp < nG)
                iXp += nCx;
            }

            value(iXp, iYp, iZ) = packed[*iCounter];
            *iCounter = *iCounter + 1;
          }
        }
      }
    }
  } catch (...) {
    DidWork = false;
  }

  return DidWork;
}

// -----------------------------------------------------------------------------
// Pack one face for the NEUTRALS (den, temp, vel)
// -----------------------------------------------------------------------------

bool Neutrals::pack_one_face(int iReceiver,
                             precision_t *buffer,
                             int nG, int iDir,
                             bool IsPole) {

  bool DidWork = true;
  int64_t iP;
  MPI_Request request;

  // Current PE is the sender, so check if receiver exists:
  if (iReceiver > -1) {
    iP = 0;

    for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      if (species[iSpecies].DoAdvect)
        DidWork = pack_border(species[iSpecies].density_scgc,
                              buffer, &iP, nG, iDir);
    }

    DidWork = pack_border(temperature_scgc, buffer, &iP, nG, iDir);

    for (int iComp = 0; iComp < 3; iComp++) {
      if (IsPole && iComp < 2)
        // Need to mirror zonal and meridional winds across pole:
        DidWork = pack_border(-velocity_vcgc[iComp], buffer, &iP, nG, iDir);
      else
        DidWork = pack_border(velocity_vcgc[iComp], buffer, &iP, nG, iDir);
    }
  }

  return DidWork;
}

// -----------------------------------------------------------------------------
// Unpack one face for the NEUTRALS (den, temp, vel)
// -----------------------------------------------------------------------------

bool Neutrals::unpack_one_face(int iSender,
                               precision_t *rbuffer,
                               int nG, int iDir,
                               bool DoReverseX,
                               bool DoReverseY,
                               bool XbecomesY) {

  bool DidWork = true;
  int64_t iP = 0;

  for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    if (species[iSpecies].DoAdvect)
      DidWork = unpack_border(species[iSpecies].density_scgc,
                              rbuffer, &iP, nG, iDir,
                              DoReverseX, DoReverseY, XbecomesY);
  }

  DidWork = unpack_border(temperature_scgc, rbuffer, &iP, nG, iDir,
                          DoReverseX, DoReverseY, XbecomesY);

  for (int iComp = 0; iComp < 3; iComp++) {
    DidWork = unpack_border(velocity_vcgc[iComp], rbuffer, &iP,
                            nG, iDir, DoReverseX, DoReverseY, XbecomesY);
  }

  return DidWork;
}


// -----------------------------------------------------------------------------
// Send for asynchronous communication. Don't touch buffer until mpi
// says that it has been received.
// -----------------------------------------------------------------------------

bool Grid::send_one_face(int64_t iFace) {

  bool DidWork = true;

  MPI_Isend(interchanges[iFace].buffer,
            interchanges[iFace].iSizeTotal,
            MPI_BYTE,
            interchanges[iFace].iProc_to,
            interchanges[iFace].iTag,
            aether_comm,
            &interchanges[iFace].requests);

  return DidWork;
}

// -----------------------------------------------------------------------------
// Receive asynchronously. Don't use the receive buffer until mpi says
// that it is ok.
// -----------------------------------------------------------------------------

bool Grid::receive_one_face(int64_t iFace) {

  bool DidWork = true;

  MPI_Recv(interchanges[iFace].rbuffer,
           interchanges[iFace].iSizeTotal,
           MPI_BYTE,
           interchanges[iFace].iProc_to,
           interchanges[iFace].iTag,
           aether_comm,
           MPI_STATUS_IGNORE);

  return DidWork;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------


Grid::messages_struct Grid::make_new_interconnection(int64_t iDir,
                                                     int64_t nVars,
                                                     int64_t iProc_to,
                                                     arma_vec edge_center,
                                                     bool IsPole,
                                                     bool DoReverseX,
                                                     bool DoReverseY,
                                                     bool XbecomesY) {

  messages_struct new_inter;

  int64_t nPtsX = nGCs * (nY - nGCs * 2) * (nZ - nGCs * 2);
  int64_t nPtsY = nGCs * (nX - nGCs * 2) * (nZ - nGCs * 2);

  new_inter.iFace = iDir;
  new_inter.DoReverseX = DoReverseX;
  new_inter.DoReverseY = DoReverseY;
  new_inter.IsPole = IsPole;
  new_inter.XbecomesY = XbecomesY;

  if (iDir == 0 || iDir == 2)
    new_inter.iSizeTotal = nVars * nPtsX * sizeof(precision_t);
  else
    new_inter.iSizeTotal = nVars * nPtsY * sizeof(precision_t);

  new_inter.buffer = static_cast<precision_t*>(malloc(new_inter.iSizeTotal));
  new_inter.rbuffer = static_cast<precision_t*>(malloc(new_inter.iSizeTotal));

  new_inter.iProc_to = iProc_to;

  new_inter.iTag =
    (edge_center(0) + 1) * 1000000 +
    (edge_center(1) + 1) * 10000 +
    (edge_center(2) + 1) * 100;

  return new_inter;
}

// -----------------------------------------------------------------------------
// Exchange messages for the NEUTRALS:
//   1. (first time) Set up the exchanging interfaces between edges / side
//      1a. Set up geometry weirdness from one edge to another for unpacking
//      1b. Set up across the pole weirdness for e/w velocities
//   2. Pack variables into buffers (all sides)
//   3. Send messages out (all sides asynchronously)
//   4. Receive messages (all sides asynchronously)
//   5. Wait for messages to be received
//   5. Unpack variables from all sides
// -----------------------------------------------------------------------------

bool Neutrals::exchange(Grid &grid) {

  std::string function = "Neutrals::exchange";
  static int iFunction = -1;
  report.enter(function, iFunction);

  bool DidWork = true;

  int iTag, iDir;
  int iSpecies;

  static int64_t iX, nX = grid.get_nX();
  static int64_t iY, nY = grid.get_nY();
  static int64_t iZ, nZ = grid.get_nZ();
  static int64_t nG = grid.get_nGCs();
  static int64_t nPtsX = nG * (nY - nG * 2) * (nZ - nG * 2);
  static int64_t nPtsY = nG * (nX - nG * 2) * (nZ - nG * 2);
  static bool IsFirstTime = true;

  static int64_t nVarsToPass = 0;

  if (IsFirstTime) {

    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      if (species[iSpecies].DoAdvect)
        nVarsToPass++;

    // Temperature + Velocities
    nVarsToPass += (1 + 3);

    bool ReverseY, ReverseX, IsPole, XbecomesY;

    // Message passing right:
    ReverseX = false;
    ReverseY = false;
    IsPole = false;
    XbecomesY = false;

    // This is for the CubeSphere Grid:
    if (grid.iRoot == 4 && grid.iRootXp == 2) {
      ReverseY = true;
      XbecomesY = true;
    }

    if (grid.iRoot == 5 && grid.iRootXp == 2) {
      ReverseX = true;
      XbecomesY = true;
    }

    grid.interchanges.push_back(grid.make_new_interconnection(0,
                                                              nVarsToPass,
                                                              grid.iProcXp,
                                                              grid.edge_Xp,
                                                              IsPole,
                                                              ReverseX,
                                                              ReverseY,
                                                              XbecomesY));

    // Message passing up:
    ReverseX = false;
    ReverseY = grid.DoesTouchNorthPole;
    IsPole = grid.DoesTouchNorthPole;
    XbecomesY = false;

    // This is for the CubeSphere Grid:
    if (grid.iRoot == 0 && grid.iRootYp == 5) {
      ReverseX = true;
      XbecomesY = true;
    }

    if (grid.iRoot == 2 && grid.iRootYp == 5) {
      ReverseY = true;
      XbecomesY = true;
    }

    if (grid.iRoot == 3 && grid.iRootYp == 5) {
      ReverseX = true;
      ReverseY = true;
    }

    if (grid.iRoot == 5 && grid.iRootYp == 3) {
      ReverseX = true;
      ReverseY = true;
    }

    grid.interchanges.push_back(grid.make_new_interconnection(1,
                                                              nVarsToPass,
                                                              grid.iProcYp,
                                                              grid.edge_Yp,
                                                              IsPole,
                                                              ReverseX,
                                                              ReverseY,
                                                              XbecomesY));

    // Message passing left:
    ReverseX = false;
    ReverseY = false;
    IsPole = false;
    XbecomesY = false;

    if (grid.iRoot == 5 && grid.iRootXm == 0) {
      ReverseY = true;
      XbecomesY = true;
    }

    if (grid.iRoot == 4 && grid.iRootXm == 0) {
      ReverseX = true;
      XbecomesY = true;
    }

    grid.interchanges.push_back(grid.make_new_interconnection(2,
                                                              nVarsToPass,
                                                              grid.iProcXm,
                                                              grid.edge_Xm,
                                                              IsPole,
                                                              ReverseX,
                                                              ReverseY,
                                                              XbecomesY));

    // Message passing down:
    ReverseX = false;
    ReverseY = grid.DoesTouchSouthPole;
    IsPole =  grid.DoesTouchSouthPole;
    XbecomesY = false;

    // This is for the CubeSphere Grid:
    if (grid.iRoot == 0 && grid.iRootYm == 4) {
      ReverseY = true;
      XbecomesY = true;
    }

    if (grid.iRoot == 2 && grid.iRootYm == 4) {
      ReverseX = true;
      XbecomesY = true;
    }

    if (grid.iRoot == 4 && grid.iRootYm == 3) {
      ReverseX = true;
      ReverseY = true;
    }

    if (grid.iRoot == 3 && grid.iRootYm == 4) {
      ReverseX = true;
      ReverseY = true;
    }

    grid.interchanges.push_back(grid.make_new_interconnection(3,
                                                              nVarsToPass,
                                                              grid.iProcYm,
                                                              grid.edge_Ym,
                                                              IsPole,
                                                              ReverseX,
                                                              ReverseY,
                                                              XbecomesY));

    IsFirstTime = false;
  }

  for (int iDir = 0; iDir < 4; iDir++) {
    if (report.test_verbose(2))
      std::cout << "packing : " << iDir << " " << iProc
                << " " << grid.interchanges[iDir].iProc_to
                << " " << grid.interchanges[iDir].iTag << "\n";

    DidWork = pack_one_face(grid.interchanges[iDir].iProc_to,
                            grid.interchanges[iDir].buffer,
                            nG, grid.interchanges[iDir].iFace,
                            grid.interchanges[iDir].IsPole);
  }

  // Send all faces asynchronously:
  for (int iDir = 0; iDir < 4; iDir++) {
    if (grid.interchanges[iDir].iProc_to >= 0)
      DidWork = grid.send_one_face(iDir);
  }

  // Receive all faces asynchronously:
  for (int iDir = 0; iDir < 4; iDir++) {
    if (grid.interchanges[iDir].iProc_to >= 0)
      DidWork = grid.receive_one_face(iDir);
  }

  // Wait for messages to get there:
  for (int iDir = 0; iDir < 4; iDir++) {
    if (grid.interchanges[iDir].iProc_to >= 0)
      MPI_Wait(&grid.interchanges[iDir].requests, MPI_STATUS_IGNORE);
  }

  // Unpack all faces:
  for (int iDir = 0; iDir < 4; iDir++) {
    if (grid.interchanges[iDir].iProc_to >= 0) {
      DidWork = unpack_one_face(grid.interchanges[iDir].iProc_to,
                                grid.interchanges[iDir].rbuffer,
                                nG, iDir,
                                grid.interchanges[iDir].DoReverseX,
                                grid.interchanges[iDir].DoReverseY,
                                grid.interchanges[iDir].XbecomesY);
    } else
      set_horizontal_bcs(iDir, grid);
  }

  // Wait for all processors to be done.
  MPI_Barrier(aether_comm);

  report.exit(function);

  return DidWork;
}

