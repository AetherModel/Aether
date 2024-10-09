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
// This is the main exchange messages for the neutrals.
//   We are exchanging densities, temperatures, and velocities
// -----------------------------------------------------------------------------


bool Neutrals::exchange_old(Grid &grid) {

  std::string function = "Neutrals::exchange";
  static int iFunction = -1;
  report.enter(function, iFunction);

  bool DidWork = true;
  int64_t nGCs = grid.get_nGCs();

  for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    if (species[iSpecies].DoAdvect)
      DidWork = exchange_one_var(grid, species[iSpecies].density_scgc, false);
  }

  DidWork = exchange_one_var(grid, temperature_scgc, false);

  // velocity components:
  // reverse east across the pole:
  DidWork = exchange_one_var(grid, velocity_vcgc[0], true);
  // reverse north across the pole:
  DidWork = exchange_one_var(grid, velocity_vcgc[1], true);
  // don't reverse vertical across the pole:
  DidWork = exchange_one_var(grid, velocity_vcgc[2], false);

  report.exit(function);
  return DidWork;
}



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
// Pack one face for one variable - this is a generic code to pack one
// variable one one face.
// -----------------------------------------------------------------------------

bool pack_one_var_on_one_face(arma_cube var_scgc,
                              int iDirToPass,
                              Grid &grid) {

  static int nG = grid.get_nGCs();
  int iDir = grid.interchangesOneVar[iDirToPass].iFace;
  int iReceiver = grid.interchangesOneVar[iDirToPass].iProc_to;
  precision_t *buffer = grid.interchangesOneVar[iDirToPass].buffer;
  bool IsPole = grid.interchangesOneVar[iDirToPass].IsPole;

  bool DidWork = true;
  int64_t iP;

  // Current PE is the sender, so check if receiver exists:
  if (iReceiver > -1) {
    iP = 0;

    DidWork = pack_border(var_scgc, buffer, &iP, nG, iDir);

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
// Send for asynchronous communication. Don't touch buffer until mpi
// says that it has been received.
// -----------------------------------------------------------------------------

bool Grid::send_one_var_one_face(int64_t iFace) {

  bool DidWork = true;

  MPI_Isend(interchangesOneVar[iFace].buffer,
            interchangesOneVar[iFace].iSizeTotal,
            MPI_BYTE,
            interchangesOneVar[iFace].iProc_to,
            interchangesOneVar[iFace].iTag,
            aether_comm,
            &interchangesOneVar[iFace].requests);

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
// Receive asynchronously. Don't use the receive buffer until mpi says
// that it is ok.
// -----------------------------------------------------------------------------

bool Grid::receive_one_var_one_face(int64_t iFace) {

  bool DidWork = true;

  MPI_Recv(interchangesOneVar[iFace].rbuffer,
           interchangesOneVar[iFace].iSizeTotal,
           MPI_BYTE,
           interchangesOneVar[iFace].iProc_to,
           interchangesOneVar[iFace].iTag,
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

  if (iDir == 0 || iDir == 2) {
    new_inter.iSizeTotal = nVars * nPtsX * sizeof(precision_t);
    new_inter.index.set_size(nGCs, nY);
    new_inter.ratio.set_size(nGCs, nY);
  } else {
    new_inter.iSizeTotal = nVars * nPtsY * sizeof(precision_t);
    new_inter.index.set_size(nGCs, nX);
    new_inter.ratio.set_size(nGCs, nX);
  }

  new_inter.buffer = static_cast<precision_t*>(malloc(new_inter.iSizeTotal));
  new_inter.rbuffer = static_cast<precision_t*>(malloc(new_inter.iSizeTotal));

  new_inter.iProc_to = iProc_to;

  new_inter.iTag =
    (edge_center(0) + 1) * 1000000 +
    (edge_center(1) + 1) * 10000 +
    (edge_center(2) + 1) * 100;

  return new_inter;
}


/*
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

bool Neutrals::exchange_really_old(Grid &grid) {

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

  for (int i = 0; i < nSpecies; ++i)
    fill_corners(species[i].density_scgc, nG);
  
  fill_corners(temperature_scgc, nG);
  for (int iDir = 0; iDir < 3; iDir++)
    fill_corners(velocity_vcgc[iDir], nG);
  
  // Wait for all processors to be done.
  MPI_Barrier(aether_comm);

  report.exit(function);

  return DidWork;
}
*/


// -----------------------------------------------------------------------------
// Initialize interfaces between horizontal sides on a grid
// -----------------------------------------------------------------------------

bool exchange_sides_init(Grid &grid, int64_t nVarsToPass) {

  std::string function = "exchange_sides_init";
  static int iFunction = -1;
  report.enter(function, iFunction);

  bool DidWork = true;

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

  grid.interchangesOneVar.push_back(
    grid.make_new_interconnection(0,
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

  grid.interchangesOneVar.push_back(
    grid.make_new_interconnection(1,
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

  grid.interchangesOneVar.push_back(
    grid.make_new_interconnection(2,
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

  grid.interchangesOneVar.push_back(
    grid.make_new_interconnection(3,
                                  nVarsToPass,
                                  grid.iProcYm,
                                  grid.edge_Ym,
                                  IsPole,
                                  ReverseX,
                                  ReverseY,
                                  XbecomesY));

  report.exit(function);
  return DidWork;
}

// -----------------------------------------------------------------------------
// Exchange messages for one generic variable:
//   1. (first time) Set up the exchanging interfaces between edges / side
//      1a. Set up geometry weirdness from one edge to another for unpacking
//      1b. Set up across the pole weirdness for e/w velocities
//   2. Pack the variable into buffers (all sides)
//   3. Send messages out (all sides asynchronously)
//   4. Receive messages (all sides asynchronously)
//   5. Wait for messages to be received
//   5. Unpack variable from all sides
// -----------------------------------------------------------------------------

bool exchange_one_var(Grid &grid,
                      arma_cube &var_to_pass,
                      bool doReverseSignAcrossPole) {

  // This function is only needed if we do interpolation, which only happens in
  // the horizontal directions
  if (!grid.get_HasXdim() & !grid.get_HasYdim()) return true;

  std::string function = "exchange_one_var";
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
  static arma_cube var_scgc;

  static int64_t nVarsToPass = 1;

  if (IsFirstTime) {
    DidWork = exchange_sides_init(grid, nVarsToPass);
    var_scgc.set_size(nX, nY, nX);
    IsFirstTime = false;
  }

  int64_t iP;
  precision_t oneSign = 1.0;

  for (int iDir = 0; iDir < 4; iDir++) {
    if (report.test_verbose(2))
      std::cout << "packing one var : " << iDir << " " << iProc
                << " " << grid.interchangesOneVar[iDir].iProc_to
                << " " << grid.interchangesOneVar[iDir].iTag << "\n";

    if (grid.interchangesOneVar[iDir].IsPole &&
        doReverseSignAcrossPole)
      var_scgc = -1.0 * var_to_pass;
    else
      var_scgc = 1.0 * var_to_pass;

    // Current PE is the sender, so check if receiver exists:
    if (grid.interchangesOneVar[iDir].iProc_to > -1) {
      iP = 0;
      report.print(2, "Packing Border");
      DidWork = pack_border(var_scgc,
                            grid.interchangesOneVar[iDir].buffer,
                            &iP,
                            nG,
                            iDir);
      report.print(2, "Done Packing Border");
    }
  }

  // Send all faces asynchronously:
  for (int iDir = 0; iDir < 4; iDir++) {
    if (grid.interchangesOneVar[iDir].iProc_to >= 0) {
      report.print(2, "Sending one face");
      DidWork = grid.send_one_var_one_face(iDir);
    }
  }

  // Receive all faces asynchronously:
  for (int iDir = 0; iDir < 4; iDir++) {
    if (grid.interchangesOneVar[iDir].iProc_to >= 0) {
      report.print(2, "Receiving one face");
      DidWork = grid.receive_one_var_one_face(iDir);
    }
  }

  // Wait for messages to get there:
  for (int iDir = 0; iDir < 4; iDir++) {
    if (grid.interchangesOneVar[iDir].iProc_to >= 0)
      MPI_Wait(&grid.interchangesOneVar[iDir].requests, MPI_STATUS_IGNORE);
  }

  // Unpack all faces:
  for (int iDir = 0; iDir < 4; iDir++) {
    if (grid.interchangesOneVar[iDir].iProc_to >= 0) {
      iP = 0;
      report.print(2, "Unpacking Border");
      DidWork = unpack_border(var_to_pass,
                              grid.interchangesOneVar[iDir].rbuffer,
                              &iP,
                              nG,
                              iDir,
                              grid.interchangesOneVar[iDir].DoReverseX,
                              grid.interchangesOneVar[iDir].DoReverseY,
                              grid.interchangesOneVar[iDir].XbecomesY);
      report.print(2, "Done Unpacking Border");
    }
  }

  // Wait for all processors to be done.
  MPI_Barrier(aether_comm);

  // If this is a cubesphere grid, interpolate ghostcells to their proper location
  if (grid.IsCubeSphereGrid & grid.gcInterpolationSet) {
    report.print(3, "Interpolating Ghostcells to Proper Location");
    var_scgc = interpolate_ghostcells(var_to_pass, grid);
    var_to_pass = var_scgc;
  }

  // Now we fill in the corners so that we don't have zero values there:
  fill_corners(var_to_pass, nG);

  report.exit(function);
  return DidWork;
}

// -----------------------------------------------------------------------------
// This tests the ghostcell interpolation on the cubesphere by message passing
// latitudes and longitudes and checking to see if they match where they are
// expected to be.
// -----------------------------------------------------------------------------

bool test_ghostcell_interpolation(Grid &grid) {

  std::string function = "test_ghostcell_interpolation";
  static int iFunction = -1;
  report.enter(function, iFunction);

  bool didWork = true;

  static int64_t iX, nX = grid.get_nX();
  static int64_t iY, nY = grid.get_nY();
  static int64_t iZ, nZ = grid.get_nZ();
  static int64_t nG = grid.get_nGCs();
  int64_t iStart, iEnd, jStart, jEnd, iDir;

  // Check the latitudes and longitudes to make sure that they map to
  // the same location after message passing

  // original lats and lons in degrees:
  arma_cube lats = grid.geoLat_scgc * cRtoD;
  arma_cube lons = grid.geoLon_scgc * cRtoD;

  // save the lats and lons that are the gold standard:
  arma_cube latsGood = lats;
  arma_cube lonsGood = lons;

  // Set the ghostcells in the test cells to 0:
  set_gcs_to_value(lats, 0.0, nG);
  set_gcs_to_value(lons, 0.0, nG);

  // Message pass to get the lats and lons from the other faces
  report.print(1, "starting exchange to pass lats and lons");
  didWork = exchange_one_var(grid, lats, false);
  didWork = exchange_one_var(grid, lons, false);

  // Now we will go through can compare the cells from the message passed cells to
  // the saved cells:
  if (report.test_verbose(1)) {

    // Loop through the different directions:
    for (iDir = 0; iDir < 4; iDir++) {

      // Iteration indices - last point is <, so it is not included:
      // Right:
      if (iDir == 0) {
        iStart = nX - nG;
        iEnd = nX;
        jStart = nG;
        jEnd = nY - nG;
      }

      // Up:
      if (iDir == 1) {
        jStart = nY - nG;
        jEnd = nY;
        iStart = nG;
        iEnd = nX - nG;
      }

      // Left:
      if (iDir == 2) {
        iStart = 0;
        iEnd = nG;
        jStart = nG;
        jEnd = nY - nG;
      }

      // Down:
      if (iDir == 3) {
        jStart = 0;
        jEnd = nG;
        iStart = nG;
        iEnd = nX - nG;
      }

      // we really only care about one slice, so just take the first non-ghostcell slice:
      iZ = nG;
      std::cout << "Looping through iDir = " << iDir << "\n";

      for (iX = iStart; iX < iEnd; iX++) {
        for (iY = jStart; iY < jEnd; iY++) {

          std::cout << "iX, iY : " << iX << " " << iY << " ";
          std::cout << " lats message passed : " << lats(iX, iY, iZ);
          std::cout << " lats pre-pass : " << latsGood(iX, iY, iZ);

          if (!compare(latsGood(iX, iY, iZ), lats(iX, iY, iZ))) {
            std::cout << " <-- lats don't match!!! ";
            didWork = false;
          }

          std::cout << "\n";
          std::cout << "       lons message passed : " << lons(iX, iY, iZ);
          std::cout << " lons pre-pass : " << lonsGood(iX, iY, iZ);

          if (!compare(lonsGood(iX, iY, iZ), lons(iX, iY, iZ))) {
            std::cout << " <-- lons don't match!!! ";
            didWork = false;
          }

          std::cout << "\n";
        }
      }
    }
  }

  report.exit(function);
  return didWork;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

arma_cube interpolate_ghostcells(arma_cube varIn, Grid &grid) {

  std::string function = "interpolate_ghostcells";
  static int iFunction = -1;
  report.enter(function, iFunction);

  bool didWork = true;

  int64_t iDir;
  static int64_t iX, ix_, nX = grid.get_nX();
  static int64_t iY, iy_, nY = grid.get_nY();
  static int64_t iG, nG = grid.get_nGCs();
  precision_t r_;
  arma_cube varOut = varIn;

  // Loop through the different directions:
  for (iDir = 0; iDir < 4; iDir++) {

    // Need to do the interpolation for each row of ghost cells
    for (iG = 0; iG < nG; iG++) {

      // Right and Left sides of the block
      if (iDir == 0 || iDir == 2) {
        if (iDir == 0)
          iX = nX - 1 - iG;

        if (iDir == 2)
          iX = iG;

        for (iY = 0; iY < nY; iY++) {
          iy_ = grid.interchangesOneVar[iDir].index(iG, iY);
          r_ = grid.interchangesOneVar[iDir].ratio(iG, iY);

          if (iy_ > -1)
            varOut.tube(iX, iY) =
              (1.0 - r_) * varIn.tube(iX, iy_) + r_ * varIn.tube(iX, iy_ + 1);
        }
      }

      // Up and Down sides of the block
      if (iDir == 1 || iDir == 3) {
        if (iDir == 1)
          iY = nY - 1 - iG;

        if (iDir == 3)
          iY = iG;

        for (iX = 0; iX < nX; iX++) {
          ix_ = grid.interchangesOneVar[iDir].index(iG, iX);
          r_ = grid.interchangesOneVar[iDir].ratio(iG, iX);

          if (ix_ > -1)
            varOut.tube(iX, iY) =
              (1.0 - r_) * varIn.tube(ix_, iY) + r_ * varIn.tube(ix_ + 1, iY);
        }
      }
    }
  }

  report.exit(function);
  return varOut;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

bool find_ghostcell_interpolation_coefs(Grid &grid) {

  // This function is only needed if we do interpolation, which only happens in
  // the horizontal directions
  if (!grid.get_HasXdim() & !grid.get_HasYdim()) return true;

  std::string function = "find_ghostcell_interpolation_coefs";
  static int iFunction = -1;
  report.enter(function, iFunction);

  bool didWork = true;

  static int64_t iX, nX = grid.get_nX();
  static int64_t iY, nY = grid.get_nY();
  static int64_t iZ, nZ = grid.get_nZ();
  static int64_t nG = grid.get_nGCs();

  // Test to see if the longitudes are the same as the original
  arma_cube yOther = grid.refy_angle * cRtoD;
  arma_cube xOther = grid.refx_angle * cRtoD;
  arma_cube yLocal = grid.refy_angle * cRtoD;
  arma_cube xLocal = grid.refx_angle * cRtoD;

  iZ = nG;

  // Set the ghostcells to 0:
  set_gcs_to_value(yOther, 0.0, nG);
  set_gcs_to_value(xOther, 0.0, nG);

  // Message pass to get the X and Y from the other faces
  report.print(2, "starting exchange to find interpolation indices");
  didWork = exchange_one_var(grid, yOther, false);
  didWork = exchange_one_var(grid, xOther, false);

  if (report.test_verbose(2)) {

    std::cout << "finished exchange";

    std::cout << "xLocal : \n" << xLocal.slice(iZ) << "\n";
    std::cout << "xOther : \n" << xOther.slice(iZ) << "\n";
    std::cout << "yLocal : \n" << yLocal.slice(iZ) << "\n";
    std::cout << "yOther : \n" << yOther.slice(iZ) << "\n";
    std::cout << "xLocal, row : \n" << yOther.slice(iZ).row(2) << "\n";
    std::cout << "xLocal, col : \n" << yOther.slice(iZ).col(2) << "\n";
  }

  // Now that we have X, Y from other faces in their ghostcells and we
  // have our local X, Y in the ghost cells, we need to find
  // interpolation indices along each direction.

  // These are just variables to make things easier to deal with:
  arma_vec ind, rat;
  arma_vec from, to;

  // For the cubesphere, we can assume that we only need to do this
  // for one slice, since the angles are constant with altitude.

  precision_t r, y;
  int64_t iy_;
  int64_t iDir = 0;

  // Loop through the different directions:
  for (iDir = 0; iDir < 4; iDir++) {

    // Need to do the interpolation for each row of ghost cells
    for (int64_t iG = 0; iG < nG; iG++) {
      //std::cout << "iG : " << iG << " " << iX << "\n";

      if (iDir == 0)
        iX = nX - 1 - iG;

      if (iDir == 1)
        iY = nY - 1 - iG;

      if (iDir == 2)
        iX = iG;

      if (iDir == 3)
        iY = iG;

      // From are the message passed ghost cells (wider),
      // To are the locations of the native grid (narrower)
      if (iDir == 0 || iDir == 2) {
        // This is for right and left...
        if (grid.interchangesOneVar[iDir].XbecomesY) {
          // moving from a row to a vector, so do a transpose:
          from = xOther.slice(iZ).row(iX).t();
        } else
          from = yOther.slice(iZ).row(iX).t();

        if (grid.interchangesOneVar[iDir].DoReverseY)
          from = reverse(from);

        to = yLocal.slice(iZ).row(iX).t();
      } else {
        // This is for up and down...
        if (grid.interchangesOneVar[iDir].XbecomesY) {
          // moving from a row to a vector, so do a transpose:
          from = yOther.slice(iZ).col(iY);
        } else
          from = xOther.slice(iZ).col(iY);

        if (grid.interchangesOneVar[iDir].DoReverseX)
          from = reverse(from);

        to = xLocal.slice(iZ).col(iY);
      }

      if (report.test_verbose(3)) {
        std::cout << "iDir : " << iDir << "\n";
        std::cout << "drx : " << grid.interchangesOneVar[iDir].DoReverseX << "\n";
        std::cout << "dry : " << grid.interchangesOneVar[iDir].DoReverseY << "\n";
        std::cout << "xby : " << grid.interchangesOneVar[iDir].XbecomesY << "\n";
        std::cout << "from : " << from.t();
        std::cout << "to : " << to.t();
      }

      // This takes in vectors and returns them (ind and rat are vectors)
      didWork = find_interpolation_coefficients(from,
                                                to,
                                                ind,
                                                rat);
      grid.interchangesOneVar[iDir].index.row(iG) = ind.t();
      grid.interchangesOneVar[iDir].ratio.row(iG) = rat.t();

      if (report.test_verbose(3)) {
        std::cout << "from :\n";
        std::cout << from;
        std::cout << "to :\n";
        std::cout << to;
        std::cout << "ind :\n";
        std::cout << ind;
        std::cout << "rat :\n";
        std::cout << rat;
      }
    }
  }

  grid.gcInterpolationSet = true;

  if (report.test_verbose(3))
    didWork = test_ghostcell_interpolation(grid);

  // Wait for all processors to be done.
  MPI_Barrier(aether_comm);

  report.exit(function);
  return didWork;
}



