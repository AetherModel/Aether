// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "aether.h"

// -----------------------------------------------------------------------------
// This is where all of the exchange messages routines will sit.
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// Pack variables for message passing
//   value is variable to pack
//   packed is the buffer to pack it in
//   iCounter is the location in the buffer to start
//   nG is the number of ghost cells
//   iDir is the direction of the message pass:
//     0 - left -> right
//     1 - bottom -> top
//     2 - right -> left
//     3 - top -> bottom
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
//     0 - left -> right
//     1 - bottom -> top
//     2 - right -> left
//     3 - top -> bottom
//   IsPole indicates whether the message passing happened across the pole,
//     in which case the longitudes (Xs) have to be mirrored
// -----------------------------------------------------------------------------

bool unpack_border(arma_cube &value,
                   precision_t *packed,
                   int64_t *iCounter,
                   int64_t nG,
                   int iDir,
                   bool IsPole) {

  bool DidWork = true;
  static int64_t nX = value.n_rows;
  static int64_t nY = value.n_cols;
  static int64_t nZ = value.n_slices;

  int64_t iXstart, iXend;
  int64_t iYstart, iYend;
  int64_t xInc = 1;

  // ----------------------------
  // left / right message passing
  if (iDir == 0 || iDir == 2) {
    iYstart = nG;
    iYend = nY - nG;

    if (iDir == 0) {
      // left -> right
      iXstart = 0;
      iXend = nG;
    } else {
      // right -> left
      iXstart = nX - nG;
      iXend = nX - 1;
    }
  }

  // ----------------------------
  // top / bottom message passing
  if (iDir == 1 || iDir == 3) {
    if (!IsPole) {
      iXstart = nG;
      iXend = nX - nG;
    } else {
      // If we are at the pole, mirror the unpacking (left-to-right)
      iXend = nG;
      iXstart = nX - nG;
      xInc = -1;
    }

    if (iDir == 1) {
      // bottom -> top
      iYstart = 0;
      iYend = nG;
    } else {
      // top -> bottom
      iYstart = nY - nG;
      iYend = nY - 1;
    }
  }

  try {
    for (int64_t iZ = nG; iZ < nZ - nG; iZ++) {
      for (int64_t iY = iYstart; iY < iYend; iY++) {
        for (int64_t iX = iXstart; iX < iXend; iX += xInc) {
          value(iX, iY, iZ) = packed[*iCounter];
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
// Exchange one face for the NEUTRALS:
//   1. pack all of the variables (den, temp, vel)
//   2. send the buffer
//   3. receive the buffer
//   4. Unpack all of the variables (den, temp, vel)
//   5. Wait for everyone to finish (technically the send...)
// -----------------------------------------------------------------------------

bool Neutrals::exchange_one_face(int iReceiver, int iSender,
                                 precision_t *buffer,
                                 int64_t iTotalSize,
                                 int nG, int iDir) {

  bool DidWork = true;
  int64_t iP;
  MPI_Request request;
  int iTag = iDir;
  bool IsPole = false;

  if (iReceiver == iSender && (iDir == 1 || iDir == 3))
    IsPole = true;

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

    MPI_Isend(buffer, iTotalSize, MPI_BYTE,
              iReceiver, iTag, aether_comm, &request);
  }

  // Now we can switch to being the receiver, so check if sender exists:
  if (iSender > -1) {
    MPI_Recv(buffer, iTotalSize, MPI_BYTE, iSender, iTag, aether_comm,
             MPI_STATUS_IGNORE);

    iP = 0;

    for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      if (species[iSpecies].DoAdvect)
        DidWork = unpack_border(species[iSpecies].density_scgc,
                                buffer, &iP, nG, iDir, IsPole);
    }

    DidWork = unpack_border(temperature_scgc, buffer, &iP, nG, iDir, IsPole);

    for (int iComp = 0; iComp < 3; iComp++) {
      DidWork = unpack_border(velocity_vcgc[iComp], buffer, &iP,
                              nG, iDir, IsPole);
    }
  }

  if (iReceiver > -1)
    MPI_Wait(&request, MPI_STATUS_IGNORE);

  return DidWork;
}

// -----------------------------------------------------------------------------
// Exchange messages for the NEUTRALS:
//   1. Set up a bunch of variables and buffers
//   2. exchange messages in the east/west directions
//   3. exchange messages in the north/south directions (but ignore poles)
//   4. exchange messages across the poles
// -----------------------------------------------------------------------------

bool Neutrals::exchange(Grid grid, Report &report) {

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
  int64_t nPtsX = nG * (nY - nG * 2) * (nZ - nG * 2);
  int64_t nPtsY = nG * (nX - nG * 2) * (nZ - nG * 2);

  int64_t nVarsToPass = 0;

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    if (species[iSpecies].DoAdvect)
      nVarsToPass++;

  // Temperature + Velocities
  nVarsToPass += (1 + 3);

  int64_t iTotalSizeX = nVarsToPass * nPtsX * sizeof(precision_t);
  int64_t iTotalSizeY = nVarsToPass * nPtsY * sizeof(precision_t);

  // Create a temporary c-array to use to message pass the variables
  static precision_t *xBuffer = static_cast<precision_t*>(malloc(iTotalSizeX));
  static precision_t *yBuffer = static_cast<precision_t*>(malloc(iTotalSizeY));

  int64_t iP;

  // -------------------------
  // Pass from left -> right
  // -------------------------
  iDir = 0;
  DidWork = exchange_one_face(grid.iProcXp, grid.iProcXm, xBuffer,
                              iTotalSizeX, nG, iDir);

  // -------------------------
  // Pass from right -> left
  // -------------------------
  iDir = 2;
  DidWork = exchange_one_face(grid.iProcXm, grid.iProcXp, xBuffer,
                              iTotalSizeX, nG, iDir);

  // -------------------------
  // Pass from bottom -> top (not poles!)
  // -------------------------
  iDir = 1;
  int iUpPe = grid.iProcYp;
  int iDownPe = grid.iProcYm;

  if (grid.DoesTouchNorthPole)
    iUpPe = -1;

  if (grid.DoesTouchSouthPole)
    iDownPe = -1;

  DidWork = exchange_one_face(iUpPe, iDownPe, yBuffer, iTotalSizeY, nG, iDir);

  // -------------------------
  // Pass from top -> bottom (not poles!)
  // -------------------------
  iDir = 3;
  DidWork = exchange_one_face(iDownPe, iUpPe, yBuffer, iTotalSizeY, nG, iDir);

  // -------------------------
  // Pass Across North Poles - symmetric message passing
  // -------------------------
  if (grid.DoesTouchNorthPole) {
    iDir = 1;
    int iUpPe = grid.iProcYp;
    int iDownPe = grid.iProcYp;
    DidWork = exchange_one_face(iUpPe, iDownPe, yBuffer, iTotalSizeY, nG, iDir);
  }

  // -------------------------
  // Pass Across South Poles - symmetric message passing
  // -------------------------
  if (grid.DoesTouchSouthPole) {
    iDir = 3;
    int iUpPe = grid.iProcYm;
    int iDownPe = grid.iProcYm;
    DidWork = exchange_one_face(iUpPe, iDownPe, yBuffer, iTotalSizeY, nG, iDir);
  }

  report.exit(function);

  return DidWork;
}

