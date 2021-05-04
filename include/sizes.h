// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_SIZES_H_
#define INCLUDE_SIZES_H_

// This is the file that defines the number of grid points in each
// direction.  The entire code is based on these numbers, so you need
// to recompile if you change these numbers.
// 
// These are temporary and will eventually be removed.

// This is for the geographic grid:

#define nGeoGhosts 2

// This is for the magnetic grid:

#define nMagGhosts 2

#define nMagAlts 60
#define nMagAltsG nMagGhosts + nMagAlts + nMagGhosts
#define iMagAltStart_ nMagGhosts // Inclusive!!!
#define iMagAltEnd_ nMagGhosts + nMagAlts - 1 // Inclusive!!!

#define nMagLons 20
#define nMagLonsG nMagGhosts + nMagLons + nMagGhosts
#define iMagLonStart_ nMagGhosts // Inclusive!!!
#define iMagLonEnd_ nMagGhosts + nMagLons - 1 // Inclusive!!!

#define nMagLats 10
#define nMagLatsG nMagGhosts + nMagLats + nMagGhosts
#define iMagLatStart_ nMagGhosts // Inclusive!!!
#define iMagLatEnd_ nMagGhosts + nMagLats - 1 // Inclusive!!!

// This is for character string lengths:

#define nCharsShort 20
#define nCharsLong 400

#endif  // INCLUDE_SIZES_H_
