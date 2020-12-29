// (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef AETHER_INCLUDE_SIZES_H_
#define AETHER_INCLUDE_SIZES_H_

// This is for the geographic grid:

#define nGeoGhosts 2

#define nGeoAlts 60
#define nGeoAltsG nGeoGhosts + nGeoAlts + nGeoGhosts
#define iGeoAltStart_ nGeoGhosts // Inclusive!!!
#define iGeoAltEnd_ nGeoGhosts + nGeoAlts - 1 // Inclusive!!!

#define nGeoLons 36
#define nGeoLonsG nGeoGhosts + nGeoLons + nGeoGhosts
#define iGeoLonStart_ nGeoGhosts // Inclusive!!!
#define iGeoLonEnd_ nGeoGhosts + nGeoLons - 1 // Inclusive!!!

#define nGeoLats 72
#define nGeoLatsG nGeoGhosts + nGeoLats + nGeoGhosts
#define iGeoLatStart_ nGeoGhosts // Inclusive!!!
#define iGeoLatEnd_ nGeoGhosts + nGeoLats - 1 // Inclusive!!!

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

#endif // AETHER_INCLUDE_SIZES_H_
