# Coordinates in Aether

There are a variety of coordinates in Aether. This document describes some of them.

## Spherical Coordinates

The easiest coordinate system to understand within Aether is the spherical system, which is a longitude, latitude, radial (LLR) coordinate system. When the planet is a pure sphere, the LLR system is orthogonal - meaning that the grid lines up perfectly with the lines of constant longitude, latitude, and radius.

In Aether, longitude and latitude are expressed in radians and are positive towards the east and towards the north. Radius is expressed in meters and is positive away from the planet (upwards).  In Aether, often Altitude is used instead of radius. When the planet is a perfect sphere, these are offset by a constant value.

If the planet is an oblate spheriod, then the equator is larger than the pole, so that the planetary radius is dependent on latitude. Aether is currently set up so that Altitude is not dependent on latitude or longitude, so that if an oblate spheriod is used, then a constant altitude would have a radius that is dependent on latitude.  This means that the coordinate system is not purely orthogonal. At this time, this is not dealt with properly.

Because Aether considers gravity to be a function of radius and explicitly includes the centrifugal acceleration, the pertubation away from a perfect sphere should mostly cancel.

## i, j, k Coordinates

As described in the grid.md file, Aether uses a logical '(i, j, k)' 3D grid structure.  Therefore, we refer to the 'primary' coordinates as the ijk coordinate system.  What this means is that the i-coordinate is in the i-direction, the j-coordinate is in the j-direction, and the k-coordinate is in the k-direction.  

For the (perfectly) spherical grid, the i-coordinate is longitude, the j-coordinate is latitude, and the k-coordinate is radius.

For the Cubedsphere grid, the i-coordinate is RIGHT, the j-coodinate is UP, and the k-coordinate is radius.  Each face of the cubedsphere has the same coordinate system, but only with reference to that face.  This means that if each face is looked at independently, the lower left corner is at (about) i = -45, j = -45 deg, while the upper right corner is at i = +45, j = +45 deg. Radius is treated the same as in a spherical grid.

Should the official coordinates be in the native coordinates (which could be different for each system), or should the coordinates be in meters, such that when gradients are taken, they are in /m?

Maybe we could have:

i_scgc, j_scgc, k_scgc - coordinates in the native coordinates (radians, meters, etc.)
im_scgc, jm_scgc, km_scgc - coordinate in meters

The question is what variables do we need?

Locations:
- Cell Centers (these are the center of each volume)
- Cell Edges in the i, j, k directions (these are the center of each area)
- Cell Corners

All locations should be described in the following coordinates:
- i, j, k
- lon, lat, radius (+alt)
- magnetic lon, invariant lat?

