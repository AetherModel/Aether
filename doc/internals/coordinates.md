
# Coordinate Systems

## GEO - Geographic coordinate system (or geodetic)
- Longitude (radians) - radians east of prime meridian (0 - 2pi)
- Latitude (radians) - angle between the equatorial plane and point. For oblate spheriod, this angle is not the same as the angle orthogonal to the surface of the planet and the equatorial plane.
- Radius (meters) - the distance to the center of the planet (Altitude is often used instead).

## PCPF - Planet-centered Planet-fixed or Geocentric coordinate system
- Cartesian coordinates of Geographic coordinate system in meters (X, Y, Z).
- X (meters) - aligned with the equator and prime meridian
- Z (meters) - aligned with rotation axis of the planet
- Y (meters) - completes the right-hand coordinate system

## PSE - Planetary Solar Ecliptic coordinates
- Cartesian coordinates tying together planet and Sun
- X (meters) - points from the center of the planet to the Sun
- Y (meters) - points from the center of the planet towards dusk and opposite planet's motion around the sun
- Z (meters) - orthogonal to the ecliptic plane

## Dipole Coordinates
- Longitude (radians) - radians east of the meridian that contains the north magnetic pole and north rotation axis
- P (meters) - Identifies the field line. This is the same as L-shell and is constant along wach field line.
- Q (dimensionless) - parameterizes the distance along the field line, related to magnetic latitude & radius. This varies along the field line, but the values are idential for all field lines within each node. q=0 at the equator, and approaches positive (negative) infinity as theta points towards the north (south) pole. Thus, q values will be negative in the southern hemisphere and the change in q "upwards" will be negative in the northern hemisphere. See [../../edu/examples/Dipole](../../edu/examples/Dipole) for more information.

## More Dipole Coordinates
- L-shell (Planetary Radii) - The distance from the planet's center at which the magnetic field encounters the dipole's equatorial plane
- Magnetic Latitude (radians) - angle between the dipole's equatorial plane and a point.
- Invariant Latitude (radians) - angle between the dipole's equatorial plane and the point at which the field-line passes through a reference radius of the planet ([specified in the inputs](../internals/grid.md#inputs)).  This is constant along the field-line and is related to the L-Shell.
- Magnetic Local Time (hours) - Angle between the sun, the north magnetic pole, and the point. Explicitly, this is done in PSE XY coordinates, ignoring the Z coorinate.

> The dipole `(i,j,k)` coordinates are (magnetic longitude, p, q).

# Coordinates in Aether

There are a variety of coordinates in Aether. This document describes some of them.

## Spherical Coordinates

The easiest coordinate system to understand within Aether is the spherical system, which is a longitude, latitude, radial (LLR) coordinate system. When the planet is a pure sphere, the LLR system is orthogonal - meaning that the grid lines up perfectly with the lines of constant longitude, latitude, and radius.

In Aether, longitude and latitude are expressed in radians and are positive towards the east and towards the north. Radius is expressed in meters and is positive away from the planet (upwards).  In Aether, often Altitude is used instead of radius. When the planet is a perfect sphere, these are offset by a constant value.

If the planet is an oblate spheriod, then the equator is larger than the pole, so that the planetary radius is dependent on latitude. Aether is currently set up so that Altitude is not dependent on latitude or longitude, so that if an oblate spheriod is used, then a constant altitude would have a radius that is dependent on latitude.  This means that the coordinate system is not purely orthogonal. At this time, this is not dealt with properly.

Because Aether considers gravity to be a function of radius and explicitly includes the centrifugal acceleration, the pertubation away from a perfect sphere should mostly cancel.

## i, j, k Coordinates

As described in the grid.md file, Aether uses a logical `(i, j, k)` 3D grid structure.  Therefore, we refer to the 'primary' coordinates as the ijk coordinate system.  What this means is that the i-coordinate is in the i-direction, the j-coordinate is in the j-direction, and the k-coordinate is in the k-direction.  

For the (perfectly) spherical grid, the i-coordinate is longitude, the j-coordinate is latitude, and the k-coordinate is radius.

For the Cubedsphere grid, the i-coordinate is RIGHT, the j-coodinate is UP, and the k-coordinate is radius.  Each face of the cubedsphere has the same coordinate system, but only with reference to that face.  This means that if each face is looked at independently, the lower left corner is at (about) i = -45, j = -45 deg, while the upper right corner is at i = +45, j = +45 deg. Radius is treated the same as in a spherical grid.

For the dipole coordinate system, the i-coordinate is magnetic longitude, the j-coordinate is L-shell, and the k-coordinate is Q: a dimensionless parameter, normalized to the planet radius, representing diatance along a magnetic field line. The dipole is orthogonal to a dipolar magnetic field.


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
- magnetic lon, invariant lat (only the dipole magnetic grid then has magnetic latitude)


