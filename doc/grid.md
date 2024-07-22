
# Grids in Aether

Aether uses a 2d domain decomposition and the grid system is basically a 3D i, j, k system, meaning that the arrays within Aether are 3D arrays.  Aether decomposes the grid in the first 2 dimensions (i and j) using a quadtree structure, while the 3rd dimension is left alone and each processor solves for the entire 3rd dimension.

Practically, what this means is that Aether uses powers of 4 to specify the grid system.  When you ask for 4x the number of processors it doubles the resolution in i and j. You can't double in i or j independently.

Aether using root nodes, which specify the smallest number of processors that can be run on.  For example, the simple "Sphere" grid has one root node that handles the entire Earth (in latitude and longitude). If the resolution needs to be doubles, 4 processors can be asked for.  If the resolution is doubled again, 16 processors are needed. Etc. However, in the altitude/radial direction, the number of points that are specified in the aether.json is never doubled when more processors are asked for. 

## Grid Types Explained

Aether has two types of grid systems - the neutral grid (NeuGrid) and the ion grid (IonGrid).  For each type of constituent (neutral or ion), their primary grid is the one where most of the equations are solved, and then they are passed to the other grid.  For example, the neutral winds are solved for on the NeuGrid, and then passed onto the IonGrid in order to calculate source terms for the ions.  As another example, the ion advection is solved for on the ion grid.  The ion densities are then passed to the neutral grid, where the source terms for the neutrals are calculated.

These grids can be identical or nearly identical. If they are, then it is best to have them on a neutral type of grid, since the stability of the neutrals along the 3rd dimension (where gravity is prime) is hard to achieve.

The neutral grid system typically has its third axis aligned (mostly) with the radial direction. This is to allow special solvers to treat gravity and the gradient in pressure in a special way. There are two issues with solving the neutrals in the third dimension: (1) often, the top of thermosphere models are supposed to be the exosphere, which means that they can't extend too far in the vertical direction; and (2) a problem with neutral solvers that can solve the full momentum equation in the third dimension is that they struggle with having too many scale heights in a domain. These limit the full height of the model domain for the neutrals.

For the ions, with systems that have magnetic fields, the plasma often moves up fieldlines on the dayside and down fieldlines on the nightside.  This transport is often above the top of the neutral grid. Further, the ions are often structured by the magnetic field, making this the natural coordinate system.  For planets without magnetic fields, a grid similar to the neutrals may be useful.  The ion grid can extend above and below the neutral grid with both a magnetic-field-aligned grid and a spherical grid.

## Grid Shapes Explained

Aether currently has three basic grid shapes: spherical, cubesphere, and dipole. The spherical grid is a (i, j, k) = (longitude, latitude, altitude) system, with these being (mostly) orthogonal to each other. This grid system can simulate a sub-region of the Earth if desired. 

### TL;DR

The user needs to specify the shape of the grid, which specifies the grid shape and the number of root nodes. Shapes include: sphere (1 root node), sphere6 (6 root nodes), cubesphere (6 root nodes), dipole (1 root node), dipole4 (4 root nodes), and dipole6 (6 root nodes).

### The Sphere Grid

The sphere grid is a normal longitude, latitude, altitude grid.

### The Cubesphere Grid

The cubesphere grid is composed of 6 different faces, similar to a cube, but where each cube "face" is pushed out to form a sphere.  The corners of the cube intersect the sphere, and all of the other grid points on the cube are pushed out to form the cube.  One cube face defines the southern polar region, one face defines the northern polar region, and the other four faces are spaced in longitude around the equatorial region. For the cubesphere grid, the (k) dimension is altitude.  The (i, j) system is set up so that (i) is considered left-right on the cube face, while (j) is up-down. For the four faces around the equator, (i, j) is roughly (lon, lat), but not exactly.  For the polar faces, the relationship between (i, j) and longitude, latitude is much more complex.

For both the spherical grid and the cubesphere grid, the altitudes (k) can be stretched or uniform. A lower boundary is set and the delta-altitude is specified as either a constant distance or a constant percentage of the bulk scale-height.

### The Dipole Grid

The dipole grid is aligned with the magnetic field. The (k) dimension is along the fieldline, (i) is magnetic longitude, and (j) is roughly latitude for the bottom of the field-line.  Each field-line starts at the lowest altitude and curves towards the equator.  In the northern hemisphere, this means that the field-lines curve south, while in the southern hemisphere, they curve north. The latitudinal spacing is such that there is a dependence on the L-shell (i.e., the equatorial radial extent of the field-line).  Along the (k) dimension, field-lines either terminate when they reach the equatorial plane, forming half of a full field-line or they terminate at the highest point specified in the aether.json file.  Any grid point with an L-shell less than the peak altitude will terminate in the equatorial plane, while any field-line that has an L-shell above the peak altitude will simply terminate.  Field-lines that terminate in the equatorial plane have corresponding field-lines in the other hemisphere, so ghostcells are used to pass information back and forth.  Field-lines that terminate at the maximum altitude have vertical boundary conditions set in the ghost cells.

The transition from the "closed" field-line region to the "open" field-line region is a natural break point in the grid. The transition between these regions can be handled with ghostcells in the "latitudinal" direction.  Therefore it makes the most sense to have have 4 distinct regions in "latitude": south open, south closed, north closed, north open.  The message passing is treated differently at the boundaries between each of these regions.  

### Root Nodes

A fundamental assumption within Aether is that each processor does computation on one and only one block. This means that each processor does not deal with multiple blocks, and therefore the distribution of blocks across processors has to match exactly.  This document uses the words "block" and "node" somewhat interchangably.  Technically, a "block" is single (i, j, k) grid, while a "node" can be multiple "blocks" that make up a section of the globe.

Aether uses a quadtree system to subdivide and distribute the grids (or blocks) across processors.  This means that when an additional level of refinement is desired, an individual block is split in 4 - the number of blocks is doubled in both (i) and (j). The question then is - how many blocks to start with?  These are the root nodes.

For the whole globe sphere shape, there is one single root node, which allows users to run the code on a single processor.  When a user asks for one processor using the sphere shape, there is only one single block, which is the root node, which spans then entire globe. When the user asks for four processors using the sphere shape, the number of blocks in latitude are doubled and the number of blocks in longitude are doubled. There is still only one root node, but the number of blocks is four, with 2 in the longitudinal direction and 2 in the latitudinal direction.  If the user asks for 16 processors using the sphere shape, the blocks are sub-divided again, with sill one single root node, and 16 blocks - four in the longitudinal direction and four in the altitudinal direction.  With a sphere grid, the number of processors that can be used to specify the grid are then: 1, 4, 16 (=4^2), 64 (=4^3), 256 (=4^4), 1024 (=4^5), etc.

With a cubesphere grid, there are six root nodes, meaning that the code needs six processors to run on just to start. Each root node is a face of the cubesphere. If the user asks for 24 processors (i.e., 6 root nodes that are each divided into four blocks each), each root node is split in half along the left-right direction and the up-down direction. For a cubesphere grid, the number of processors that can be used to specify the grid are then: 6, 24 (= 6 * 4), 96 (6 * 4^2), 384 (6 * 4^3), etc.

The root nodes indicate the span of the grid that they cover. This is done in a header file.  In order to accomplish this, the lower-left corner location (ORGINS) is specified as well the span of the root node in the left-to-right (i) direction (RIGHTS) and in the (j) down-to-up direction (UPS). The easiest example is here:
```bash
namespace Sphere {
  /// The normalized origins of each face of the cube (i.e. corner)
  static const arma_mat ORIGINS = {
				   { 0.0, -0.5, 0.0}
  };
  /// Normalized right steps in cube
  static const arma_mat RIGHTS = {
				  {2.0, 0.0, 0.0}
  };
  /// Normalized right steps in cube
  static const arma_mat UPS = {
			       {0.0, 1.0, 0.0}
  };
};
```

Since the sphere goes from -90 deg to +90 deg in latitude, and 0 deg to 360 deg in longitude, and pi is the normalizer, then the grid should go from -0.5 to 0.5 in the UPS direction, so the ORIGIN is placed at -0.5 and the span is 1.0.  In the longitudinal direction, the grid should go from 0 - 2, so the ORIGIN is placed at 0.0 and the span is 2.0.

This could be altered to have two root nodes.  If someone wanted the root node to be "square" in that the lat and lon spans are the same, this could be done:
```bash
namespace Sphere2 {
  /// The normalized origins of each face of the cube (i.e. corner)
  static const arma_mat ORIGINS = {
				   { 0.0, -0.5, 0.0},
				   { 1.0, -0.5, 0.0}
  };
  /// Normalized right steps in cube
  static const arma_mat RIGHTS = {
				  {1.0, 0.0, 0.0},
				  {1.0, 0.0, 0.0}
  };
  /// Normalized right steps in cube
  static const arma_mat UPS = {
			       {0.0, 1.0, 0.0},
			       {0.0, 1.0, 0.0}
  };
};
```
Notice that the namespace is different, so that it can be unique. In this case, there are two ORGINS (offset by 1.0 in longitude), two RIGHTS (which are the same), and two UPS (which are the same).

In both of these examples, the third dimension doesn't change.  This is because a single altitude in a spherical grid can be fully described with two variables (lat and lon).  For a cubesphere grid, on the other hand, there are three variables that are needed - each face is in an XY, XZ, or YZ plane, so all three (X, Y, and Z) are needed. The 6 root nodes for the cubesphere are specified in the cubesphere.h header file.

## Specifying the Grid

There are many different components to specifying the actual grid that is desired, namely:
- Min and Max latitude
- Min and Max longitude
- Min Altitude, whether a stretched altitude is desired, and the altitudinal spacing

In addition, the number of grid points that should be used in each block are specified:
- nLons or nX - number of grid cells per block in the i direction
- nLats or nY - number of grid cells per block in the j direction
- nAlts or nZ - total number of grid cells in the block in the k direction

### Horizontal Resolution

For some grid shapes (Sphere and Dipole), the total number of grid cells in the i and j direction can be determined by, for example, multiplying the number of blocks in the i direction by the number of cells in each block in the i direction. So, with a spherical grid with one root node, and 256 processors used, the number of blocks in the i and j direction is (256 = 4 * 4 * 4 * 4. breaking it into both directions - (2*2) * (2*2) * (2*2) * (2*2) or (2 * 2 * 2 * 2) * (2 * 2 * 2 * 2) or 16 x 16) 16 and 16. So, the total number of blocks in the i direction is 16 * nLons and in the j direction is 16 * nLats.

For the cubesphere grid, the nX and nY are the number of grid cells in the i and j direction.  At this time, these have to be identical in order to have the grid cells match up along the boundaries to the top and bottom nodes. The resolution of the Cubesphere grid is roughly 360 deg / (4 * nX * sqrt(nProc/6)).  For example, if nX = 18, and 24 processors are requested, then the resolution = 360 / (4 * 18 * sqrt(4)) = 360 / (72 * 2) = 2.5 deg. As another example, to make a grid with 1 deg resolution, with 96 processors, nX would have to be 1 = 360 / (4 * nX * 4) = 22.5 / nX, so nX has to be around 22. (If nX were 22, and nProc = 96, then the resolution would be 1.02 deg).

### Vertical Resolution

In all grids, the nAlts or nZ are not parallelized, so the number of points in the k direction is what is specified. For the sphere and cubesphere grids, this is the number of altitude points.  On the dipole grid, this is the number of points along the dipole flux tube.




```bash
    "NeuGrid" : {
        "Shape" : 
	    "MinLat" : -90.0,
	    "MaxLat" :  90.0,
	    "MinLon" :   0.0,
	    "MaxLon" : 360.0,
	    "MinAlt" : 100.0,
	    "dAlt"   : 5.0,
        "AltFile" : "",
	    "IsUniformAlt" : true},
```