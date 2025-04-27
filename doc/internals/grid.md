# Grids in Aether

Aether uses a 2d domain decomposition and the grid system is basically a 3D `(i,
j, k)` system, meaning that the arrays within Aether are 3D arrays. Aether
decomposes the grid in the first 2 dimensions (`i` and `j`) using a quadtree
structure, while the 3rd dimension is left alone and each processor solves for
the entire 3rd dimension.

Practically, what this means is that Aether uses powers of 4 to specify the grid
system. When you ask for 4x the number of processors it doubles the resolution
in `i` and `j`. You can't double in `i` or `j` independently.

Aether using root nodes, which specify the smallest number of processors that
can be run on. For example, the simple "Sphere" grid has one root node that
handles the entire Earth (in latitude and longitude). If the resolution needs
to be doubles, 4 processors can be asked for. If the resolution is doubled
again, 16 processors are needed, etc. However, in the altitude/radial direction,
the number of points that are specified in the aether.json is unchanged, as it
does not rely on the number of processors used.

- [Grid Types Explained](#grid-types-explained)
- [Grid Shapes Explained](#grid-shapes-explained)
  - [TL;DR](#tldr)
  - [The Sphere Grid](#the-sphere-grid)
  - [The Cubesphere Grid](#the-cubesphere-grid)
  - [The Dipole Grid](#the-dipole-grid)
    - [Inputs:](#inputs)
  - [Root Nodes](#root-nodes)
    - [Sphere](#sphere)
    - [Cubesphere](#cubesphere)
    - [Dipole](#dipole)
    - [Specifying Root Nodes](#specifying-root-nodes)
- [Specifying the Grid](#specifying-the-grid)
  - [Horizontal Resolution](#horizontal-resolution)
  - [Vertical Resolution](#vertical-resolution)

## Grid Types Explained

Aether has two types of grid systems - the neutral grid (`neuGrid`) and the ion
grid (`ionGrid`). For each type of constituent (neutral or ion), their primary
grid is the one where most of the equations are solved, and then they are passed
to the other grid. For example, the neutral winds are solved for on the
`neuGrid`, and then passed onto the `ionGrid` in order to calculate source terms
for the ions. As another example, the ion advection is solved for on the
`ionGrid`. The ion densities are then passed to the `neuGrid`, where the source
terms for the neutrals are calculated.

These grids can be identical or nearly identical. If they are, then it is best
to have them on a neutral type of grid, since the stability of the neutrals
along the 3rd dimension (where gravity is prime) is hard to achieve.

The neutral grid system typically has its third axis aligned (mostly) with the
radial direction. This is to allow special solvers to treat gravity and the
gradient in pressure in a special way. There are two issues with solving the
neutrals in the third dimension: (1) often, the top of thermosphere models are
supposed to be the exosphere, which means that they can't extend too far in the
vertical direction; and (2) neutral solvers struggle with having too many scale
heights in a domain when solving the full momentum equation. These limit the
full height of the model domain for the neutrals.

For the ions, with systems that have magnetic fields, the plasma often moves up
field-lines on the dayside and down field-lines on the nightside. This transport
is often above the top of the neutral grid. Further, the ions are often
structured by the magnetic field, making this the natural coordinate system. For
planets without magnetic fields, a grid similar to the neutrals may be useful.
The ion grid can extend above and below the neutral grid with both a
magnetic-field-aligned grid and a spherical grid.

## Grid Shapes Explained

Aether currently has three basic grid shapes: `spherical`, `cubesphere`, and
`dipole`. The spherical grid is an (`i`, `j`, `k`) = (longitude, latitude,
altitude) system, with these being (mostly) orthogonal to each other. This grid
system can simulate a sub-region of the Earth if desired.

### TL;DR

The user needs to specify the shape of the grid, which specifies the grid shape
and the number of root nodes. Shapes include: `sphere` (1 root node), `sphere6`
(6 root nodes), `cubesphere` (6 root nodes),`dipole4` (4 root nodes), and
`dipole6` (6 root nodes).

### The Sphere Grid

The sphere grid is a normal longitude, latitude, altitude grid.

### The Cubesphere Grid

The cubesphere grid is composed of 6 different faces, similar to a cube, but
where each cube "face" is pushed out to form a sphere. The corners of the cube
are first set to intersect the sphere, then all of the other grid points on the
cube are pushed out until they intersect the sphere.

| <a title="A2569875, CC BY-SA 4.0 &lt;https://creativecommons.org/licenses/by-sa/4.0&gt;, via Wikimedia Commons" href="https://commons.wikimedia.org/wiki/File:Cube_with_spherical_cube.gif"><img width="512" alt="Cube with spherical cube" src="https://upload.wikimedia.org/wikipedia/commons/f/f6/Cube_with_spherical_cube.gif?20200830172901"></a> |
|:--:|
| A graphical representation of the cubesphere grid. <br> (*Source: A2569875, CC BY-SA 4.0 <https://creativecommons.org/licenses/by-sa/4.0>, via Wikimedia Commons*) |

One cube face defines the southern polar region, one face defines the northern
polar region, and the other four faces are spaced in longitude around the
equatorial region. For the cubesphere grid, the `k` dimension is altitude. The
(`i, j`) system is set up so that `i` is considered left-right on the cube face,
while `j` is up-down. For the four faces around the equator, (`i, j`) is roughly
(longitude, latitude), but not exactly. For the polar faces, the relationship
between (`i, j`) and (longitude, latitude) is much more complex.

For both the spherical grid and the cubesphere grid, the altitudes, `k`, can be
stretched or uniform. A lower boundary is set and the delta-altitude is
specified as either a constant distance or a constant percentage of the bulk
scale-height.

### The Dipole Grid

The dipole grid is aligned with the magnetic field. The `k` dimension is along
the fieldline, `i` is magnetic longitude, and `j` is roughly latitude for the
bottom of the field-line. 

Each fieldline starts at the lowest modeled altitude
and curves towards the equator. In the northern hemisphere, this means that the
fieldlines curve south, while in the southern hemisphere they curve north.

The dipole grid is evenly spaced in **invariant latitude** (where the field line
passes the minumum altitude) and **q** (the dipole coordinate
specifying how far along the field line a point lies). Q is dimensionless and defined 
to be $-\infty$ at the south pole, $+\infty$ at the north pole, and 0 at the
magnetic equator. The equations for p (L-shell) and q are the following,
where r is the distance from the origin and $\theta$ is *colatitude*:

```math
p = \frac{r}{\sin^2\theta}
```

```math
q = \frac{\cos{\theta}}{r^2}
```

Here is how the dipole grid is generated:

1. Receive latitude range of this block from the quadtree. This will look 
something like `lower_left_norm=(0.0, -0.5, 0.0)` and `size_up_norm=(0, 0.25, 0)`
for the node nearest the south pole in dipole4. From this, determine if we are
in the southern hemisphere. If we are, everything will be done as if it was the north
hemisphere and then reversed & negated at the end.
2. Store the latitude (j) component of `lower_left_norm` as `lat_origin`. If this
node is in the southern hemisphere, store the top of the node's extent as lat_origin.
3. Scale this node's portion of the quadtree to be limited by the user-provided
`lat_range`. These for the invariant latitudes, which are evenly spaced between 
the latitude range provided and dictate where each field line passes through the minimum
altitude provided.
   - At the poles, put the last corner at $89.9^\circ$ magnetic latitude, or
$0.1^\circ$ and $179.9^\circ$ magnetic ***co***latitude. Add
another corner 1/2 way between this point and the last "real" corner, and put
cell centers between these corners.
4. Determine if this node will have closed or open field lines. There are two conditions:
   - If the node is touching the equator
   - If the lowest L-shell is below the maximum altitude. This is rare, but prevents unexpected behavior.
5. Determine the limits, then values, of the q-coordinate for all points along each field
lines on this node. The q-values on each node are identical, and the p-value is
constant along each field line (by definition). To solve for q, use the p-values
from step 3 and the altitude, as described below and $q=\sqrt{(1-r/p)/r^4}$.
   - If the field line closes, $q_{min}=0$. There will be a corner/edge at the 
magnetic equator and two ghost cell centers across the equator for message passing.
   - If the field line does not close, $q_{min}$ is calculated from the highest
altitude point on the lowest latitude field line. This is the point farthest
from the planet on the most equatorward field line (and since q=0 at the
equator, it has the lowest allowed q-value).
   - The maximum q-value is solved for identically in open & closed blocks with
the lower altitude limit and the highest latitude field line. The point closest
to the planet on the highest latitude field line has the highest allowed q-value
(q=$\pm$infinity at the poles).
6. We now have `p` (step 3) and `q` (step 5) for all points on the grid. From this
we solve for $(r, \theta)$, and any other coordinates we need.

See [edu/examples/Dipole](../../edu/examples/Dipole) for more detailed information
and to experiment with the available options in a Python script.

#### Inputs:

- ***Shape***: either `dipole4` or `dipole6`. Cannot be run on a single core.
- ***nLonsPerBlock***: number of magnetic longitudes
- ***nLatsPerBlock***: number of field lines (invariant latitudes)
- ***nAlts***: Number of points along each field line. A number of these will 
be discarded for being at too low of altitude.
- ***AltRange***: (`min_alt`, `max_alt`) - the altitude (in km) range to bound 
cells by. 
- ***LatRange***: (`min_lat`, `max_lat`) - the limits on invariant latitudes 
(in degrees). Sets the limits on the latitudes where field lines cross `min_alt`.


### Root Nodes

> This document uses the words "block" and "node" somewhat
interchangably. Technically, a "block" is single (`i, j, k`) grid, while a
"node" can be multiple "blocks" that make up a section of the globe.

A fundamental assumption within Aether is that each processor does computation
on one and only one block. This means that each processor does not deal with
multiple blocks, and therefore the distribution of blocks across processors has
to match exactly.

Aether uses a quadtree system to subdivide and distribute the grids (or blocks)
across processors. This means that when an additional level of refinement is
desired, an individual block is split in 4 - the number of blocks is doubled in
both `i` and `j`. The question then is *how many blocks to start with*? These
are the root nodes.

#### Sphere

For the whole globe `sphere` shape, there is one single root node, which allows
users to run the code on a single processor. When a user asks for one processor
using this sphere shape, there is only one single block, which is the root node,
and spans the entire globe. When the user asks for four processors using the
sphere shape, the number of blocks in latitude are doubled and the number of
blocks in longitude are doubled. There is still only one root node, but the
number of blocks is four, with 2 in the longitudinal direction and 2 in the
latitudinal direction. If the user asks for 16 processors using the sphere
shape, the blocks are sub-divided again, with sill one single root node, and 16
blocks - four in the longitudinal direction and four in the altitudinal
direction. With a sphere grid, the number of processors that can be used to
specify the grid are then: 1, 4, 16 (=4^2), 64 (=4^3), 256 (=4^4), 1024 (=4^5),
etc.

#### Cubesphere

With a cubesphere grid, there are six root nodes, meaning that the code needs
six processors to run on just to start. Each root node is a face of the
cubesphere. If the user asks for 24 processors (i.e., 6 root nodes that are each
divided into four blocks each), each root node is split in half along the
left-right direction and the up-down direction. For a cubesphere grid, the
number of processors that can be used to specify the grid are then: 6, 24 (6
\* 4), 96 (6 \* 4^2), 384 (6 \* 4^3), etc.

#### Dipole

The dipole grid requires >4 root nodes to ensure the coordinates are 
mutually orthogonal. The available shapes are `dipole4` and `dipole6`, for 
compatibility with the neutral grid being a sphere or cubesphere. In both cases,
each root node covers the entire longitude range and given a portion of the latitude 
range. So in the case of `dipole4`, the four nodes each cover 1/4 of the available
latitudes and all of the longitudes. The available latitudes are scaled to the latitude
limits specified in the input file, so the divisions will not be at $\pm45^\circ$ and
$0^\circ$ latitude, rather will be offset to evenly divide the entire range across the
blocks. Dividing the root nodes works identically to the spherical grid, for example
`dipole4` can be used with 16 MPI tasks and each root node is divided into four blocks,
forming a 2x2 grid.

#### Specifying Root Nodes

The root nodes indicate the span of the grid that they cover. This is done in a
header file. In order to accomplish this, the lower-left corner location
(ORGINS) is specified as well the span of the root node in the left-to-right (i)
direction (RIGHTS) and in the (j) down-to-up direction (UPS). The easiest
example is here:

```cpp
namespace Sphere {
 /// The normalized origins of each face of the cube (i.e. corner)
 static const arma_mat ORIGINS = {
                 {0.0, -0.5, 0.0}
 };
 /// Normalized right steps in cube
 static const arma_mat RIGHTS = {
                 {2.0,  0.0, 0.0}
 };
 /// Normalized right steps in cube
 static const arma_mat UPS = {
                 {0.0,  1.0, 0.0}
 };
};
```

Since the sphere goes from -90$^\circ$ to +90$^\circ$ in latitude, and 0$^\circ$
to 360$^\circ$ in longitude, and pi is the normalizer, then the grid should go
from -0.5 to 0.5 in the UPS direction, so the ORIGIN is placed at -0.5 and the
span is 1.0. In the longitudinal direction, the grid should go from 0 to 2, so
the ORIGIN is placed at 0.0 and the span is 2.0.

This could be altered to have two root nodes. If someone wanted the root node to
be "square", in that the latitude and longitude spans are the same, this could
be done with:

```cpp
namespace Sphere2 {
 /// The normalized origins of each face of the cube (i.e. corner)
 static const arma_mat ORIGINS = {
                 {0.0, -0.5, 0.0},
                 {1.0, -0.5, 0.0}
 };
 /// Normalized right steps in cube
 static const arma_mat RIGHTS = {
                 {1.0,  0.0, 0.0},
                 {1.0,  0.0, 0.0}
 };
 /// Normalized right steps in cube
 static const arma_mat UPS = {
                 {0.0,  1.0, 0.0},
                 {0.0,  1.0, 0.0}
 };
};
```

Notice that the namespace is different, so that it can be unique. In this case,
there are two ORGINS (offset by 1.0 in longitude), two RIGHTS (which are the
same), and two UPS (which are the same).

In both of these examples, the third dimension doesn't change. This is because a
single altitude in a spherical grid can be fully described with two variables
(lat and lon). For a cubesphere grid, on the other hand, there are three
variables that are needed - each face is in an XY, XZ, or YZ plane, so all three
(X, Y, and Z) are needed. The 6 root nodes for the cubesphere are specified in
the [cubesphere.h header file](../../include/cubesphere.h).

## Specifying the Grid

There are many different components to specifying the actual grid that is
desired, namely:

- Min and Max latitude
- Min and Max longitude
- Min Altitude, whether a stretched altitude is desired, and the altitudinal
  spacing

In addition, the number of grid points that should be used in each block are
specified:

- nLons or nX - number of grid cells per block in the i direction
- nLats or nY - number of grid cells per block in the j direction
- nAlts or nZ - total number of grid cells in the block in the k direction

### Horizontal Resolution

For some grid shapes (`sphere` and `dipole`), the total number of grid cells in
the `i` and `j` direction can be determined by, for example, multiplying the
*number of blocks* in the `i` direction by the *number of cells in each block*
in the `i` direction. So, with a spherical grid with one root node, and 256
processors used, the number of blocks in the `i` and `j` direction is 256 ( = 4
\* 4 \* 4 \* 4). Breaking it into both directions - (2\*2) \* (2*2) \* (2\*2) \*
(2\*2) or (2 \* 2 \* 2 \* 2) \* (2 \* 2 \* 2 \* 2) or (16 \* 16) - 16 and 16.
So, the total number of blocks in the `i` direction is 16 \* nLons and in the
`j` direction is 16 \* nLats.

For the cubesphere grid, the `nX` and `nY` are the number of grid cells in the
`i` and `j` direction. At this time, these have to be identical in order to have
the grid cells match up along the boundaries between the top and bottom nodes.
The resolution of the Cubesphere grid is roughly (360$^\circ$ / (4 \* nX \*
sqrt(nProc/6))). For example, if `nX` = 18, and 24 processors are requested,
then the resolution = 360 / (4 \* 18 \* sqrt(4)) = 360 / (72 \* 2) =
2.5$^\circ$. As another example, to make a grid with 1$^\circ$ resolution with
96 processors, we get (1 = 360 / (4 \* nX \* 4) = 22.5 / nX), so `nX` has to be
around 22. (If `nX` were 22, and `nProc` = 96, then the resolution would be
1.02$^\circ$).

### Vertical Resolution

In all grids, the nAlts (`nZ`) are not parallelized, so the number of points
in the `k` direction is what is specified by the user. For the `sphere` and
`cubesphere` grids, this is the number of altitude points. On the `dipole` grid,
this is the number of points along the dipole flux tube.

```json
  "neuGrid" : {
    "Shape" : "sphere",
    "LatRange" : [-90.0, 90.0],
    "nLatsPerBlock" : 18,
    "LonRange" : [0.0, 360.0],
    "nLonsPerBlock" : 36,
    "nAlts" : 50,
    "MinAlt" : 100.0,
    "dAltkm" : 5.0,
    "dAltScale" : 0.25,
    "IsUniformAlt" : true,
    "AltFile" : ""},
```

```json
    "ionGrid": {
        "Shape": "dipole4",
        "nLonsPerBlock": 36,
        "nLatsPerBlock": 18,
        "nAlts": 100,
        "LatRange": [10, 80],
        "AltRange": [80.0, 1000],
	      "LonRange": [0.0, 360.0]},
```

The dipole grid has both open field-lines and closed field-lines. The closed
field-lines are near the equator, while the open field-lines are near the poles.
The variable `MaxAlt` sets where the differentiation occurs - if the apex
height of all field-lines on this block are above this altitude, then it is open. 