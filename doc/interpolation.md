
# How to use the built in interpolator in Aether

Before attempting interpolation, you need to have the geogrid, some grids where you want to perform the interpolation (the size of the grid should match exactly with the geogrid, and the values in the grid should represent the quantity at the center of cells in the geogrid), and a set of points at which you want to interpolate the values.

First, place the longitude, latitude, and altitude of the points into three different vectors. For example, if you have points p1 (lon 180, lat 0, alt 10000), p2 (lon 90, lat -30, alt 15000), p3 (lon 270, lat 40, alt 10000), and p4 (lon 0, lat 0, alt 15000), then you should create three vectors as follows:

```bash
std::vector<precision_t> Lons = {180, 90, 270, 0};
std::vector<precision_t> Lats = {0, -30, 40, 0};
std::vector<precision_t> Alts = {10000, 15000, 10000, 15000};
```
Second, set the interpolation coefficients. Continuing with the previous example, assuming you have a geogrid named "geo_grid," you should call:

```bash
geo_grid.set_interpolation_coefs(Lons, Lats, Alts);
```

All subsequent calls to geo_grid.get_interpolation_values will use these coefficients until geo_grid.set_interpolation_coefs is called again (and the coefficients are updated to the newly set ones).

Finally, obtain the values at the set of points. Continuing with the previous example, if you want to perform interpolation on "arma_cube data1" and "arma_cube data2" with the set of points, you can execute:

```bash
std::vector<precision_t> ans1 = geo_grid.get_interpolation_values(data1);
std::vector<precision_t> ans2 = geo_grid.get_interpolation_values(data2);
```

The first element in "ans1" (i.e., ans1[0]) represents the interpolated value for the quantity "data1" at point p1 (lon 180, lat 0, alt 10000), and the second element is the result for point p2, and so on. This is the same for "ans2", except that it is for the quantity "data2".