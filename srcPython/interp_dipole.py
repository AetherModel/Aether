import numpy as np
import postAether

#----------------------------------------------------------------------------
# Earth dipole constants (matching src/tools.cpp and share/run/UA/inputs/orbits.csv)
#----------------------------------------------------------------------------

_DIPOLE_TILT = np.radians(10.0)       # co-latitude of magnetic pole [rad]
_DIPOLE_ROTATION = np.radians(270.0)  # geographic longitude of magnetic pole [rad]

#----------------------------------------------------------------------------
# Dipole-to-geographic interpolation helpers
#----------------------------------------------------------------------------

def dip_clean_coords(coords_by_block, use_magnetic=False):
    """Pack block list into a (3, nBlocks, nX, nY, nZ) coordinate array.

    Parameters
    ----------
    coords_by_block : list of dicts
        Block data as returned by read_block_files.
    use_magnetic : bool
        If True, use magnetic coordinates (mlon, invLat) for the horizontal
        axes instead of geographic (lon, lat).  The dipole grid is organized
        by magnetic coordinates, so searchsorted only works with mlon/invLat.
    """
    vars_list = coords_by_block[0]['vars']
    if use_magnetic:
        iX = postAether.find_var_index(vars_list, 'mlon')
        iY = postAether.find_var_index(vars_list, 'invLat')
    else:
        iX = postAether.find_var_index(vars_list, 'lon')
        iY = postAether.find_var_index(vars_list, 'lat')
    iAlt = postAether.find_var_index(vars_list, 'z')
    if iAlt < 0:
        iAlt = postAether.find_var_index(vars_list, 'alt')
    if iX < 0 or iY < 0 or iAlt < 0:
        needed = ('mlon, invLat' if use_magnetic else 'lon, lat') + ', z/alt'
        raise ValueError(
            'Bfield file is missing required variables (' + needed +
            '). Found vars: ' + str(vars_list))
    ishape = coords_by_block[0][iX].shape
    grid = np.zeros([3, len(coords_by_block), ishape[0], ishape[1], ishape[2]])
    for n, blk in enumerate(coords_by_block):
        grid[0, n] = blk[iX]
        grid[1, n] = blk[iY]
        grid[2, n] = blk[iAlt]
    return grid

def dip_get_block_lims(grid, nGCs=2):
    """Return per-block [(lon_min, lon_max), (lat_min, lat_max)] limits.

    Uses the full block extent (including ghost cells) for the search region.
    Interior limits are also returned so dip_find_one_pt can prefer the block
    where the point is deepest in the interior, avoiding ghost-cell data.
    """
    lims = []
    int_lims = []
    for iB in range(grid.shape[1]):
        lims.append([
            (np.min(grid[0, iB, :, 0, 0]), np.max(grid[0, iB, :, 0, 0])),
            (np.min(grid[1, iB, 0, :, 0]), np.max(grid[1, iB, 0, :, 0])),
        ])
        int_lims.append([
            (np.min(grid[0, iB, nGCs:-nGCs, 0, 0]),
             np.max(grid[0, iB, nGCs:-nGCs, 0, 0])),
            (np.min(grid[1, iB, 0, nGCs:-nGCs, 0]),
             np.max(grid[1, iB, 0, nGCs:-nGCs, 0])),
        ])
    return lims, int_lims

def dip_find_one_pt(pt, src_grid, src_lims, src_int_lims=None, nGCs=2):
    """Locate which source block contains pt; return (found, block, left, down).

    When multiple blocks contain the point (overlap / ghost-cell region),
    prefer the block where the point is deepest in the interior to avoid
    using potentially inconsistent ghost-cell data.

    Returns
    -------
    found : bool
    onblock, left, down : int  (meaningful only when found=True)
    """
    onblock, found = 0, False
    best_depth = -1.0
    for n, lim in enumerate(src_lims):
        if (lim[0][0] < pt[0] < lim[0][1]) and (lim[1][0] < pt[1] < lim[1][1]):
            if src_int_lims is not None:
                ilim = src_int_lims[n]
                # Depth = min distance to any interior edge (normalised)
                sx = ilim[0][1] - ilim[0][0]
                sy = ilim[1][1] - ilim[1][0]
                dx = min(pt[0] - ilim[0][0], ilim[0][1] - pt[0]) / sx if sx > 0 else 0.5
                dy = min(pt[1] - ilim[1][0], ilim[1][1] - pt[1]) / sy if sy > 0 else 0.5
                depth = min(dx, dy)
            else:
                depth = 0.0  # no preference info; keep last-match behaviour
            if depth > best_depth:
                best_depth = depth
                onblock, found = n, True
    if not found:
        return False, 0, 0, 0
    left = int(np.searchsorted(src_grid[0, onblock, :, 0, 0], pt[0]))
    down = int(np.searchsorted(src_grid[1, onblock, 0, :, 0], pt[1]))
    return True, onblock, left, down

def dip_find_all_corners(dst_grid, src_grid, src_lims, src_int_lims=None, nGCs=2):
    """For every point in dst_grid find the 8 surrounding src_grid corners.

    Because altitude varies dramatically across field lines (different Y-indices),
    we compute the altitude bracket independently for each of the four horizontal
    positions (l,d), (lm,d), (l,dm), (lm,dm).  Using a single (bd,bu) for all
    four positions would give badly wrong altitude weights near the equator.

    Points that fall outside the dipole domain (dip_find_one_pt returns
    found=False) get all corners set to zero so they can be detected later.
    """
    _, nBlocks, nX, nY, nZ = dst_grid.shape
    nX_s, nY_s, nZ_s = src_grid.shape[2], src_grid.shape[3], src_grid.shape[4]
    corners  = np.zeros([3, nBlocks, nX, nY, nZ, 8], dtype=int)
    onblocks = np.zeros([nBlocks, nX, nY, nZ], dtype=int)
    inside   = np.zeros([nBlocks, nX, nY, nZ], dtype=bool)

    def _alt_bracket(blk, il, jd, target_alt):
        """Return (below, above, in_range) for the altitude column at (blk,il,jd).

        in_range is False when target_alt is outside [col[0], col[-1]], meaning
        the point is outside the dipole grid's altitude domain at this field line.
        Ghost cells are excluded from both ends of the altitude column.
        """
        col = src_grid[2, blk, il, jd, nGCs:-nGCs]
        above = int(np.searchsorted(col, target_alt))
        below = max(0, above - 1)
        above = min(above, len(col) - 1)
        in_range = (target_alt >= col[0]) and (target_alt <= col[-1])
        # Offset indices back to full-array coordinates
        return below + nGCs, above + nGCs, in_range

    # Interior index range (excludes ghost cells)
    lo_x, hi_x = nGCs, nX_s - nGCs - 1
    lo_y, hi_y = nGCs, nY_s - nGCs - 1

    for iB in range(nBlocks):
        for iX in range(nX):
            for iY in range(nY):
                for iZ in range(nZ):
                    pt = dst_grid[:, iB, iX, iY, iZ]
                    found, blk, l, d = dip_find_one_pt(pt, src_grid, src_lims, src_int_lims)
                    if not found:
                        # Leave corners as zeros; detected downstream as "outside"
                        continue
                    onblocks[iB, iX, iY, iZ] = blk
                    l  = max(lo_x, min(l,  hi_x))
                    d  = max(lo_y, min(d,  hi_y))
                    lm, dm = max(lo_x, l - 1), max(lo_y, d - 1)
                    # Per-column altitude brackets; if any column is altitude-out-of-range,
                    # treat the whole point as outside (leave corners at zeros).
                    bd_ld,   bu_ld,   ok_ld   = _alt_bracket(blk, l,  d,  pt[2])
                    bd_lmd,  bu_lmd,  ok_lmd  = _alt_bracket(blk, lm, d,  pt[2])
                    bd_lmdm, bu_lmdm, ok_lmdm = _alt_bracket(blk, lm, dm, pt[2])
                    bd_ldm,  bu_ldm,  ok_ldm  = _alt_bracket(blk, l,  dm, pt[2])
                    if not (ok_ld and ok_lmd and ok_lmdm and ok_ldm):
                        continue  # altitude out of range at ≥1 column → outside
                    corners[:, iB, iX, iY, iZ, 0] = [l,  d,  bd_ld]
                    corners[:, iB, iX, iY, iZ, 1] = [lm, d,  bd_lmd]
                    corners[:, iB, iX, iY, iZ, 2] = [lm, dm, bd_lmdm]
                    corners[:, iB, iX, iY, iZ, 3] = [l,  dm, bd_ldm]
                    corners[:, iB, iX, iY, iZ, 4] = [l,  d,  bu_ld]
                    corners[:, iB, iX, iY, iZ, 5] = [lm, d,  bu_lmd]
                    corners[:, iB, iX, iY, iZ, 6] = [lm, dm, bu_lmdm]
                    corners[:, iB, iX, iY, iZ, 7] = [l,  dm, bu_ldm]
                    inside[iB, iX, iY, iZ] = True
    return corners, onblocks, inside

def dip_compute_weights_and_idx(src_grid, corners, onblocks, dst_grid, inside=None):
    """Vectorised trilinear weights and flat source indices from cached corners.

    Corner order (set by dip_find_all_corners):
      0:(l,d,bd_ld)    1:(lm,d,bd_lmd)    2:(lm,dm,bd_lmdm)  3:(l,dm,bd_ldm)
      4:(l,d,bu_ld)    5:(lm,d,bu_lmd)    6:(lm,dm,bu_lmdm)  7:(l,dm,bu_ldm)

    Because altitude varies significantly across field lines, each horizontal
    column (l,d), (lm,d), (lm,dm), (l,dm) has its own altitude bracket and
    therefore its own altitude fraction dz_col.  This makes the interpolation
    exact for any field that is separable in (lon, lat, alt).

    Parameters
    ----------
    src_grid : np.ndarray, shape (3, nB_src, nX_s, nY_s, nZ_s)
    corners  : np.ndarray, shape (3, nB_dst, nX, nY, nZ, 8), int
    onblocks : np.ndarray, shape (nB_dst, nX, nY, nZ), int
    dst_grid : np.ndarray, shape (3, nB_dst, nX, nY, nZ)
    inside   : np.ndarray, shape (nB_dst, nX, nY, nZ), bool, optional
        Mask of points inside the source domain.  Outside points get zero
        weights so the interpolated result is 0.0 there.

    Returns
    -------
    weights  : np.ndarray, shape (nB_dst, nX, nY, nZ, 8), float64
    flat_idx : np.ndarray, shape (nB_dst, nX, nY, nZ, 8), int64
    """
    _, nX_s, nY_s, nZ_s = src_grid.shape[1:]

    # Flat index into the ravelled source array for each of the 8 corners
    ob       = onblocks[:, :, :, :, np.newaxis]  # broadcast over 8 corners
    flat_idx = (ob         * nX_s * nY_s * nZ_s +
                corners[0] *        nY_s * nZ_s +
                corners[1] *               nZ_s +
                corners[2])                       # (nB_dst, nX, nY, nZ, 8)

    # Coordinates of the 8 source corners for every destination point
    src_flat      = src_grid.reshape(3, -1)       # (3, nB_src*nX_s*nY_s*nZ_s)
    corner_coords = src_flat[:, flat_idx]          # (3, nB_dst, nX, nY, nZ, 8)

    # Longitude fraction dx: lon only varies with l, so min/max over all 8 is correct
    def _frac_horiz(dim):
        lo   = corner_coords[dim].min(axis=-1)
        hi   = corner_coords[dim].max(axis=-1)
        span = hi - lo
        return np.where(span > 1e-12,
                        np.clip((dst_grid[dim] - lo) / span, 0.0, 1.0),
                        0.5)

    dx = _frac_horiz(0)   # (nB, nX, nY, nZ)
    dy = _frac_horiz(1)

    # Per-column altitude fraction:
    # corners 0-3 are the LOWER brackets; corners 4-7 are the UPPER brackets
    # for the four horizontal positions (l,d), (lm,d), (lm,dm), (l,dm).
    alt_lo  = corner_coords[2, ..., :4]           # (nB, nX, nY, nZ, 4)
    alt_hi  = corner_coords[2, ..., 4:]           # (nB, nX, nY, nZ, 4)
    span_z  = alt_hi - alt_lo
    tgt_alt = dst_grid[2, ..., np.newaxis]        # broadcast over 4 columns
    dz_col  = np.where(span_z > 1e-12,
                       np.clip((tgt_alt - alt_lo) / span_z, 0.0, 1.0),
                       0.5)                        # (nB, nX, nY, nZ, 4)

    # Horizontal weights for the 4 column positions: (l,d), (lm,d), (lm,dm), (l,dm)
    h_w = np.stack([
        dx     * dy,        # column 0: (l,  d )
        (1-dx) * dy,        # column 1: (lm, d )
        (1-dx) * (1-dy),    # column 2: (lm, dm)
        dx     * (1-dy),    # column 3: (l,  dm)
    ], axis=-1)             # (nB, nX, nY, nZ, 4)

    # Final 8 weights: lower half uses (1-dz_col), upper half uses dz_col
    weights = np.concatenate([
        h_w * (1.0 - dz_col),    # corners 0-3 (lower)
        h_w *        dz_col,     # corners 4-7 (upper)
    ], axis=-1)                   # (nB_dst, nX, nY, nZ, 8)

    # Zero out weights for points outside the source domain
    if inside is not None:
        weights[~inside] = 0.0

    return weights, flat_idx

def dip_do_interpolate_fast(weights, flat_idx, src_data):
    """Interpolate one variable using precomputed weights and flat indices.

    Parameters
    ----------
    weights  : np.ndarray, shape (nB_dst, nX, nY, nZ, 8)
    flat_idx : np.ndarray, shape (nB_dst, nX, nY, nZ, 8), int
    src_data : np.ndarray, shape (nB_src, nX_src, nY_src, nZ_src)

    Returns
    -------
    np.ndarray, shape (nB_dst, nX, nY, nZ)
    """
    corner_vals = src_data.flatten()[flat_idx]    # (nB_dst, nX, nY, nZ, 8)
    return (weights * corner_vals).sum(axis=-1)

def dip_get_corners(dst_grid, src_grid):
    """Return corner indices, onblocks, and inside mask."""
    print('  --> Computing interpolation corners...')
    src_lims, src_int_lims = dip_get_block_lims(src_grid)
    # The source (dipole) grid stores mlon in [0, 360] while the destination
    # (geographic) grid may use [-180, 180].  Wrap dst mlon into the source
    # grid's range so the point-finding search works correctly.
    src_mlon_lo = min(lim[0][0] for lim in src_lims)
    src_mlon_hi = max(lim[0][1] for lim in src_lims)
    if dst_grid[0].min() < src_mlon_lo or dst_grid[0].max() > src_mlon_hi:
        dst_grid = dst_grid.copy()
        dst_grid[0] = src_mlon_lo + (dst_grid[0] - src_mlon_lo) % 360.0
    corners, onblocks, inside = dip_find_all_corners(dst_grid, src_grid, src_lims, src_int_lims)
    n_inside = inside.sum()
    n_total = inside.size
    print(f'  --> Points inside dipole domain: {n_inside}/{n_total} ({100*n_inside/n_total:.1f}%)')
    return corners, onblocks, inside


def prefetch_weights(basefiles_info=[], isVerbose=False):
    """Precompute interpolation weights & indices for dip->geo conversion.

    Parameters
    ----------
    basefiles_info : list
        Info from get_base_files().  Must include both dipole and geographic
        background field files (3DBF).
    isVerbose : bool

    Returns
    -------
    dict
        Contains 'weights', 'flat_idx', 'inside', 'geoGridGeo' — everything
        needed by apply_cached_interp() to interpolate a dipole file.
    """
    bfield_dip = None
    bfield_geo = None
    for fileInfo in basefiles_info:
        if '3DBF' not in fileInfo['coreFile']:
            continue
        gs = postAether.get_gridshape(fileInfo['coreFile'], fileInfo['isNetCDF'])
        if gs == 'dipole' and bfield_dip is None:
            bfield_dip, _ = postAether.read_block_files(fileInfo['coreFile'],
                                             fileInfo['isNetCDF'],
                                             isVerbose=isVerbose)
        elif gs != 'dipole' and bfield_geo is None:
            bfield_geo, _ = postAether.read_block_files(fileInfo['coreFile'],
                                             fileInfo['isNetCDF'],
                                             isVerbose=isVerbose)

    if bfield_dip is None or bfield_geo is None:
        print('  Warning: Need both dipole and geographic BF files to '
              'precompute interpolation weights.')
        return None

    dipGrid = dip_clean_coords(bfield_dip, use_magnetic=True)
    geoGrid = dip_clean_coords(bfield_geo, use_magnetic=True)
    geoGridGeo = dip_clean_coords(bfield_geo, use_magnetic=False)

    corners, onblocks, inside = dip_get_corners(geoGrid, dipGrid)

    if isVerbose:
        print('  --> Building interpolation weight matrix...')
    weights, flat_idx = dip_compute_weights_and_idx(
        dipGrid, corners, onblocks, geoGrid, inside)

    return {
        'weights': weights,
        'flat_idx': flat_idx,
        'inside': inside,
        'geoGridGeo': geoGridGeo,
        'geoBFData': bfield_geo,
    }



#----------------------------------------------------------------------------
# apply cached interpolation weights
#----------------------------------------------------------------------------

def apply_cached_interp(cache, allBlockData, isVerbose=True):
    """Interpolate dipole block data using precomputed weights from prefetch_weights().

    Parameters
    ----------
    cache : dict
        From prefetch_weights(): weights, flat_idx, inside, geoGridGeo, geoBFData.
    allBlockData : list of dicts
        Dipole block data.

    Returns
    -------
    list of dicts
        Block data on the geographic grid (same format as interp_dip2geo output).
    """
    weights = cache['weights']
    flat_idx = cache['flat_idx']
    geoGridGeo = cache['geoGridGeo']
    geoBFData = cache['geoBFData']

    nBlocks_src = len(allBlockData)
    nBlocks_dst = geoGridGeo.shape[1]
    nVars = len(allBlockData[0]['vars'])

    iLon = postAether.find_var_index(allBlockData[0]['vars'], 'lon')
    iLat = postAether.find_var_index(allBlockData[0]['vars'], 'lat')
    iAlt = postAether.find_var_index(allBlockData[0]['vars'], 'alt')
    if iAlt < 0:
        iAlt = postAether.find_var_index(allBlockData[0]['vars'], 'z')
    coord_idxs = {i for i in [iLon, iLat, iAlt] if i >= 0}

    outBlockData = []
    for iB in range(nBlocks_dst):
        outBlockData.append(
            {key: allBlockData[0][key]
             for key in allBlockData[0] if isinstance(key, str)})
        outBlockData[-1]['gridshape'] = geoBFData[0].get('gridshape', 'latlon')

    if isVerbose:
        print('  --> Interpolating ', nVars, ' variables (cached weights)...')

    for var_idx in range(nVars):
        if var_idx in coord_idxs:
            continue
        src_data = np.array([allBlockData[b][var_idx]
                             for b in range(nBlocks_src)])
        result = dip_do_interpolate_fast(weights, flat_idx,
                                                       src_data)
        for iB in range(nBlocks_dst):
            outBlockData[iB][var_idx] = result[iB]

    for iB in range(nBlocks_dst):
        if iLon >= 0:
            outBlockData[iB][iLon] = geoGridGeo[0, iB]
        if iLat >= 0:
            outBlockData[iB][iLat] = geoGridGeo[1, iB]
        if iAlt >= 0:
            outBlockData[iB][iAlt] = geoGridGeo[2, iB]

    return outBlockData


def interp_dip2geo(allBlockData, isVerbose=True):
    """Interpolate dipole-grid block data onto the geographic grid.

    Looks in the current directory for background field files (3DBF) and
    reads the gridShape attribute to identify dipole vs geographic grids,
    then trilinearly interpolates every variable in allBlockData onto
    the geographic grid.  Interpolation weights are cached to disk so
    repeated calls on the same grids are fast.

    Parameters
    ----------
    allBlockData : list of dicts
        Block data as returned by read_block_files, on the dipole grid.

    Returns
    -------
    list of dicts
        Block data on the geographic grid, compatible with write_netcdf /
        write_hdf5.  Coordinate variables (lon, lat, alt/z) are replaced
        with geographic coordinates from the geographic background field.
    """
    filesInfo = postAether.get_base_files()
    bfield_files = [f for f in filesInfo if '3DBF' in f['coreFile']]

    # Classify background field files by gridShape attribute
    bfDip = []
    bfGeo = []
    for fInfo in bfield_files:
        gridshape = postAether.get_gridshape(fInfo['coreFile'],
                                             fInfo['isNetCDF'])
        if gridshape == 'dipole':
            bfDip.append(fInfo)
        else:
            bfGeo.append(fInfo)

    if not bfDip:
        print('  Warning: No dipole background field file found; '
              'skipping interpolation.')
        return allBlockData
    if not bfGeo:
        print('  Warning: No geographic background field file found; '
              'skipping interpolation.')
        return allBlockData

    if isVerbose:
        print('  --> Reading dipole grid   :', bfDip[0]['coreFile'])
    dipBFData, _ = postAether.read_block_files(bfDip[0]['coreFile'],
                                    bfDip[0]['isNetCDF'], isVerbose=False)
    if isVerbose:
        print('  --> Reading geographic grid:', bfGeo[0]['coreFile'])
    geoBFData, _ = postAether.read_block_files(bfGeo[0]['coreFile'],
                                    bfGeo[0]['isNetCDF'], isVerbose=False)

    # Search grids use magnetic coordinates (mlon, invLat) because the dipole
    # grid is organised by magnetic coordinates — geographic lon is NOT monotonic
    # along the dipole X-axis, so searchsorted would give wrong results.
    dipGrid = dip_clean_coords(dipBFData, use_magnetic=True)
    geoGrid = dip_clean_coords(geoBFData, use_magnetic=True)
    # Geographic coordinates for the output (to replace lon/lat/z in the result)
    geoGridGeo = dip_clean_coords(geoBFData, use_magnetic=False)

    corners, onblocks, inside = dip_get_corners(geoGrid, dipGrid)

    nBlocks_src = len(allBlockData)
    nBlocks_dst = geoGrid.shape[1]
    nVars = len(allBlockData[0]['vars'])

    # Identify coordinate variable indices so we don't interpolate them
    iLon = postAether.find_var_index(allBlockData[0]['vars'], 'lon')
    iLat = postAether.find_var_index(allBlockData[0]['vars'], 'lat')
    iAlt = postAether.find_var_index(allBlockData[0]['vars'], 'alt')
    if iAlt < 0:
        iAlt = postAether.find_var_index(allBlockData[0]['vars'], 'z')
    coord_idxs = {i for i in [iLon, iLat, iAlt] if i >= 0}
    # # Also skip magnetic coordinate variables from interpolation
    # for mvar in ['mlon', 'invLat', 'mlat', 'mlt', 'radius']:
    #     idx = postAether.find_var_index(allBlockData[0]['vars'], mvar)
    #     if idx >= 0:
    #         coord_idxs.add(idx)

    # Build output blocks: copy string metadata
    outBlockData = []
    for iB in range(nBlocks_dst):
        outBlockData.append(
            {key: allBlockData[0][key]
             for key in allBlockData[0] if isinstance(key, str)})
        outBlockData[-1]['gridshape'] = geoBFData[0].get('gridshape',
                                                           'latlon')

    if isVerbose:
        print('  --> Building interpolation weight matrix...')
    weights, flat_idx = dip_compute_weights_and_idx(dipGrid, corners, onblocks, geoGrid, inside)

    if isVerbose:
        print('  --> Interpolating', nVars, 'variables onto geographic grid...')

    for var_idx in range(nVars):
        if var_idx in coord_idxs:
            continue  # replaced below with geographic coordinates
        src_data = np.array([allBlockData[b][var_idx] for b in range(nBlocks_src)])
        result = dip_do_interpolate_fast(weights, flat_idx, src_data)
        for iB in range(nBlocks_dst):
            outBlockData[iB][var_idx] = result[iB]

    # Replace coordinate variables with geographic grid coordinates
    # (geoGridGeo has lon/lat/z; geoGrid has mlon/invLat/z used for search)
    for iB in range(nBlocks_dst):
        if iLon >= 0:
            outBlockData[iB][iLon] = geoGridGeo[0, iB]
        if iLat >= 0:
            outBlockData[iB][iLat] = geoGridGeo[1, iB]
        if iAlt >= 0:
            outBlockData[iB][iAlt] = geoGridGeo[2, iB]

    return outBlockData


def convert_mag2geo(mlon, mlat):
    """Convert magnetic coordinates to geographic coordinates.

    Applies Earth's dipole rotation (270 deg) and tilt (10 deg) to
    transform from magnetic (mlon, mlat) to geographic (glon, glat).
    Follows the same convention as mag_to_geo() in src/tools.cpp:
      1. Rotate around Y by +tilt   (undo dipole tilt)
      2. Rotate around Z by +rotation (undo dipole rotation)

    Parameters
    ----------
    mlon, mlat : array_like
        Magnetic longitude and latitude in radians. Any shape.

    Returns
    -------
    glon, glat : np.ndarray
        Geographic longitude [0, 2pi) and latitude [-pi/2, pi/2] in radians.
    """
    mlon = np.asarray(mlon, dtype=float)
    mlat = np.asarray(mlat, dtype=float)

    # Magnetic (lon, lat) -> unit vector in magnetic Cartesian frame
    x = np.cos(mlat) * np.cos(mlon)
    y = np.cos(mlat) * np.sin(mlon)
    z = np.sin(mlat)

    ct, st = np.cos(_DIPOLE_TILT), np.sin(_DIPOLE_TILT)
    cr, sr = np.cos(_DIPOLE_ROTATION), np.sin(_DIPOLE_ROTATION)

    # Ry(+tilt)
    x1 = x * ct - z * st
    y1 = y
    z1 = x * st + z * ct

    # Rz(+rotation)
    x2 = x1 * cr + y1 * sr
    y2 = -x1 * sr + y1 * cr
    z2 = z1

    # Cartesian -> (glon, glat)
    glat = np.arcsin(np.clip(z2, -1.0, 1.0))
    glon = np.arctan2(y2, x2) % (2.0 * np.pi)
    return glon, glat


def convert_geo2mag(glon, glat):
    """Convert geographic coordinates to magnetic coordinates.

    Inverse of convert_mag2geo.  Follows the same convention as
    get_dipole() in src/dipole.cpp:
      1. Rotate around Z by -rotation
      2. Rotate around Y by -tilt

    Parameters
    ----------
    glon, glat : array_like
        Geographic longitude and latitude in radians. Any shape.

    Returns
    -------
    mlon, mlat : np.ndarray
        Magnetic longitude [0, 2pi) and latitude [-pi/2, pi/2] in radians.
    """
    glon = np.asarray(glon, dtype=float)
    glat = np.asarray(glat, dtype=float)

    # Geographic (lon, lat) -> unit vector in geographic Cartesian frame
    x = np.cos(glat) * np.cos(glon)
    y = np.cos(glat) * np.sin(glon)
    z = np.sin(glat)

    ct, st = np.cos(_DIPOLE_TILT), np.sin(_DIPOLE_TILT)
    cr, sr = np.cos(_DIPOLE_ROTATION), np.sin(_DIPOLE_ROTATION)

    # Rz(-rotation)
    x1 = x * cr - y * sr
    y1 = x * sr + y * cr
    z1 = z

    # Ry(-tilt)
    x2 = x1 * ct + z1 * st
    y2 = y1
    z2 = -x1 * st + z1 * ct

    # Cartesian -> (mlon, mlat)
    mlat = np.arcsin(np.clip(z2, -1.0, 1.0))
    mlon = np.arctan2(y2, x2) % (2.0 * np.pi)
    return mlon, mlat


def geo_to_invlat(glon_deg, glat_deg, alt_m, R_planet=6371.0e3):
    """Compute magnetic longitude and invariant latitude from geographic coords.

    Combines the dipole rotation (geo->mag) with the L-shell calculation
    to produce the invariant latitude, which is the Y-axis coordinate of
    Aether's dipole grid.

    Follows get_dipole() in src/dipole.cpp:
      invLat = sign(mlat) * acos(1/sqrt(L))
    where L = r / (R * cos^2(mlat)).

    Parameters
    ----------
    glon_deg, glat_deg : array_like
        Geographic longitude and latitude in degrees.
    alt_m : array_like
        Altitude above the surface in metres.
    R_planet : float
        Mean planet radius in metres (default: Earth).

    Returns
    -------
    mlon_deg : np.ndarray
        Magnetic longitude in degrees [0, 360).
    invlat_deg : np.ndarray
        Invariant latitude in degrees, signed by hemisphere.
    """
    mlon, mlat = convert_geo2mag(np.radians(glon_deg), np.radians(glat_deg))

    r = np.asarray(alt_m, dtype=float) + R_planet
    cos2 = np.cos(mlat) ** 2
    cos2_safe = np.where(cos2 > 1e-20, cos2, 1e-20)
    L = r / (R_planet * cos2_safe)
    invlat = np.sign(mlat) * np.arccos(np.clip(1.0 / np.sqrt(L), -1.0, 1.0))

    return np.degrees(mlon), np.degrees(invlat)


#----------------------------------------------------------------------------
# Command-line interface
#----------------------------------------------------------------------------

def unpack_to_blocks(data):
    """Convert read_aether_file output to list-of-block-dicts format.

    Postprocessed netcdf files may be consolidated (3D arrays) or
    multi-block (4D arrays with a leading block dimension).  This function
    detects which case applies and returns a list of block dicts with 3D
    arrays, matching the format expected by dip_clean_coords and friends.
    """
    # Find the first spatial variable (3D or 4D) to determine format
    for i in range(len(data['vars'])):
        arr = data[i]
        if arr.ndim == 4:
            # Multi-block: split along first axis
            nBlocks = arr.shape[0]
            blocks = []
            for iB in range(nBlocks):
                blk = {k: data[k] for k in data if isinstance(k, str)}
                for j in range(len(data['vars'])):
                    if data[j].ndim == 4:
                        blk[j] = data[j][iB]
                    else:
                        blk[j] = data[j]
                blocks.append(blk)
            return blocks
        elif arr.ndim == 3:
            # Single block / consolidated — wrap in a list
            return [data]
    # Fallback (no spatial vars found)
    return [data]


if __name__ == '__main__':
    import argparse
    import sys
    import os

    parser = argparse.ArgumentParser(
        description='Interpolate Aether dipole-grid output to geographic grid')
    parser.add_argument('-geo', default=None, type=str,
                        help='A (postprocessed) file with the target geographic grid.'
                        ' At minimum this needs to be a netcdf file with coordinates.')
    parser.add_argument('-outdir', default='',
                        help='Path to directory where interpolated files will be saved.'
                        ' Default is to put outputs next to input dipole_files.')
    parser.add_argument('-v', action='store_true',
                        help='Verbose output')
    parser.add_argument('-rm', action='store_true',
                        help='Delete source files?')
    parser.add_argument('-hdf5',
                        help='output HDF5 files?',
                        action="store_true")
    parser.add_argument('dipole_files', nargs='+',
                        help='Path to the file(s) to interpolate')
    args = parser.parse_args()

    if args.geo is None:
        print("Error: -geo argument is required (geographic grid file)")
        sys.exit(1)

    # ------------------------------------------------------------------
    # 1. Read geographic grid (target) and first dipole file (source grid)
    # ------------------------------------------------------------------
    if args.v:
        print(f"Reading geographic grid: {args.geo}")
    geo_data = postAether.read_aether_file(args.geo)
    geo_blocks = unpack_to_blocks(geo_data)

    if args.v:
        print(f"Reading dipole grid from: {args.dipole_files[0]}")
    dip_data_first = postAether.read_aether_file(args.dipole_files[0])
    dip_blocks_first = unpack_to_blocks(dip_data_first)

    # ------------------------------------------------------------------
    # 2. Build grids in magnetic coordinates and compute weights (once)
    #    Geographic files only have geographic coords (lon, lat, z), so
    #    we convert to magnetic (mlon, mlat) via the dipole rotation.
    #    For a centered dipole, mlat == invLat.
    # ------------------------------------------------------------------
    src_grid_geo = dip_clean_coords(dip_blocks_first, use_magnetic=False)
    src_mlon, src_invlat = geo_to_invlat(src_grid_geo[0], src_grid_geo[1],
                                         src_grid_geo[2])
    src_grid = np.array([src_mlon, src_invlat, src_grid_geo[2]])

    dst_grid_geo = dip_clean_coords(geo_blocks, use_magnetic=False)
    dst_mlon, dst_invlat = geo_to_invlat(dst_grid_geo[0], dst_grid_geo[1],
                                         dst_grid_geo[2])
    dst_grid = np.array([dst_mlon, dst_invlat, dst_grid_geo[2]])

    corners, onblocks, inside = dip_get_corners(dst_grid, src_grid)
    weights, flat_idx = dip_compute_weights_and_idx(
        src_grid, corners, onblocks, dst_grid, inside)

    # Variables that should not be interpolated (coordinates / metadata)
    coord_vars = {'lon', 'lat', 'z', 'alt', 'mlon', 'invLat',
                  'mlat', 'mlt', 'radius', 'time'}

    nBlocks_dst = dst_grid_geo.shape[1]

    # ------------------------------------------------------------------
    # 3. Process each dipole file
    # ------------------------------------------------------------------
    for dip_file in args.dipole_files:
        print(f"Interpolating: {dip_file}")

        dip_data = postAether.read_aether_file(dip_file)
        dip_blocks = unpack_to_blocks(dip_data)
        nBlocks_src = len(dip_blocks)
        nVars = len(dip_blocks[0]['vars'])

        # Build list of spatial variables to write (skip 'time' and
        # non-3D vars — write_netcdf handles 'time' separately via
        # data['time'] and expects integer-indexed 3D arrays starting at 0)
        spatial_vars = []    # (src_var_idx, var_name)
        for var_idx in range(nVars):
            var_name = dip_blocks[0]['vars'][var_idx]
            if var_name == 'time':
                continue
            if dip_blocks[0][var_idx].ndim != 3:
                continue
            spatial_vars.append((var_idx, var_name))

        # Map original var names to new output indices
        out_var_names = [name for _, name in spatial_vars]
        src_to_out = {src_idx: out_idx
                      for out_idx, (src_idx, _) in enumerate(spatial_vars)}

        # Initialize output blocks with metadata
        out_blocks = []
        for iB in range(nBlocks_dst):
            blk = {'vars': list(out_var_names),
                    'units': [dip_blocks[0]['units'][si]
                              for si, _ in spatial_vars],
                    'time': dip_blocks[0]['time']}
            if 'long_name' in dip_blocks[0]:
                blk['long_name'] = [dip_blocks[0]['long_name'][si]
                                    for si, _ in spatial_vars]
            out_blocks.append(blk)

        # Interpolate each non-coordinate variable; copy coords from geo grid
        for out_idx, (src_idx, var_name) in enumerate(spatial_vars):
            if var_name in coord_vars:
                # Replace coordinates with geographic grid values
                if var_name == 'lon':
                    for iB in range(nBlocks_dst):
                        out_blocks[iB][out_idx] = dst_grid_geo[0, iB]
                elif var_name == 'lat':
                    for iB in range(nBlocks_dst):
                        out_blocks[iB][out_idx] = dst_grid_geo[1, iB]
                elif var_name in ('z', 'alt'):
                    for iB in range(nBlocks_dst):
                        out_blocks[iB][out_idx] = dst_grid_geo[2, iB]
                else:
                    # Other coord vars (mlon, invLat, etc.) — skip
                    # Fill with zeros so write_netcdf has something
                    dst_shape = dst_grid_geo[0, 0].shape
                    for iB in range(nBlocks_dst):
                        out_blocks[iB][out_idx] = np.zeros(dst_shape)
                continue
            src_data = np.array([dip_blocks[b][src_idx]
                                 for b in range(nBlocks_src)])
            result = dip_do_interpolate_fast(weights, flat_idx, src_data)
            for iB in range(nBlocks_dst):
                out_blocks[iB][out_idx] = result[iB]

        # Determine output filename: insert _interp before the extension
        base = os.path.basename(dip_file)
        name, ext = os.path.splitext(base)
        out_name = name + '_interp' + ext
        if args.outdir:
            os.makedirs(args.outdir, exist_ok=True)
            out_path = os.path.join(args.outdir, out_name)
        else:
            out_path = os.path.join(os.path.dirname(dip_file) or '.', out_name)

        # Write output
        if args.hdf5:
            out_path = os.path.splitext(out_path)[0] + '.h5'
            postAether.write_hdf5(out_blocks, out_path, isVerbose=args.v)
        else:
            is_consolidated = (nBlocks_dst == 1)
            if is_consolidated:
                postAether.write_netcdf(out_blocks[0], out_path,
                                        isVerbose=args.v,
                                        isConsolidated=True)
            else:
                postAether.write_netcdf(out_blocks, out_path,
                                        isVerbose=args.v,
                                        isConsolidated=False)

        if args.rm:
            os.remove(dip_file)
            if args.v:
                print(f"  Deleted: {dip_file}")

    print("Done.")

