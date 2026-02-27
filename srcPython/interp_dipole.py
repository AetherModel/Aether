import numpy as np
import postAether

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
    ishape = coords_by_block[0][0].shape
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
        Dipole atmosphere block data.

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
