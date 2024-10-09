import matplotlib.pyplot as plt
import numpy as np

####        Set inputs          ####

# Set to None or '' to not save, just show (pan/zoom capabilities).
fig_save_path = None


nf = 100  # number of field lines
nz = 200  # number of grid cells on each field line
alt_min = 90  # in km, min altitude

# in degrees, start/end latitudes of field lines. between (90,0)
max_blat = 89.9
min_blat = 2.25

gams = 0.1
# point distribution along field lines.
# higher puts more pts at high altitudes.

baselat_spacing_factor = 6  # exponential factor for spacing baselats (uses cos() too)

# consts:
Re = 6371  # in km


# leave empty for "auto" (aaron b choose), or put in custom limits here.
# [left, right, bottom, top]
limits_left_plot = None  # default is whole image
limits_right_plot = None  # default is [-0.1, 2.0, -0.3, 1.75]


# ------------------------------------------------------------------------
# Main code is here:
# ------------------------------------------------------------------------

####        Useful Functions:          ####


def calc_baselats(blat_min, blat_max, n_f, factor):
    """Lay down base latitudes

    Args:
        blat_min (float): min latitude to start (positive), in degrees.
        blat_max (float): max latitude in degrees, positiv & less than 90.0
        nf (int): Number of field lines. Must be even!
        factor (float): Factor to use to space latitudes.

    Returns:
        numpy array: starting latitudes

    Notes: follows very similar approach to doi:10.1029/2000JA000035
    - Not exactly exponential, but exponential and a cosine!
    - Spaced according to B-field strength.

    """
    # Space out base latitudes:

    baselats = np.linspace(
        np.cos(np.deg2rad(blat_min ** (1 / factor))),
        np.cos(np.deg2rad(blat_max ** (1 / factor))),
        num=n_f,
    )

    # baselats are in southern hemisphere
    baselats = -1 * np.flip(np.rad2deg(np.arccos(baselats)) ** factor)

    return baselats


def calc_q(alt, theta, re=None, is_alt=False, isnt_re=False):
    """Calculate q (distance along field line) for a point

    Args:
        alt (float, array-like): altitude (kn or Re)
        theta (float or array-like): magnetic latitude (in degrees).
            0 at mag equator (measured from North Pole).
        re (float): radius of earth in km. only used if alt is altitude and/or in km
        is_alt (bool): altitude is altitude? False means it's radius. Default is False
        isnt_re (bool): altitude is in km? False means altitude has units [Re]. Default is False

    Returns:
        (array or float): q-value for the given (alt, theta) point.

    Notes:
    - make sure theta has units degrees and is measured from north pole
    - See <doi:10.1029/2000JA000035> for more information.
    - is_alt & isnt_re are from debugging & aren't used anymore.

    """

    if is_alt:  # convert altitude to radius
        alt = alt + re
    if isnt_re:  # convert km to re
        alt = alt / re

    return np.cos(np.deg2rad(90 - theta)) / (alt**2)


def calc_p(alt, theta, re=None, is_alt=False, isnt_re=False):
    """Calculate p-value (l-shell) for a given altitude & latitude

    Args:
        alt (float): altitude (in km or Re), (above surface or from origin)
        theta (float): latitude, in degrees
        re (float): earth radius in km. optional.
        is_alt (bool): is alt altitude? if it's radius, set to False (default).
        isnt_re (bool): altitude is in km? False (default) means altitude has units [Re].

    Returns:
        (array or float): p-value. Same as l-Shell, in Re.

    Notes:
    - make sure theta has units degrees and is measured from north pole
    - See <doi:10.1029/2000JA000035> for more information.
    - is_alt & isnt_re are from debugging & aren't used anymore.

    """
    if is_alt:
        alt = alt + re
    if isnt_re:
        alt = alt / re
    return alt / (np.sin(np.deg2rad(90 - theta)) ** 2)


def qp_solve(q, p):
    """Solve for radius given (q,p) dipole coordinates
    Methodology from <https://arxiv.org/pdf/physics/0606044?

    Args:
        q (float): Distance from mag equator of dipole. Negative in south hemisphere.
        p (float): L-shell, basically. Units of Re

    Returns:
        float: R, distance from origin (of coord system) for point (q, p)
    """

    term0 = 256.0 / 27.0 * q**2 * p**4
    term1 = (1.0 + np.sqrt(1.0 + term0)) ** (2.0 / 3.0)
    term2 = term0 ** (1.0 / 3.0)
    term3 = 0.5 * ((term1**2 + term1 * term2 + term2**2) / term1) ** (3.0 / 2.0)
    new_r = p * (4.0 * term3) / (1.0 + term3) / (1.0 + np.sqrt(2.0 * term3 - 1.0))

    return new_r


####        The main stuff:          ####


def calc_exp_grid(nf0, nz0, altminre, baselats, pvals, gamma):
    """Exponential Grid laydown, keeps grid parallel & perpendicular to B

    Args:
        nf0 (int): number of field lines
        nz0 (int): number of points along each field line. MUST BE EVEN!
        altminre (float): min altitude (in Re from center of earth) to trace from
        baselats (list/array): latitudes to start at (from calc_baselats)
        pvals (flaot): p-values (dipole coords, so it's L-shells) for all field lines
        gamma (float): factor used to space points along field line.
            Lower values increase point density @ low altitudes

    Returns:
        [np.array, np.array]: (nf, nz) dimensionsl arrays of:
            - Latitudes (in deg)
            - Radii (in re) of grid

    ** This can be vectorized & shortened A LOT.
    Left in this state for readability and in case things need to be changed.

    """

    lats_2d = []
    rs_2d = []

    nzh = int(nz0 / 2)

    for f_iter in range(nf0):
        # q in sh & nh
        q_S = calc_q(altminre, baselats[f_iter], is_alt=False, isnt_re=False)
        q_N = calc_q(altminre, -baselats[f_iter], is_alt=False, isnt_re=False)

        # linear spacing - made it really readable, could be done "cleaner"
        delqp = (q_N - q_S) / nz0
        qp0 = []
        for i in range(nz0):
            qp0.append(q_S + i * delqp)

        # exp grid laydown, (same for all calls here, speed not an issue though)
        delqp = altminre * delqp
        f00s = []
        for i in range(nz0):
            f00s.append(gamma + (1 - gamma) * np.exp(-(((i - nzh) / (nz0 / 10)) ** 2)))

        # spacing according to sinh function
        ft = []
        for i in range(nz0):
            fb0 = (1 - f00s[i]) / np.exp(-q_S / delqp - 1)
            fa = f00s[i] - fb0
            ft.append(fa + fb0 * np.exp(-(qp0[i] - q_S) / delqp))

        # q values, from south -> equator
        qpnew = []
        for i in range(nzh):
            delq = qp0[i] - q_S
            qpnew.append(q_S + ft[i] * delq)

        # qpnew is from south-equator. extend it to north pole,
        # so *-1 & reverse order so it is ascending
        qpnew.extend(np.flip(np.array(qpnew)) * -1)

        ilats = []
        irs = []

        for i in range(nz0):
            # use qpsolve to get r from (q,p)
            irs.append(qp_solve(qpnew[i], pvals[f_iter]))
            # Use dipole equations to get lat from q and r
            # q = cos(theta)/r**2
            ilats.append(np.rad2deg(np.arcsin(qpnew[i] * irs[-1] ** 2)))

        # Put into 2-D; lists easier & then return numpy
        lats_2d.append(ilats)
        rs_2d.append(irs)

    return np.array(lats_2d), np.array(rs_2d)


####        Actual computation:          ####


AltMinRe = (alt_min + Re) / Re  # alt min in Re

# baselats
baselats = calc_baselats(min_blat, max_blat, nf, baselat_spacing_factor)

# l-shells
pvals = calc_p(AltMinRe, baselats)

# take those, make lats & radii
lats, rs = calc_exp_grid(nf, nz, AltMinRe, baselats, pvals, gams)


print("making plot")

# change variables in case anyone wants to make different plots
xs = lats[:, :]
ys = rs[:, :]


fig, ax = plt.subplots(1, 2, figsize=(7, 7))

# scatter points, same color is same field line
for x, y in zip(xs, ys):
    ax[0].scatter(y * np.cos(np.deg2rad(x)), y * np.sin(np.deg2rad(x)))
    ax[1].scatter(y * np.cos(np.deg2rad(x)), y * np.sin(np.deg2rad(x)))

# change variables again, overwrite previous xs, ys
# Take every 4th field line so it's more clear to see things
xs = lats[:, :-8:4]
ys = rs[:, :-8:4]
# black dashed lines, same nz value at different nf's
for x, y in zip(xs.T, ys.T):
    ax[0].plot(
        y * np.cos(np.deg2rad(x)), y * np.sin(np.deg2rad(x)), linestyle="--", color="k"
    )
    ax[1].plot(
        y * np.cos(np.deg2rad(x)), y * np.sin(np.deg2rad(x)), linestyle="--", color="k"
    )

# make square-ish
ax[0].set_aspect(1)
ax[1].set_aspect(1)

# custom limits?
if limits_left_plot:
    ax[1].set_xlim(limits_left_plot[0], limits_left_plot[1])
    ax[1].set_ylim(limits_left_plot[2], limits_left_plot[3])


if limits_right_plot:
    ax[1].set_xlim(limits_right_plot[0], limits_right_plot[1])
    ax[1].set_ylim(limits_right_plot[2], limits_right_plot[3])
else:
    ax[1].set_xlim(-0.1, 2)
    ax[1].set_ylim(-0.3, 1.75)

# save or show:
if fig_save_path:
    plt.savefig(fig_save_path)
else:
    plt.show()
plt.close("all")
