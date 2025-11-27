import matplotlib.pyplot as plt
import numpy as np

####        Set inputs          ####

# To save plot, change last line of the code, otherwise it is just shown

# Number of lats/alts (without ghost cells)
nLatsPerBlock_in = 12 
nAltsPerBlock_in = 12 

# in degrees, where to begin & end grid between (90,0)
# Grid is mirrored across N/S hemisphere
max_blat = 85
min_blat = 12

# In km (above surface)
min_alt = 80
max_alt = 800

# Number of "blocks" to simulate - i.e. # of processors in Aether run (must be >4 & even)
nBlocks = 6

nGCs = 2

# consts:
Re_KM = 6371  # in km
cPI = np.pi


########                  ~~~~~~~~~~~~~~~~~~~~                  ########
# OUTLINE:
# - constants & inputs (above)
# - def main (the function to create the grid)
# - def all conversions
# - def the ploting function
# - Run script:
#   - make the grid as Aether would
#   - call the plotting function

def main(alt_minRE, alt_maxRE, lat_min, lat_max, origins, extent, nLatsPerBlock, nAltsPerBlock):

    pcenters = np.zeros([len(origins), nLatsPerBlock, nAltsPerBlock])
    qcenters = np.zeros([len(origins), nLatsPerBlock, nAltsPerBlock])

    # Loop through QT's origins
    for n, origin in enumerate(origins):
        close_this_block = False
        isSouth = False
        # If we're in south hemisphere, flip everything. WIll be undone later.
        if origin < -0.01:
            isSouth = True
            origin  = -1*origin - extent

        lat0 = (2*(lat_max - lat_min))*origin    
        dlat = extent * (2*(lat_max - lat_min)) / (nLatsPerBlock - nGCs*2)
        
        # Put latitudes down evenly (centers & corners)
        # - This forms the invariant latitudes which field lines must pass thru
        lat1d = []
        lat1d_co = []   
        pcenters1d = np.empty(nLatsPerBlock)
        qs = np.empty((nLatsPerBlock, nAltsPerBlock))
        pcenters2d = np.empty((nLatsPerBlock, nAltsPerBlock))

        for i in range(nLatsPerBlock):
            lat1d.append(lat0 + (i - nGCs + 0.5) * dlat + lat_min)

        # IF touching pole, put last ghost cell at 89.9 degrees & the 2nd to last 1/2 way there.
        if origin + extent > 0.49:
            lat1d[-1] = 89.9
            lat1d[-2] = (lat1d[-1] + lat1d[-2]) /2

        for i in range(nLatsPerBlock):
            pcenters1d[i] = alt_minRE / (np.sin(cPI/2 - np.deg2rad(lat1d[i]))**2)
            # Easier to save later if we get this:
            pcenters2d[i, :] = alt_minRE / (np.sin(cPI/2 - np.deg2rad(lat1d[i]))**2)

        # Corners (only used here to determine if we need to close > 1 block per hemisphere)
        for i in range(nLatsPerBlock+1): 
            lat1d_co.append(lat0 + (i - nGCs) * dlat + lat_min)
        pcorners = alt_minRE / np.sin(cPI/2 - np.deg2rad(lat1d_co)) **2

        ## Determine if field lines should close. There are two conditions:
        # - If the lowest l-shell in this block < altMin
        if np.min(pcorners) < alt_maxRE: 
            close_this_block = True
        # - Or if we are touching the equator
        if origin < 0.01: # NH equator
            close_this_block + True
        if np.abs(extent + origin) < 0.01: # SH equator
            close_this_block = True

        ## Setting up the q-values...

        # The idea here is that we either want the field line to close (wrap over equator)
        # or to have its boundaries entirely within min/max alt. 
        # We do not want field lines ending before max_alt, and vice-versa.
        # By definition, q=0 at the equator and +/- infinity at the N/S poles, so:
        # - Q_max is obtained from the minimum altitude point on the highest latitude field line
        q_max_center = rp2q(alt_minRE, pcenters1d[-1])
        if close_this_block:
            q_min_center = 0 # if the block is closed, q_min = 0. This is the equator!
        else:
            # If open, q_min is the highest altitude point on the lowest latitude field line
            q_min_center = rp2q(alt_maxRE, pcenters1d[0])

        
        delQ = (q_max_center - q_min_center) / (nAltsPerBlock - nGCs*2) 
        for iAlt in range(nAltsPerBlock):
            qs[:, iAlt] = ((q_min_center + (iAlt - nGCs + 0.5) * delQ))
        
        # If we were in South hemisphere, multiply by -1
        # And put data in the same order as we get back from Aether
        if isSouth:
            qs = -1.0*qs
            pcenters2d = np.flip(pcenters2d, axis=0)

        qcenters[n,:] = np.flip(qs, axis=1)
        pcenters[n,:] = pcenters2d

    return qcenters, pcenters


####        Useful Functions for conversions:          ####
####        - Not all used... 
####        - Format: in2out, in as few letters as necessary
####          example: cart to geo "xy2rt": (x,y) --> (r, theta)

## NOTE: theta for the dipole coordinate system is defined as co-latitude, not latitude
## Thus, we do (cPI-theta) for rt2(q/p).
## Then things are kept as-is, until conversion back to spherical when 
## colatitude is again considered

def rt2q(r, t):
    return np.cos(cPI/2 - t)/r**2

def rt2p(r, t):
    return r/(np.sin(cPI/2 - t)**2)

def rt2qp(r, t):
    q = rt2q(r, t)
    p = rt2p(r, t)
    return q, p

def rt2xy(r, t):
    x = r*np.cos(t)
    y = r*np.sin(t)
    return x, y


def rq2t(r,q):
    return np.arcsin(q * r**2)
def rp2t(r, p):
    return np.arccos(np.sqrt(r/p))

def rp2q(r, p):
    return np.sqrt((1-r/p)/r**4)

def qp2xy(q, p):
    r_ = qp_solve(q, p)
    t_ = rq2t(r_, q)
    return rt2xy(r_, t_)

def tp2r(t, p):
    return p * (np.cos(t)**2)

def alt2r(alt, re):
    return (alt + re)/re

def r2alt(r, re):
    return r*re - re


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


def make_plot(qs, ps, alt_min_RE, Re_km=6371,
              abs_bot = False # Take the absolute value of latitude on bottom plot?
              ):

    rs = qp_solve(qs, ps)
    ts = rq2t(rs, qs)
    
    fig = plt.figure(figsize=(8,11))
    
    gs = plt.GridSpec(8,9)
    
    ax0 = fig.add_subplot(gs[:5,:3])
    for x,y in zip(*qp2xy(qs, ps)):
        ax0.scatter(x,y, s=5)
    
    xlim, ylim = ax0.get_xlim(), ax0.get_ylim()
    
    circle1 = plt.Circle((0, 0), 1, color='k', alpha = .7)
    ax0.add_patch(circle1)
    
    ax0.set_ylim(ylim)
    ax0.set_xlim(xlim)
    ax0.set_aspect(1)
    ax0.set_title('in Re:')

    ax1 = fig.add_subplot(gs[:2,4:])
    counts, _, _ = ax1.hist(rs.flatten(), bins=60)
    ax1.vlines(alt_min_RE, 0, max(counts)*1.1, linestyle = '--', alpha=.7, color='k')
    ax1.set_title(f"{np.sum(rs < alt_min_RE) / np.prod(rs.shape)*100:.2f}% of points below min_alt\n"
                 f"{np.sum(rs < 1) / np.prod(rs.shape)*100:.2f}% of points below 0 Re")
    ax1.set_xlabel('Each cell altitude in Re')
    ax1.set_ylabel('bin count')

    ax1p2 = fig.add_subplot(gs[2,4:])
    alt_min_KM = r2alt(alt_min_RE, Re_km)
    counts, bins, _ = ax1p2.hist(r2alt(rs, Re_km).flatten(), bins=200)
    ax1p2.vlines(alt_min_KM, 0, max(counts)*1.1, linestyle = '--', alpha=.7, color='k')
    ax1p2.set_xlim(-100, 1000)
    ax1p2.set_xlabel('altitude in km')

    another_hist_ax = fig.add_subplot(gs[3:5, 4:])
    another_hist_ax.hist(np.rad2deg(ts.flatten()), bins=90)
    another_hist_ax.set_xlabel('Magnetic Latitude (deg)')

    ax2 = fig.add_subplot(gs[5:,:])
    for x,y in zip(np.rad2deg(ts), r2alt(rs, Re_km)):
        if abs_bot:
            ax2.scatter(np.abs(x),y)
        else:
            ax2.scatter(x,y)
            
    ax2.hlines(100, 0 if abs_bot else -90, 90, color='k', alpha=.8)
    ax2.set_ylim(0,1500)
    ax2.set_xlabel('Magnetic Latitude (deg)')
    ax2.set_ylabel('Altitude (km)')

    plt.tight_layout()

    return fig


def generate_sym_quadtree(nBlocks):
    """
    Makes the latitude portion of the quadtree
    input: nBlocks
    outputs:
        origins (normed y-coordinate of lower-left)
        extent (size_up_norm)
    """
    origins = np.linspace(-0.5, 0.5, num=nBlocks, endpoint=False)
    extent = 1/nBlocks

    return origins, extent



# ------------------------------------------------------------------------
# Main code is here:
# ------------------------------------------------------------------------
if __name__ == "__main__":

    alt_maxRE = alt2r(max_alt, Re_KM)
    alt_minRE = alt2r(min_alt, Re_KM)
    
    nLatsPerBlock = nLatsPerBlock_in + nGCs*2
    nAltsPerBlock = nAltsPerBlock_in + nGCs*2

    origins, extent = generate_sym_quadtree(nBlocks)

    qs, ps = main(alt_minRE, alt_maxRE, min_blat, max_blat, origins, extent, nLatsPerBlock, nAltsPerBlock)

    # if we want r&theta now:
    # rs = qp_solve(qs, ps)
    # ts = rq2t(rs, qs)

    # Otherwise the plotting function does it:

    fig = make_plot(qs, ps, alt_minRE, Re_KM, abs_bot=False)

    plt.show()

