"""
Produce an image of the Gaia DR2 sky overlaid with the star trails from a hypothetical long exposure. The trails are
based on Gaia EDR3 proper motions and radial velocities. After an idea originally form Stefan Jordan.

Anthony Brown Nov 2020 - Nov 2020
"""
import numpy as np
import matplotlib.pyplot as plt
import argparse

import cartopy.crs as ccrs

from astropy.table import Table

from pygaia.astrometry.coordinates import EpochPropagation, CoordinateTransformation, Transformations


def make_plot(args):
    """
    Take the steps to make the plot.

    Parameters
    ----------

    args: array-like
        Command line arguments

    Returns
    -------

    Nothing
    """
    infile = args['inputFile']
    delta_epoch = args['exposure']*1.0e+5
    n_epochs = 100
    t0=2016.0
    epp = EpochPropagation()
    epochs = np.linspace(t0, t0+delta_epoch, n_epochs)

    if args['highres']:
        fdpi=240
        dr2sky = plt.imread('./sky-images/ESA_Gaia_DR2_AllSky_Brightness_Colour_Cartesian_4000x2000.png')
    else:
        fdpi=120
        dr2sky = plt.imread('./sky-images/ESA_Gaia_DR2_AllSky_Brightness_Colour_Cartesian_2000x1000.png')

    # Read input file (for now the Gaia EDR3 6D-gold sample from which a random selection of stars within 100 pc is
    # selected.
    data = Table.read(infile, format='fits')

    slice_dist = (data['parallax']>10) & (data['parallax_over_error']>=10)
    nsample = args['nstars_max']
    print(f"A random selection of {nsample} out of {data['ra'][slice_dist].size} selected stars will be plotted")

    rng = np.random.default_rng(53949896)
    raninds = rng.choice(np.arange(data['ra'][slice_dist].size), nsample)

    ra0 = np.deg2rad(data['ra'][slice_dist][raninds])
    dec0 = np.deg2rad(data['dec'][slice_dist][raninds])
    plx0 = data['parallax'][slice_dist][raninds]
    pmra0 = data['pmra'][slice_dist][raninds]
    pmdec0 = data['pmdec'][slice_dist][raninds]
    vrad0 = data['dr2_radial_velocity'][slice_dist][raninds]
    mag = data['phot_g_mean_mag'][slice_dist][raninds]

    ct = CoordinateTransformation(Transformations.ICRS2GAL)

    l1 = np.zeros((nsample, n_epochs))
    b1 = np.zeros((nsample, n_epochs))

    for i, t1 in zip(range(n_epochs), epochs):
        ra1, dec1 = epp.propagate_pos(ra0, dec0, plx0, pmra0, pmdec0, vrad0, t0, t1)
        l1[:,i], b1[:,i] = ct.transformSkyCoordinates(ra1, dec1)

    l1 = np.rad2deg(l1)
    b1 = np.rad2deg(b1)

    defaultProj = ccrs.PlateCarree()
    skyProj = ccrs.Mollweide()

    magrange = mag.max()-mag.min()
    magscaled_stars = 0.2+0.8*(mag.max()-mag)/magrange
    min_alpha=0.1
    magscaled_trails = min_alpha + (args['max_alpha']-min_alpha)*(mag.max()-mag)/magrange

    basename = f"star-trails-{nsample}-{delta_epoch/1000:.0f}kyr"

    fig=plt.figure(figsize=(16,9), dpi=fdpi, frameon=False, tight_layout={'pad':0.01})

    ax0 = fig.add_subplot(projection=skyProj)
    ax0.imshow(np.fliplr(dr2sky), transform=defaultProj, zorder=-1, origin='upper')
    for i in range(nsample):
        ax0.plot(l1[i,0], b1[i,0], 'o', ms=1.0, color='w', transform=defaultProj, alpha=magscaled_stars[i])
    ax0.set_global()
    ax0.invert_xaxis()

    if args['pdfOutput']:
        plt.savefig(basename+'-A.pdf')
    elif args['pngOutput']:
        plt.savefig(basename+'-A.png')
    plt.close(fig)

    fig=plt.figure(figsize=(16,9), dpi=fdpi, frameon=False, tight_layout={'pad':0.01})
    trailwidth=0.5

    ax1 = fig.add_subplot(projection=skyProj)
    ax1.imshow(np.fliplr(dr2sky), transform=defaultProj, zorder=-1, origin='upper')
    for i in range(nsample):
        l = l1[i,:]
        b = b1[i,:]
        l[l>180] = l[l>180]-360.0
        diffs = l[1:]-l[0:-1]
        if np.all(diffs>0) or np.all(diffs<0):
            ax1.plot(l, b, c='w', lw=trailwidth, alpha=magscaled_trails[i], transform=defaultProj)
        else:
            indices=(l>=0.0)
            if np.any(indices):
                xplot=l[indices]
                yplot=b[indices]
                ax1.plot(xplot, yplot, c='w', lw=trailwidth, alpha=magscaled_trails[i], transform=defaultProj)
            indices=(l<0.0)
            if np.any(indices):
                xplot=l[indices]
                yplot=b[indices]
                ax1.plot(xplot, yplot, c='w', lw=trailwidth, alpha=magscaled_trails[i], transform=defaultProj)

    ax1.set_global()
    ax1.invert_xaxis()

    if args['pdfOutput']:
        plt.savefig(basename+'-B.pdf')
    elif args['pngOutput']:
        plt.savefig(basename+'-B.png')
    plt.close(fig)


def parseCommandLineArguments():
    """
    Set up command line parsing.
    """
    parser = argparse.ArgumentParser("Produce all-sky star trail map.")
    parser.add_argument('inputFile', type=str, help="""File with EDR3 astrometry and radial velocity data""")
    parser.add_argument('--exposure', type=float, nargs='?', default=4, help="""Exposure time in units of 100 kyr (default 4)""")
    parser.add_argument('--nstars_max', type=int, nargs='?', default=2000, help="""Max number of stars to plot (default 2000)""")
    parser.add_argument('--max_alpha', type=float, nargs='?', default=0.8, help="""Max transparency of star trails (>0.1, default 0.8)""")
    parser.add_argument('-4k', action="store_true", dest="highres", help="""Generate frames for a 4K video""", default=False)
    parser.add_argument("-p", action="store_true", dest="pdfOutput", help="Make PDF plot")
    parser.add_argument("-b", action="store_true", dest="pngOutput", help="Make PNG plot", default=True)
    args = vars(parser.parse_args())
    return args

if __name__ in ('__main__'):
    args=parseCommandLineArguments()
    make_plot(args)
