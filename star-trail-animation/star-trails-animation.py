"""
Produce an animation of the Gaia DR2 sky overlaid with the star trails from a hypothetical long exposure. The trails are
based on Gaia EDR3 proper motions and radial velocities. After an idea originally form Stefan Jordan.

Anthony Brown Nov 2020 - Nov 2020
"""
from multiprocessing import Pool, cpu_count

import numpy as np
import matplotlib.pyplot as plt
import argparse

import cartopy.crs as ccrs
from astropy.table import Table

from pygaia.astrometry.coordinates import EpochPropagation, CoordinateTransformation, Transformations


def init(args, default_proj, sky_proj):
    """
    Initialize the arrays of (l,b) values for the star trails, the background image, and the start frame of the
    animation.

    Parameters
    ----------

    args: array-like
        Command line arguments
    default_proj: object
        Cartopy projection that is the default for the calculations of the star positions and trails.
    sky_proj: object
        Cartopy projection that is used for the sky map

    Returns
    -------

    Arrays with l and b coordinates per source and per epoch, the array for scaling the brightness of the stars and star
    trails, the background image, maximum number of epochs to plot, the image DPI, and the folder in which to store the
    images.

    Example:
    l, b, magscaling_stars, magscaling_trails, backgr_image, max_plot_epochs, fdpi, imfolder = init(args)
    """
    infile = args['inputFile']
    delta_epoch = args['exposure']*1.0e+5
    n_epochs = args['num_epochs']
    t0=2016.0
    epp = EpochPropagation()
    epochs = np.linspace(t0, t0+delta_epoch, n_epochs)
    max_plot_epochs = args['max_plot_epochs']
    if args['highres']:
        fdpi=240
        dr2sky = plt.imread('./sky-images/ESA_Gaia_DR2_AllSky_Brightness_Colour_Cartesian_4000x2000.png')
        imfolder = "images-4k"
    else:
        fdpi=120
        dr2sky = plt.imread('./sky-images/ESA_Gaia_DR2_AllSky_Brightness_Colour_Cartesian_2000x1000.png')
        imfolder = "images-2k"


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

    magrange = mag.max()-mag.min()
    magscaling_stars = 0.2+0.8*(mag.max()-mag)/magrange
    min_alpha=0.1
    magscaling_trails = min_alpha + (args['max_alpha']-min_alpha)*(mag.max()-mag)/magrange

    framename = imfolder+"/startframe.png"

    fig=plt.figure(figsize=(16,9), dpi=fdpi, frameon=False, tight_layout={'pad':0.01})

    ax0 = fig.add_subplot(projection=sky_proj)
    ax0.imshow(np.fliplr(dr2sky), transform=default_proj, zorder=-1, origin='upper')
    for i in range(nsample):
        ax0.plot(l1[i,0], b1[i,0], 'o', ms=1.0, color='w', transform=default_proj, alpha=magscaling_stars[i])
    ax0.set_global()
    ax0.invert_xaxis()

    plt.savefig(framename)
    plt.close(fig)
    return l1, b1, magscaling_stars, magscaling_trails, dr2sky, max_plot_epochs, fdpi, imfolder


def make_frame(l1, b1, magscaling_stars, magscaling_trails, backgr, defaultProj, skyProj, end_epoch_num, max_epochs,
        fdpi, imfolder):
    """
    Make a frame for the star trail animation.

    Parameters
    ----------

    l1: array (nsources, nepochs) 
        Galactic longitude coordinate values (degree).
    b1: array (nsources, nepochs) 
        Galactic latitude coordinate values (degree).
    magscaling_stars: array (nsources)
        Brightness (transparency) scaling of star symbols
    magscaling_trails: array (nsources)
        Brightness (transparency) scaling of star trails.
    backgr: matplotlib image
        Background image for video in default (Cartesian) projection.
    defaultProj: object
        Cartopy projection that is the default (Cartesian) for the calculations of the star positions and trails.
    sky_proj: object
        Cartopy projection that is used for the sky map.
    end_epoch_num: int
        Last epoch to consider in drawing star trail.
    max_epochs: int
        Maximum number of epochs to plot, counting back from end_epoch_num.
    fdpi: int
        Resolution parameter in dots per inch.
    imfolder: str
        Name of folder in which to store images

    Returns
    -------

    Nothing
    """
    nsample = l1.shape[0]
    start_epoch = max(0,end_epoch_num-max_epochs)

    fig=plt.figure(figsize=(16,9), dpi=fdpi, frameon=False, tight_layout={'pad':0.01})

    framename = imfolder+"/frame{0:04d}.png"
    trailwidth = 1.0
    star_fade_frames = 25

    ax1 = fig.add_subplot(projection=skyProj)
    ax1.imshow(np.fliplr(backgr), transform=defaultProj, zorder=-1, origin='upper')
    if end_epoch_num-2 <= star_fade_frames:
        fade_factor = (star_fade_frames - (end_epoch_num-2))/star_fade_frames
        for i in range(nsample):
            ax1.plot(l1[i,0], b1[i,0], 'o', ms=1.0, color='w', transform=defaultProj,
                    alpha=magscaling_stars[i]*fade_factor)
    for i in range(nsample):
        l = l1[i,start_epoch:end_epoch_num+1]
        b = b1[i,start_epoch:end_epoch_num+1]
        l[l>180] = l[l>180]-360.0
        diffs = l[1:]-l[0:-1]
        if np.all(diffs>0) or np.all(diffs<0):
            ax1.plot(l, b, c='w', lw=trailwidth, alpha=magscaling_trails[i], transform=defaultProj)
        else:
            indices=(l>=0.0)
            if np.any(indices):
                xplot=l[indices]
                yplot=b[indices]
                ax1.plot(xplot, yplot, c='w', lw=trailwidth, alpha=magscaling_trails[i], transform=defaultProj)
            indices=(l<0.0)
            if np.any(indices):
                xplot=l[indices]
                yplot=b[indices]
                ax1.plot(xplot, yplot, c='w', lw=trailwidth, alpha=magscaling_trails[i], transform=defaultProj)
    ax1.set_global()
    ax1.invert_xaxis()
    plt.savefig(framename.format(end_epoch_num-2))
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
    parser.add_argument('--num_epochs', type=int, nargs='?', default=100, help="""Number of time steps (default 100)""")
    parser.add_argument('--max_plot_epochs', type=int, nargs='?', default=np.iinfo(np.int).max, help="""Maximum number
            of time steps to plot for one trail (default inf)""")
    parser.add_argument('--max_cores', type=int, nargs='?', default=np.iinfo(np.int).max, help="""Maximum number of CPU
            cores to use (default inf)""")
    parser.add_argument('-4k', action="store_true", dest="highres", help="""Generate frames for a 4K video""", default=False)
    args = vars(parser.parse_args())
    return args

if __name__ in ('__main__'):
    args=parseCommandLineArguments()

    defaultProj = ccrs.PlateCarree()
    skyProj = ccrs.Mollweide()

    l1, b1, magscaling_stars, magscaling_trails, dr2sky, max_epochs, fdpi, imfolder = init(args, defaultProj, skyProj)
    n_epochs = l1.shape[1]

    max_cores = min(cpu_count(), args['max_cores'])
    print(f'Using {max_cores} cores...')

    with Pool(processes = max_cores) as pool:
        pool.starmap(make_frame, [(l1, b1, magscaling_stars, magscaling_trails, dr2sky, defaultProj, skyProj, n,
            max_epochs, fdpi, imfolder) for n in range(2,n_epochs)])
