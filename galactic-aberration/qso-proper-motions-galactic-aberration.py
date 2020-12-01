"""
Visualize the apparent quasar proper motions caused by the acceleration of the solar sytem barycentre

Use the results from Gaia Collaboration, Klioner, et al. (2020) to visualize the apparent proper motions of quasars
caused by the acceleration of the solar system barycentre (also known as 'Galactic aberration'). Use equation (4) from
the paper to calculate the apparent proper motions.

Anthony Brown Nov 2020 - Dec 2020
"""

import argparse

import astropy.units as u
import astropy_healpix.healpy as hp
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import Galactic
from astropy_healpix import HEALPix
from matplotlib.gridspec import GridSpec
from matplotlib.patches import ArrowStyle


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
    basename = 'PMmap-qso-galactic-aberration'

    gx = 5.04
    gy = -0.10
    gz = -0.29

    if args['quiver']:
        hplevel = 3
    else:
        hplevel = 5
    nside = hp.order2nside(hplevel)
    npix = hp.nside2npix(nside)
    ahp = HEALPix(nside=nside, order='nested', frame=Galactic())

    hpindices = np.arange(npix)
    skycoords = ahp.healpix_to_skycoord(hpindices)

    pm_l_cosb = -gx * np.sin(skycoords.l.to(u.rad)) + gy * np.cos(skycoords.l.to(u.rad))
    pm_b = -gx * np.sin(skycoords.b.to(u.rad)) * np.cos(skycoords.l.to(u.rad)) \
           - gy * np.sin(skycoords.b.to(u.rad)) * np.sin(skycoords.l.to(u.rad)) \
           + gz * np.cos(skycoords.b.to(u.rad))
    pmtot = np.sqrt(pm_l_cosb ** 2 + pm_b ** 2)

    backgr = plt.imread('../star-trail-animation/sky-images/GaiaSky-colour-2k.png')

    default_proj = ccrs.PlateCarree()
    sky_proj = ccrs.Mollweide()

    fig = plt.figure(figsize=(16, 9), dpi=120, frameon=False, tight_layout={'pad': 0.01})
    gs = GridSpec(1, 1, figure=fig)
    ax = fig.add_subplot(gs[0, 0], projection=sky_proj)
    ax.imshow(np.fliplr(backgr), transform=default_proj, zorder=-1, origin='upper')
    veccolor = plt.cm.get_cmap('tab10').colors[9]
    linecolor = 'w'  # plt.cm.get_cmap('tab10').colors[9]

    if args['quiver']:
        vscale = np.median(pmtot) / 50
        ax.quiver(skycoords.l.value, skycoords.b.value, pm_l_cosb, pm_b, transform=default_proj, angles='xy',
                  scale=vscale, scale_units='dots', color=veccolor, headwidth=4, headlength=4, headaxislength=3.5)
    else:
        if args['colourstreams']:
            ax.streamplot(skycoords.l.value, skycoords.b.value, pm_l_cosb, pm_b, transform=default_proj, linewidth=2.0,
                          density=2, color=pmtot, cmap='viridis', maxlength=0.5, arrowsize=1,
                          arrowstyle=ArrowStyle.Fancy(head_length=1.0, head_width=.4, tail_width=.4))
        elif args['lwcode'] > 0:
            ax.streamplot(skycoords.l.value, skycoords.b.value, pm_l_cosb, pm_b, transform=default_proj,
                          linewidth=args['lwcode'] * pmtot / np.median(pmtot), density=2, color=linecolor,
                          maxlength=0.5,
                          arrowsize=1, arrowstyle=ArrowStyle.Fancy(head_length=1.0, head_width=.4, tail_width=.4))
        else:
            ax.streamplot(skycoords.l.value, skycoords.b.value, pm_l_cosb, pm_b, transform=default_proj, linewidth=1.5,
                          density=2.5, color=linecolor, maxlength=0.5, arrowsize=1,
                          arrowstyle=ArrowStyle.Fancy(head_length=1.0, head_width=.4, tail_width=.4))
    # ax.gridlines()
    ax.invert_xaxis()

    if args['pdfOutput']:
        plt.savefig(basename + '.pdf')
    elif args['pngOutput']:
        plt.savefig(basename + '.png')
    else:
        plt.show()


def parse_command_line_arguments():
    """
    Set up command line parsing.
    """
    parser = argparse.ArgumentParser("Produce all-sky proper motion map.")
    parser.add_argument('--vectors', action="store_true", dest="quiver", help="Plot vectors instead of streamlines")
    parser.add_argument('--colourcode', action='store_true', dest='colourstreams', help="""Plot streamlines colour coded
    by magnitude of proper motion""")
    parser.add_argument('--lwcode', type=float, default=0.0, help="""Plot streamlines with the width indicating the
    magnitude of proper motion. Scale the widths by the factor provided""")
    parser.add_argument("-p", action="store_true", dest="pdfOutput", help="Make PDF plot")
    parser.add_argument("-b", action="store_true", dest="pngOutput", help="Make PNG plot")
    args = vars(parser.parse_args())
    return args


if __name__ in '__main__':
    cmdargs = parse_command_line_arguments()
    make_plot(cmdargs)
