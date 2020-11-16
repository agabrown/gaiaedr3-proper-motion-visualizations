"""
Plot an all-sky average proper motion map, using statistics downloaded from the Gaia archive with a query similar to the
following:

select
  gaia_healpix_index(5, source_id) as healpix_5,
  avg(pmra) as avg_pmra,
  avg(pmdec) as avg_pmdec
from gaiaedr3.gaia_source
where parallax_over_error>=10
and parallax*parallax - 2*parallax - parallax_error*parallax_error < -1
group by healpix_5

Anthony Brown Oct 2020 - Nov 2020
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, rc_file_defaults
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize, LogNorm
from matplotlib.patches import ArrowStyle
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.gridspec import GridSpec
import argparse

import cartopy.crs as ccrs

from astropy.coordinates import ICRS, Galactic
from astropy.table import Table
import astropy.units as u
import astropy_healpix.healpy as hp

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
    infile = './data/'+args['inputFile']
    basename = 'PMmap-' + args['inputFile'].split('.')[0]

    defaultProj = ccrs.PlateCarree()
    skyProj = ccrs.Mollweide()

    backgr = plt.imread('../star-trail-animation/sky-images/' +
            'ESA_Gaia_DR2_AllSky_Brightness_Colour_Cartesian_2000x1000.png')

    nside = hp.order2nside(args['hplevel'])
    hpcol = 'healpix_{0}'.format(args['hplevel'])
    edr3data = Table.read(infile)

    alpha, delta = hp.pix2ang(nside, edr3data[hpcol], lonlat=True, nest=True)
    pmra = edr3data['avg_pmra']
    pmdec = edr3data['avg_pmdec']

    icrs = ICRS(ra=alpha*u.degree, dec=delta*u.degree, pm_ra_cosdec=pmra*u.mas/u.yr, pm_dec=pmdec*u.mas/u.yr)
    galactic = icrs.transform_to(Galactic)

    fig=plt.figure(figsize=(16,9), dpi=120, frameon=False, tight_layout={'pad':0.01})
    gs = GridSpec(1, 1, figure=fig)
    ax = fig.add_subplot(gs[0,0], projection=skyProj)

    ax.imshow(np.fliplr(backgr), transform=defaultProj, zorder=-1, origin='upper')
    pmtot = (galactic.pm_l_cosb.value ** 2 + galactic.pm_b.value ** 2) ** 0.5
    pmcmap = cm.summer
    if args['quiver']:
        vscale = np.median(pmtot)/10
        ax.quiver(galactic.l.value, galactic.b.value, galactic.pm_l_cosb.value, galactic.pm_b.value,
                transform=defaultProj, angles='xy', scale=vscale, scale_units='dots', color='w',
                headwidth=1, headlength=3, headaxislength=2.5)
    else:
        ax.streamplot(galactic.l.value, galactic.b.value, galactic.pm_l_cosb.value, galactic.pm_b.value,
                transform=defaultProj, linewidth=1.5, density=2, color='w', maxlength=0.5, arrowsize=1,
                arrowstyle=ArrowStyle.Fancy(head_length=1.0, head_width=.4, tail_width=.4))
    ax.invert_xaxis()

    if args['pdfOutput']:
        plt.savefig(basename+'.pdf')
    elif args['pngOutput']:
        plt.savefig(basename+'.png')
    else:
        plt.show()

def parseCommandLineArguments():
    """
    Set up command line parsing.
    """
    parser = argparse.ArgumentParser("Produce all-sky proper motion map.")
    parser.add_argument('inputFile', type=str, help="""VOT file with proper motion stats by Healpix.""")
    parser.add_argument('hplevel', type=int, nargs='?', default=4, help="""Healpix level of input table.""")
    parser.add_argument('--vectors', action="store_true", dest="quiver", help="Plot vectors instead of streamlines")
    parser.add_argument("-p", action="store_true", dest="pdfOutput", help="Make PDF plot")
    parser.add_argument("-b", action="store_true", dest="pngOutput", help="Make PNG plot")
    args = vars(parser.parse_args())
    return args

if __name__ in ('__main__'):
    args=parseCommandLineArguments()
    make_plot(args)
