"""
Produce an animation of the Gaia DR3 sky overlaid with the star trails from a hypothetical long exposure. The trails are based on Gaia DR3 proper motions and radial velocities. After an idea originally form Stefan Jordan.
Anthony Brown Nov 2020 - Jan 2023
"""
from multiprocessing import Pool, cpu_count

import numpy as np
import matplotlib.pyplot as plt
import argparse

import cartopy.crs as ccrs
from astropy.table import Table

from pygaia.astrometry.coordinates import (
    EpochPropagation,
    CoordinateTransformation,
    Transformations,
)


def init(args):
    """
    Initialize the arrays of (l,b) values for the star trails, the background image, and other configuration variables.

    Parameters
    ----------

    args: array-like
        Command line arguments

    Returns
    -------

    Dictionary with the requested configuration. The l and b coordinates per source and per epoch for the animation
    frames and the end frame, the array for scaling the brightness of the stars and star trails, the background image,
    maximum number of epochs to plot, the image DPI, the folder in which to store the images, and the sky projection
    information.

    Example:
    config = init(args)
    """
    infile = args["inputFile"]
    delta_epoch = args["exposure"] * 1.0e5
    delta_epoch_endframe = args["exposure_endframe"] * 1.0e5
    n_epochs = args["num_epochs"]
    n_epochs_endframe = args["num_epochs_endframe"]
    t0 = 2016.0
    epp = EpochPropagation()
    epochs = np.linspace(t0, t0 + delta_epoch, n_epochs)
    epochs_endframe = np.linspace(t0, t0 + delta_epoch_endframe, n_epochs_endframe)
    max_trail_epochs = args["max_trail_epochs"]
    if args["highres"]:
        fdpi = 240
        dr3sky = plt.imread("./sky-images/GaiaSky-colour-4k.png")
        imfolder = "images-4k"
    else:
        fdpi = 120
        dr3sky = plt.imread("./sky-images/GaiaSky-colour-2k.png")
        imfolder = "images-2k"

    data = Table.read("./data/" + infile, format="fits")

    nsample = args["nstars_max"]
    print(
        f"A random selection of {nsample} out of {data['ra'].size} stars will be plotted"
    )

    rng = np.random.default_rng(args["rngseed"])
    raninds = rng.choice(np.arange(data["ra"].size), nsample)

    ra0 = np.deg2rad(data["ra"][raninds])
    dec0 = np.deg2rad(data["dec"][raninds])
    plx0 = data["parallax"][raninds]
    pmra0 = data["pmra"][raninds]
    pmdec0 = data["pmdec"][raninds]
    vrad0 = data["radial_velocity"][raninds]
    gmag = data["phot_g_mean_mag"][raninds]

    ct = CoordinateTransformation(Transformations.ICRS2GAL)

    l1 = np.zeros((nsample, n_epochs))
    b1 = np.zeros((nsample, n_epochs))
    l1ef = np.zeros((nsample, n_epochs_endframe))
    b1ef = np.zeros((nsample, n_epochs_endframe))

    for i, t1 in zip(range(n_epochs), epochs):
        ra1, dec1 = epp.propagate_pos(ra0, dec0, plx0, pmra0, pmdec0, vrad0, t0, t1)
        l1[:, i], b1[:, i] = ct.transform_sky_coordinates(ra1, dec1)
    l1 = np.rad2deg(l1)
    b1 = np.rad2deg(b1)

    for i, t1 in zip(range(n_epochs_endframe), epochs_endframe):
        ra1, dec1 = epp.propagate_pos(ra0, dec0, plx0, pmra0, pmdec0, vrad0, t0, t1)
        l1ef[:, i], b1ef[:, i] = ct.transform_sky_coordinates(ra1, dec1)
    l1ef = np.rad2deg(l1ef)
    b1ef = np.rad2deg(b1ef)

    magrange = gmag.max() - gmag.min()
    magscaling_stars = 0.2 + 0.8 * (gmag.max() - gmag) / magrange
    min_alpha = 0.1
    magscaling_trails = (
        min_alpha + (args["max_alpha"] - min_alpha) * (gmag.max() - gmag) / magrange
    )

    return {
        "l1": l1,
        "b1": b1,
        "l1ef": l1ef,
        "b1ef": b1ef,
        "magscaling_stars": magscaling_stars,
        "magscaling_trails": magscaling_trails,
        "backgr": dr3sky,
        "max_trail_epochs": max_trail_epochs,
        "dpi": fdpi,
        "imfolder": imfolder,
        "default_projection": ccrs.PlateCarree(),
        "sky_projection": ccrs.Mollweide(),
    }


def make_mwpanorama_frame(config):
    """
    Produce a frame containing only the Milky Way panorama.

    Parameters
    ----------

    config: dictionary
        Configuration of the animation

    Returns
    -------

    Nothing
    """
    framename = config["imfolder"] + "/panoramaframe.png"

    fig = plt.figure(
        figsize=(16, 9), dpi=config["dpi"], frameon=False, tight_layout={"pad": 0.01}
    )

    ax = fig.add_subplot(projection=config["sky_projection"])
    ax.imshow(
        np.fliplr(config["backgr"]),
        transform=config["default_projection"],
        zorder=-1,
        origin="upper",
    )
    ax.set_global()
    ax.invert_xaxis()

    plt.savefig(framename)
    plt.close(fig)


def make_start_frame(config):
    """
    Produce the start frame of the animation.

    Parameters
    ----------

    config: dictionary
        Configuration of the animation

    Returns
    -------

    Nothing
    """
    framename = config["imfolder"] + "/startframe.png"
    nsample = config["l1"].shape[0]

    fig = plt.figure(
        figsize=(16, 9), dpi=config["dpi"], frameon=False, tight_layout={"pad": 0.01}
    )

    ax = fig.add_subplot(projection=config["sky_projection"])
    ax.imshow(
        np.fliplr(config["backgr"]),
        transform=config["default_projection"],
        zorder=-1,
        origin="upper",
    )
    for i in range(nsample):
        ax.plot(
            config["l1"][i, 0],
            config["b1"][i, 0],
            "o",
            ms=1.0,
            color="w",
            transform=config["default_projection"],
            alpha=config["magscaling_stars"][i],
        )
    ax.set_global()
    ax.invert_xaxis()

    plt.savefig(framename)
    plt.close(fig)


def make_end_frame(config):
    """
    Produce the end frame of the animation.

    Parameters
    ----------

    config: dictionary
        Configuration of the animation

    Returns
    -------

    Nothing
    """
    framename = config["imfolder"] + "/endframe.png"
    nsample = config["l1ef"].shape[0]

    fig = plt.figure(
        figsize=(16, 9), dpi=config["dpi"], frameon=False, tight_layout={"pad": 0.01}
    )
    trailwidth = 0.5

    ax = fig.add_subplot(projection=config["sky_projection"])
    ax.imshow(
        np.fliplr(config["backgr"]),
        transform=config["default_projection"],
        zorder=-1,
        origin="upper",
    )
    for i in range(nsample):
        galon = config["l1ef"][i, :]
        galat = config["b1ef"][i, :]
        galon[galon > 180] = galon[galon > 180] - 360.0
        diffs = galon[1:] - galon[0:-1]
        if np.all(diffs > 0) or np.all(diffs < 0):
            ax.plot(
                galon,
                galat,
                c="w",
                lw=trailwidth,
                alpha=config["magscaling_trails"][i],
                transform=config["default_projection"],
            )
        else:
            indices = galon >= 0.0
            if np.any(indices):
                xplot = galon[indices]
                yplot = galat[indices]
                ax.plot(
                    xplot,
                    yplot,
                    c="w",
                    lw=trailwidth,
                    alpha=config["magscaling_trails"][i],
                    transform=config["default_projection"],
                )
            indices = galon < 0.0
            if np.any(indices):
                xplot = galon[indices]
                yplot = galat[indices]
                ax.plot(
                    xplot,
                    yplot,
                    c="w",
                    lw=trailwidth,
                    alpha=config["magscaling_trails"][i],
                    transform=config["default_projection"],
                )
    ax.set_global()
    ax.invert_xaxis()

    plt.savefig(framename)
    plt.close(fig)


def make_frame(config, end_epoch_num):
    """
    Make a frame for the star trail animation.

    Parameters
    ----------

    config: dictionary
        Configuration of the animation
    end_epoch_num: int
        Last epoch to consider in drawing star trail.

    Returns
    -------

    Nothing
    """
    nsample = config["l1"].shape[0]
    start_epoch = max(0, end_epoch_num - config["max_trail_epochs"])

    fig = plt.figure(
        figsize=(16, 9), dpi=config["dpi"], frameon=False, tight_layout={"pad": 0.01}
    )

    framename = config["imfolder"] + "/frame{0:04d}.png"
    trailwidth = 1.0
    star_fade_frames = 25

    ax = fig.add_subplot(projection=config["sky_projection"])
    ax.imshow(
        np.fliplr(config["backgr"]),
        transform=config["default_projection"],
        zorder=-1,
        origin="upper",
    )
    if end_epoch_num - 2 <= star_fade_frames:
        fade_factor = (star_fade_frames - (end_epoch_num - 2)) / star_fade_frames
        for i in range(nsample):
            ax.plot(
                config["l1"][i, 0],
                config["b1"][i, 0],
                "o",
                ms=1.0,
                color="w",
                transform=config["default_projection"],
                alpha=config["magscaling_stars"][i] * fade_factor,
            )
    for i in range(nsample):
        galon = config["l1"][i, start_epoch : end_epoch_num + 1]
        galat = config["b1"][i, start_epoch : end_epoch_num + 1]
        galon[galon > 180] = galon[galon > 180] - 360.0
        diffs = galon[1:] - galon[0:-1]
        if np.all(diffs > 0) or np.all(diffs < 0):
            ax.plot(
                galon,
                galat,
                c="w",
                lw=trailwidth,
                alpha=config["magscaling_trails"][i],
                transform=config["default_projection"],
            )
        else:
            indices = galon >= 0.0
            if np.any(indices):
                xplot = galon[indices]
                yplot = galat[indices]
                ax.plot(
                    xplot,
                    yplot,
                    c="w",
                    lw=trailwidth,
                    alpha=config["magscaling_trails"][i],
                    transform=config["default_projection"],
                )
            indices = galon < 0.0
            if np.any(indices):
                xplot = galon[indices]
                yplot = galat[indices]
                ax.plot(
                    xplot,
                    yplot,
                    c="w",
                    lw=trailwidth,
                    alpha=config["magscaling_trails"][i],
                    transform=config["default_projection"],
                )
    ax.set_global()
    ax.invert_xaxis()
    plt.savefig(framename.format(end_epoch_num - 2))
    plt.close(fig)


def parse_command_line_arguments():
    """
    Set up command line parsing.
    """
    parser = argparse.ArgumentParser("Produce all-sky star trail map.")
    parser.add_argument(
        "inputFile",
        type=str,
        help="""File with EDR3 astrometry and radial velocity data (located in ./data/)""",
    )
    parser.add_argument(
        "--exposure",
        type=float,
        nargs="?",
        default=4,
        help="""Exposure time in units of 100 kyr (default 4)""",
    )
    parser.add_argument(
        "--exposure_endframe",
        type=float,
        nargs="?",
        default=4,
        help="""Exposure time for final frame in units of 100 kyr (default 4)""",
    )
    parser.add_argument(
        "--nstars_max",
        type=int,
        nargs="?",
        default=2000,
        help="""Max number of stars to plot (default 2000)""",
    )
    parser.add_argument(
        "--max_alpha",
        type=float,
        nargs="?",
        default=0.4,
        help="""Max opacity of star trails (>0.1, default 0.4)""",
    )
    parser.add_argument(
        "--num_epochs",
        type=int,
        nargs="?",
        default=100,
        help="""Number of time steps for (default 100)""",
    )
    parser.add_argument(
        "--num_epochs_endframe",
        type=int,
        nargs="?",
        default=100,
        help="""Number of time steps for end
            frame trails (default 100)""",
    )
    parser.add_argument(
        "--max_trail_epochs",
        type=int,
        nargs="?",
        default=np.iinfo(np.int32).max,
        help="""Maximum number
            of time steps to plot for one trail (default inf)""",
    )
    parser.add_argument(
        "--max_cores",
        type=int,
        nargs="?",
        default=np.iinfo(np.int32).max,
        help="""Maximum number of CPU
            cores to use (default inf)""",
    )
    parser.add_argument(
        "--rngseed",
        type=int,
        nargs="?",
        default=53949896,
        help="""Random number generator seed (default 53949896)""",
    )
    parser.add_argument(
        "-4k",
        action="store_true",
        dest="highres",
        help="""Generate frames for a 4K UHD video (default is Full HD)""",
        default=False,
    )
    parser.add_argument(
        "--start_end_only",
        action="store_true",
        dest="startendonly",
        help="""Generate start and end frames only""",
        default=False,
    )
    args = vars(parser.parse_args())
    return args


if __name__ in "__main__":
    cmdargs = parse_command_line_arguments()

    the_config = init(cmdargs)
    make_mwpanorama_frame(the_config)
    make_start_frame(the_config)
    make_end_frame(the_config)
    if cmdargs["startendonly"]:
        print("Only generated start and end frames!")
        exit()

    num_epochs = the_config["l1"].shape[1]
    max_cores = min(cpu_count(), cmdargs["max_cores"])
    print(f"Using {max_cores} cores...")

    with Pool(processes=max_cores) as pool:
        pool.starmap(make_frame, [(the_config, n) for n in range(2, num_epochs)])
