# star-trail-animation

Python code to reproduce the visualization of proper motions in Gaia EDR3 with a video showing as trails the
displacements of stars on the sky, over a fixed time interval.

## Dependencies

* [numpy](https://numpy.org/)
* [matplotlib](https://matplotlib.org/)
* [Astropy](https://www.astropy.org/)
* [Cartopy](https://scitools.org.uk/cartopy/docs/latest/)
* [PyGaia](https://github.com/agabrown/PyGaia)

## Input data needed

The animation can be done for any sample of stars from Gaia EDR3. The catalogue fields that should be present are:
```
ra
dec
parallax
pmra
pmdec
dr2_radial_velocity
phot_g_mean_mag
```

## Jupyter notebook

The [notebook](StarTrailsOnSky.ipynb) explains the basics of rendering the frames in the animation with [matplotlib](https://matplotlib.org/).

## Generating the animation frames

This is done by running [`star-trails-animation.py`](start-trails-animation.py) from the command line. A few examples follow.

#### Usage information

```
python star-trails-animation.py -h
```
#### Small set of frames with a relatively small number of stars for testing.

Star positions are predicted for 0.4 million
years into the future and the star trails are shown over 80000 years each (fraction 10/50 of 1.6 million years). The end
frame shows the star trails over 400000 years. The frames are intended for video at Full HD resolution.

```
python star-trails-animation.py StarTrailVideoSample.fits --exposure 4 --exposure_endframe 4 --nstars_max=2000 --max_alpha=0.4 --num_epochs=50 --max_trail_epochs 10
```

#### Reproduce the animation made for the Gaia EDR3 release

Star positions are predicted for 1.6 million years into the
future and the star trails are shown over 80000 years each (fraction 20/400 of 1.6 million years). The end frame shows
the star trails over 400000 years. The frames are intended for video at 4K UHD resolution.

```
python star-trails-animation.py StarTrailVideoSample.fits --exposure 16 --exposure_endframe 4 --nstars_max=40000 --max_alpha=0.4 --num_epochs=400 --max_plot_epochs 20 -4k
```

The input data was retrieved from the Gaia archive with the following query:
```
select ra, dec, parallax, pmra, pmdec, dr2_radial_velocity, phot_g_mean_mag
from gaiaedr3.gaia_source
where astrometric_params_solved != 3
and parallax_over_error>10 and parallax>10
and ruwe<1.4
and dr2_radial_velocity is not null
and dr2_rv_nb_transits > 3
and bp_rp is not null
```

## Generating the video

The video is generated with ffmpeg (of course any other software that can turn a collection of frames into a video will
do), using the Bash script [`makevideo.sh`](makevideo.sh).

For Full HD video run:
```
makevideo.sh
```

For 4K UHD video run:
```
makevideo.sh -k
```
