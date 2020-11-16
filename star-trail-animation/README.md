# star-trail-animation

Python code to reproduce the visualization of proper motions in Gaia EDR3 with a video showing as trails the
displacements of stars on the sky, over a fixed time interval.

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

This is done by running `star-trails-animation.py` from the command line. A few examples follow.

Usage information
```
python star-trails-animation.py -h
```
Small set of frames with a relatively small number of stars for testing: star positions are predicted for 0.4 million
years into the future and the star trails are shown over 80000 years each (fraction 10/50 of 1.6 million years). The end
frame shows the star trails over 400000 years. The frames are intended for video at Full HD resolution.
```
python star-trails-animation.py StarTrailVideoSample.fits --exposure 4 --exposure_endframe 4 --nstars_max=2000 --max_alpha=0.4 --num_epochs=50 --max_trail_epochs 10
```
Reproduce the animation made for the Gaia EDR3 release: star positions are predicted for 1.6 million years into the
future and the star trails are shown over 80000 years each (fraction 20/400 of 1.6 million years). The end frame shows
the star trails over 400000 years. The frames are intended for video at 4K UHD resolution.
```
python star-trails-animation.py StarTrailVideoSample.fits --exposure 16 --exposure_endframe 4 --nstars_max=40000 --max_alpha=0.4 --num_epochs=400 --max_plot_epochs 20 -4k
```

## Generating the video
