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

## Generating the start and end frames

## Generating the animation frames
