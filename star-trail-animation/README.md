# star-trail-animation

Python code to reproduce the visualization of proper motions in Gaia (E)DR3. Produces a video showing the displacements of
stars on the sky over a certain time interval T, starting from the Gaia (E)DR3 reference epoch of J2016.0. The
displacements are animated as short trails, which illustrate the displacement over a short interval t&lt;T, and which
flow along the full trajectory of the star on the sky, as defined by T.

> **2022.12.30**
> The code and data were updated to make use of Gaia DR3 (with the new radial velocities) instead of Gaia EDR3.

## Credits

Original idea by Stefan Jordan. Python code and video by Anthony Brown. Ideas for improvement contributed by Tineke
Roegiers, Xavier Luri, Eduard Masana, Timo Prusti.

## Dependencies

* [numpy](https://numpy.org/)
* [matplotlib](https://matplotlib.org/)
* [Astropy](https://www.astropy.org/)
* [Cartopy](https://scitools.org.uk/cartopy/docs/latest/)
* [PyGaia](https://github.com/agabrown/PyGaia)
* [ffmpeg](https://ffmpeg.org/)
* [GNU Bash](https://www.gnu.org/software/bash/)

## Input data needed

The animation can be done for any sample of stars from Gaia DR3. The catalogue fields that should be present are:
```
ra
dec
parallax
pmra
pmdec
radial_velocity
phot_g_mean_mag
```

## Jupyter notebook

The [notebook](StarTrailsOnSky.ipynb) explains the basics of rendering the frames in the animation with [matplotlib](https://matplotlib.org/).

## Generating the animation frames

This is done by running [`star-trails-animation.py`](star-trails-animation.py) from the command line. A few examples follow.

#### Usage information

```
python star-trails-animation.py -h
```
#### Small set of frames with a relatively small number of stars for testing.

Star positions are predicted for 0.4 million years into the future and the star trails are shown over 80000 
years each (fraction 10/50 of 1.6 million years). The end frame shows the star trails over 400000 years. The 
frames are intended for video at Full HD resolution.

```
python star-trails-animation.py StarTrailVideoSample.fits --exposure 4 --exposure_endframe 4 --nstars_max=2000 --max_alpha=0.4 --num_epochs=50 --max_trail_epochs 10
```

#### Reproduce the animation made for the Gaia EDR3 release

Star positions are predicted for 1.6 million years into the future and the star trails are shown over 80000 
years each (fraction 20/400 of 1.6 million years). The end frame shows the star trails over 400000 years. 
The frames are intended for video at 4K UHD resolution.

```
python star-trails-animation.py StarTrailVideoSample.fits --exposure 16 --exposure_endframe 4 --nstars_max=40000 --max_alpha=0.4 --num_epochs=400 --max_trail_epochs 20 -4k
```

The input data was retrieved from the Gaia archive with the following query:
```
select ra, dec, parallax, pmra, pmdec, radial_velocity, phot_g_mean_mag, bp_rp, grvs_mag
from gaiadr3.gaia_source
where astrometric_params_solved != 3
and parallax_over_error>10 and parallax>10
and ruwe<1.4
and radial_velocity is not null
and rv_expected_sig_to_noise > 5
and bp_rp is not null
and grvs_mag is not null
```
Data quality criteria select the best astrometry and radial velocities.

The background images are expected to be stores in the [sky-images](./sky-images/) folder. The
names are hard-coded in the Python programme as `GaiaSky-colour-2k.png` and `GaiaSky-colour-4k.png`.

## Generating the video

The video is generated with ffmpeg (of course any other software that can turn a collection of frames into a video will
do), using the Bash script [`makevideo.sh`](makevideo.sh).

For Full HD video run, including title and explanation frames:
```
makevideo.sh -t
```

For 4K UHD video run, including title and explanation frames:
```
makevideo.sh -t -k
```

Videos without the title frames:
```
makevideo.sh [-k]
```

Videos in HAP format (produces `.mov` files):
```
makevideo.sh -a [-t] [-k]
```

Videos cropped to image format requested for projection in [Rabo Studio](https://zakelijk.forum.nl/nl/onze-ruimtes/rabo-studio) (3720x1104 pixels for center screen):
```
makevideo.sh -c -k [-a] [-t]
```
Without the `-k` option the video will be cropped to 1860x552 pixels.
