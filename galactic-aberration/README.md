# galactic-aberration

Visualize the apparent proper motions of quasars caused by the acceleration of the solar system barycentre. See the
 paper by [Gaia Collaboration, Klioner, et al. (2020)](https://www.cosmos.esa.int/web/gaia/edr3-papers).

## Dependencies

* [numpy](https://numpy.org/)
* [matplotlib](https://matplotlib.org/)
* [Astropy](https://www.astropy.org/)
* [astropy-healpix](https://astropy-healpix.readthedocs.io/)
* [Cartopy](https://scitools.org.uk/cartopy/docs/latest/)

## Notebooks

The notebook [GalacticAberration_QSOProperMotions.ipynb](./GalacticAberration_QSOProperMotions.ipynb) shows how to
 calculate the apparent proper motions of the QSOs and visualize them. This is not done for actual quasars but for a
  regular grid of points on the sky.
  
## Python command line code

The programme [`qso-proper-motions-galactic-aberration.py`](./qso-proper-motions-galactic-aberration.py) provides a
 command line interface to make the visualization, which offers a bit more flexibility in plot configurations.