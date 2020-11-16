# all-sky-average-proper-motions

Show the average proper motions of stars projected as streamlines on the sky. The goal is to illustrate the effects of
galactic differential rotation and the sun's peculiar motion as seen in the data.

## Dependencies

* [numpy](https://numpy.org/)
* [matplotlib](https://matplotlib.org/)
* [Astropy](https://www.astropy.org/)
* [astropy-healpix](https://astropy-healpix.readthedocs.io/)
* [Cartopy](https://scitools.org.uk/cartopy/docs/latest/)
* [healpy](https://github.com/healpy/healpy)

## Extracting the mean proper motions from the Gaia archive

The mean proper motions can be calculated on the fly as part of an ADQL query in the Gaia archive. The following example
query calculates the mean proper motions per healpix (level 5) for stars located in a narrow slice around 1 mas in
parallax (nominally 1000 pc in distance).

```
select
  gaia_healpix_index(5, source_id) as healpix_5,
  avg(pmra) as avg_pmra,
  avg(pmdec) as avg_pmdec
from gaiaedr3.gaia_source
where parallax_over_error>=10
and parallax*parallax - 2*parallax - parallax_error*parallax_error < -1
group by healpix_5
```

## Notebooks

The [AllSkyProperMotionMap](AllSkyProperMotionMap.ipynb) notebook shows how to create the streamline plot visualizing the
average proper motions.

The [ExpectedKinematicSkyMaps](ExpectedKinematicSkyMaps.ipynb) notebook presents a very simplistic model for the Milky
Way disk kinematics that can predict the gross features seen in the proper motion (and radial velocity) data.

## Python code

The [`proper-motion-map.py`](proper-motion-map.py) python code can be invoked from the command line for creating the all
sky proper motion maps.
