# Folder for input data from Gaia EDR3

This folder is expected to contain the input data from Gaia EDR3 needed for the production of the animation. They are
not stored in this Github repository. 

The Gaia EDR3 catalogue fields that should be present in the input data are:
```
ra
dec
parallax
pmra
pmdec
dr2_radial_velocity
phot_g_mean_mag
```

The input data for the Gaia EDR3 release animation was retrieved from the Gaia archive with the following query:
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
