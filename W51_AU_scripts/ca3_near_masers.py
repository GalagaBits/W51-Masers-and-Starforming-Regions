# Converted from notebooks/ca3_near_masers.ipynb (exact code; one '# %%' block per notebook cell)
#
# Paper values this script produces: the maser separations (in arcsec and au) on main.tex
# lines 161 (CA3), 166 (W51 North, incl. position angles) and 169 (d2).
#
# Results print to the terminal; figures are saved to files, never opened.
#
# Suggestions are commented out. Uncomment a whole tagged block to apply it:
#   [A2]   use the adopted distance, 5.41 kpc = 5410 pc (Sato et al. 2010). au = arcsec x pc.
#   [FIX]  correct a bug found while converting.

# Save figures to files instead of opening plot windows (must run before pyplot is imported):
import matplotlib
matplotlib.use('Agg')

# %% cell 0
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.coordinates import SkyCoord


masers = pd.read_csv('/Volumes/GalagaResearch/Developer/masers/W51/W51-Masers-and-Starforming-Regions/notebooks/w51_irs2_ammonia_maser_spectral_spatial_fits.csv')
sources = pd.read_csv('/Volumes/GalagaResearch/Developer/masers/W51/W51-Masers-and-Starforming-Regions/notebooks/FK5_ICRS_offsets/offset_sources.csv')

print(masers)

# %% cell 1
print(masers.columns)

# %% cell 2
print(masers[(masers['source'] == 'CA3')])

# %% cell 3
ra_ca3_masers = masers[(masers['source'] == 'CA3')]['x_mean'].values

# %% cell 4
dec_ca3_masers = masers[(masers['source'] == 'CA3')]['y_mean'].values

# %% cell 5
coords_ca3 = SkyCoord(ra_ca3_masers, dec_ca3_masers, unit=(u.deg, u.deg), frame='icrs')

# %% cell 6
ra_d2_masers = masers[(masers['source'] == 'd2')]['x_mean'].values
dec_d2_masers = masers[(masers['source'] == 'd2')]['y_mean'].values
coords_d2 = SkyCoord(ra_d2_masers, dec_d2_masers, unit=(u.deg, u.deg), frame='icrs')

ra_w51north_masers = masers[(masers['source'] == 'W51North')]['x_mean'].values
dec_w51north_masers = masers[(masers['source'] == 'W51North')]['y_mean'].values
coords_w51north = SkyCoord(ra_w51north_masers, dec_w51north_masers, unit=(u.deg, u.deg), frame='icrs')

# %% cell 7
print(masers[(masers['source'] == 'W51North')])

# %% cell 8


# %% cell 9
# from Kband_separation.ipynb. 

ra_ca3_Kband_source = 290.9170697
dec_ca3_Kband_source = 14.51823865

ca3_Kband_source = SkyCoord(ra_ca3_Kband_source, dec_ca3_Kband_source, unit=(u.deg, u.deg), frame='icrs')

# %% cell 10
distances = []
separations = []
for i in zip(ra_ca3_masers, dec_ca3_masers):
    # print(ra_ca3_Kband_source - i[0])
    # print(dec_ca3_Kband_source - i[1])

    distance = np.sqrt((ra_ca3_Kband_source - i[0])**2 + (dec_ca3_Kband_source - i[1])**2)
    print(distance * u.deg.to(u.arcsec))
    distances.append(distance * u.deg.to(u.arcsec))

distance_from_ca3 = np.mean(distances) * u.arcsec

for i in range(len(distances)):
    coords_ca3[i].separation(ca3_Kband_source).to(u.arcsec)
    separation = coords_ca3[i].separation(ca3_Kband_source).to(u.arcsec)
    separations.append(separation.value)

separations = np.array(separations)
print('CA3 max separation (arcsec), au @ 5128 pc:', np.max(separations), np.max(separations) * 5128)
# [A2] CA3 max separation at 5.41 kpc (new value for line 161):
print('CA3 max separation (arcsec), au @ 5410 pc:', np.max(separations), np.max(separations) * 5410)

# %% cell 11
print('CA3 mean separation:', distance_from_ca3)

# %% cell 12
distance = 5128 * u.pc
# [A2] uncomment to use the adopted distance for the line below:
distance = 5410 * u.pc

print('CA3 mean separation (au) at `distance`:', (distance_from_ca3 * distance).value)

# %% cell 13
# from Kband_separation.ipynb. 

ra_d2_Kband_source = 290.91590672
dec_d2_Kband_source = 14.51800645
d2_Kband_source = SkyCoord(ra_d2_Kband_source, dec_d2_Kband_source, unit=(u.deg, u.deg), frame='icrs')

separations = []

for i in range(len(distances)):
    coords_d2[i].separation(d2_Kband_source).to(u.arcsec)
    separation = coords_d2[i].separation(d2_Kband_source).to(u.arcsec)
    separations.append(separation.value)

separations = np.array(separations)
print('d2 max separation (arcsec), au @ 5128 pc:', np.max(separations), np.max(separations) * 5128)
# [FIX] The loop above runs range(len(distances)) = 7 (the number of CA3 masers), so it skips
# the last 2 of the 9 d2 masers, (8,5) and (9,8). Uncomment to use all 9. The maximum does not
# change (0.0630", the second (6,3) spot), but this makes the calculation correct:
# separations = coords_d2.separation(d2_Kband_source).to(u.arcsec).value
# [A2] d2 max separation at 5.41 kpc (new value for line 169):
print('d2 max separation (arcsec), au @ 5410 pc:', np.max(separations), np.max(separations) * 5410)

# %% cell 14
print(coords_w51north)

# %% cell 15
ra_w51north_Kband_source = 290.91687133
dec_w51north_Kband_source = 14.51818835
w51north_Kband_source = SkyCoord(ra_w51north_Kband_source, dec_w51north_Kband_source, unit=(u.deg, u.deg), frame='icrs')

separations = []
pa = []

for i in range(0,2):
    coords_w51north[i].separation(w51north_Kband_source).to(u.arcsec)
    separation = coords_w51north[i].separation(w51north_Kband_source).to(u.arcsec)
    separations.append(separation.value)

    position_angle = coords_w51north[i].position_angle(w51north_Kband_source).to(u.deg)
    # [FIX] The line above gives the direction of the CENTER as seen from the maser, which is off by
    # 180 deg. Uncomment to get the P.A. of the maser measured from the center (east of north):
    # position_angle = w51north_Kband_source.position_angle(coords_w51north[i]).to(u.deg)
    pa.append(position_angle.value)



separations = np.array(separations)
print('W51 North farther maser: separation (arcsec), au @ 5128 pc:', np.max(separations), np.max(separations) * 5128)
# [A2] the same at 5.41 kpc (new value for line 166):
print('W51 North farther maser: separation (arcsec), au @ 5410 pc:', np.max(separations), np.max(separations) * 5410)

# %% cell 16


# %% cell 17
print('W51 North nearer maser: separation (arcsec), au @ 5128 pc:', separations[0], separations[0] * 5128)
# [A2] the same at 5.41 kpc (new value for line 166):
print('W51 North nearer maser: separation (arcsec), au @ 5410 pc:', separations[0], separations[0] * 5410)

# %% cell 18
print('W51 North mean separation (arcsec):', np.mean(separations))

# %% cell 19
print('W51 North P.A. (deg):', pa)

# %% cell 20


# %% cell 21

