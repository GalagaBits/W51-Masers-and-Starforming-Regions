# Converted from notebooks/SED_modeling.ipynb (exact code; one '# %%' block per notebook cell)
#
# Paper values this script produces: the CA3 radius r = 48 au (main.tex lines 263, 279, 288) with its
# density and Nlyc, and the two SED panels combined into Plots/SED_modeling.001.png (Figure 4).
#
# Results print to the terminal; every figure is saved to a PNG file, never opened. The 2D Gaussian
# fits in cell 2 save one check plot each: gaussian2D_fit_<source>_<freq>GHz.png.
#
# Suggestions are commented out. Uncomment a whole tagged block to apply it:
#   [A2]   use the adopted distance, 5.41 kpc (Sato et al. 2010).
#   [NOTE] information only; nothing to uncomment.

# Save figures to files instead of opening plot windows (must run before pyplot is imported):
import matplotlib
matplotlib.use('Agg')

# %% cell 0
import os
import sys
sys.path.append('/Volumes/GalagaResearch/Developer/masers/W51/W51-Masers-and-Starforming-Regions/notebooks/W51CA3_Fitting/hiimodel')

from HII_model import HIIregion
import pandas as pd
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np


import numpy as np
import warnings
import matplotlib.pyplot as plt
from spectral_cube import SpectralCube, Projection
from astropy.modeling import models, fitting
from astropy.convolution import convolve_models
from astropy import units as u
from astropy.io import fits
import radio_beam
from radio_beam import Beam
import pandas as pd
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import pixel_to_skycoord

plt.rcParams['lines.linewidth'] = 3
plt.rcParams['lines.markersize'] = 50

# %% cell 1
#Accessing Cube Data
def fit_2dgaussian(directory, maser_table, freq, amplitude, x, y, x_stddev, y_stddev, 
    theta=np.pi/3, size=13, channel=0, cube=True, format='fits', table=True, source_name='None',
    fixed_xstddev=False, fixed_ystddev=False, fixed_theta=False):

    """
    Fit a 2D Gaussian to a cutout of a spectral cube or image.

    """

    if cube:
        
        cube = SpectralCube.read(directory, format=format)

        if channel is None:
            pass
        else:
            cube = cube[channel,:,:]

        cube_cutout = cube[y - size:y + size, x - size:x + size]
        yy, xx = cube_cutout.spatial_coordinate_map

        cord = cube.world[y,x]

        x_guess, y_guess = cord[1], cord[0]

    else:
        
        hdul = fits.open(directory)

        hdu = hdul[0]

        data_cube = hdul[0].data

        if data_cube.ndim == 4:
            data = data_cube[0, 0, :, :] * u.Jy/u.beam

        else:
            data = hdu.data * u.Jy/u.beam

        header = hdu.header

        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            ww = WCS(header)

        data_cutout = data[y - size:y + size, x - size:x + size]
        ypixel_grid, xpixel_grid = np.mgrid[y - size:y + size, x - size:x + size]
        
        # xx, yy = ww.pixel_to_world_values(xpixel_grid, ypixel_grid)

        # xx = (xx * u.deg)
        # yy = (yy * u.deg)

        # x_guess, y_guess = xx[size, size], ypixel_grid[size, size]

    p_init_gauss2D = models.Gaussian2D(x_mean=x, y_mean=y, amplitude=amplitude * (u.Jy/u.beam),
                                    x_stddev=(x_stddev), y_stddev=(y_stddev), theta=(theta)*u.rad)

    if fixed_xstddev:
        p_init_gauss2D.x_stddev.fixed = True

    if fixed_ystddev:
        p_init_gauss2D.y_stddev.fixed = True

    if fixed_theta:
        p_init_gauss2D.theta.fixed = True

    fit_p = fitting.LevMarLSQFitter()
    data_cutout[np.isnan(data_cutout)] = 0.

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        p_gauss2D = fit_p(p_init_gauss2D, xpixel_grid, ypixel_grid, data_cutout)

  
    cov = np.diag(fit_p.fit_info['param_cov'])
    errors = np.sqrt(cov)

    amplitude_error=np.format_float_scientific(errors[0],precision=9)
    x_mean_error=np.format_float_scientific(errors[1],precision=9)
    y_mean_error=np.format_float_scientific(errors[2],precision=9)
    x_stddev_error=np.format_float_scientific(errors[3],precision=9)
    y_stddev_error=np.format_float_scientific(errors[4],precision=9)
    theta_error=np.format_float_scientific(errors[5],precision=9)

    amplitude = np.format_float_scientific(p_gauss2D.amplitude[0],precision=9)
    x_mean = np.format_float_scientific(p_gauss2D.x_mean[0],precision=9)
    y_mean = np.format_float_scientific(p_gauss2D.y_mean[0],precision=9)
    x_stddev = np.format_float_scientific(p_gauss2D.x_stddev[0],precision=9)
    y_stddev = np.format_float_scientific(p_gauss2D.y_stddev[0],precision=9)
    theta = np.format_float_scientific(p_gauss2D.theta[0],precision=9)

                
    plt.figure(figsize=(20, 6))
    plt.suptitle(f"Source-{source_name} freq-{freq}", fontsize=18)

    plt.subplot(1, 3, 1)
    plt.title("Image")
    plt.imshow(data_cutout.value, origin='lower', cmap='inferno')
    plt.colorbar()

    plt.subplot(1, 3, 2)
    plt.title("Model")
    plt.imshow(p_gauss2D(xpixel_grid, ypixel_grid).value, origin='lower', cmap='inferno')
    plt.colorbar()

    plt.subplot(1, 3, 3)
    plt.title("Residual")
    plt.imshow(data_cutout.value - p_gauss2D(xpixel_grid, ypixel_grid).value, origin='lower', cmap='inferno')
    plt.colorbar(label='S (Jy)')

    plt.savefig(f"gaussian2D_fit_{source_name}_{freq}GHz.png")  # was plt.show(): save instead of opening a window
    plt.close()

    print("Parameter" + ' ' + 'Error' +"\n")
    print("Amplitude" +"\n")
    print(amplitude + ' ' + amplitude_error +"\n")

    print("x_mean"  +"\n")
    print(x_mean + ' ' + x_mean_error +"\n")

    print("y_mean" +"\n")
    print(y_mean + ' ' + y_mean_error +"\n")

    print("x_stddev" +"\n")
    print(x_stddev + ' ' + x_stddev_error +"\n")

    print("y_stddev" +"\n")
    print(y_stddev + ' ' + y_stddev_error +"\n")

    print("theta" +"\n")
    print(theta + ' ' + theta_error +"\n")

    # print('\n\nSkyCoord Conversion \n----------------- \n')
    # print("Parameter" + ' ' + 'Error' +"\n")

    skycoord = pixel_to_skycoord(float(x_mean), float(y_mean), ww)
    skycoord_err = SkyCoord(ra=x_mean_error, dec=y_mean_error, unit=(u.deg, u.deg), frame=ww.wcs.radesys.lower())

    # error propagation? # Not sure if this is correct
    skycoord_error = pixel_to_skycoord(float(x_mean_error), float(y_mean_error), ww)


    ra_hms_str = skycoord.ra.to_string(unit=u.hour, sep=':', precision=4)
    dec_dms_str = skycoord.dec.to_string(unit=u.deg, sep=':', precision=4)


    print(skycoord_error)

    if table:
        new_row = pd.DataFrame({'source': [source_name] , 
            'freq': [freq], 'x_pix': [x], 'y_pix': [y], 'amplitude': [amplitude], 
            'amplitude_error': [amplitude_error], 'x_mean': [x_mean], 'x_mean_error': [x_mean_error], 
            'y_mean': [y_mean], 'y_mean_error': [y_mean_error], 
            'x_mean_sky': [skycoord.ra.deg], 'x_mean_sky_error': [skycoord_error.ra.deg],
            'y_mean_sky': [skycoord.dec.deg], 'y_mean_sky_error': [skycoord_error.dec.deg],
            'ra_hms': [ra_hms_str], 'dec_dms': [dec_dms_str],
            'SkyCoord' : [skycoord],
            'frame': [skycoord.frame.name], 
            'header': [ww]})
        maser_table = pd.concat([maser_table, new_row], ignore_index=True)

        return maser_table

# %% cell 2
ca3_positions_1mm = pd.DataFrame(columns=['source','freq', 'x_pix', 'y_pix', 'amplitude', 'amplitude_error', 
'x_mean', 'x_mean_error', 'y_mean', 'y_mean_error', 'x_mean_sky', 'x_mean_sky_error', 
'y_mean_sky', 'y_mean_sky_error', 
'ra_hms', 'dec_dms',
'SkyCoord', 'frame', 'header'])

directory3 = '/Volumes/GalagaResearch/Developer/masers/continuum/W51n_cont_briggsSC_tclean.image.fits'

ca3_positions_1mm = fit_2dgaussian(table=True, cube=False, directory=directory3, 
                                          maser_table=ca3_positions_1mm,
                                          amplitude=1, x=4654, y=3311, x_stddev=2, y_stddev=3, theta=np.pi/3, 
                                          size=6, freq='226.4187', source_name='d2')

ca3_positions_1mm = fit_2dgaussian(table=True, cube=False, directory=directory3, 
                                          maser_table=ca3_positions_1mm,
                                          amplitude=1, x=3295, y=3580, x_stddev=3, y_stddev=3, theta=np.pi/4, 
                                          size=8, freq='226.4187', source_name='CA3')

ca3_positions_3mm = pd.DataFrame(columns=['source','freq', 'x_pix', 'y_pix', 'amplitude', 'amplitude_error', 
'x_mean', 'x_mean_error', 'y_mean', 'y_mean_error', 'x_mean_sky', 'x_mean_sky_error', 
'y_mean_sky', 'y_mean_sky_error', 
'ra_hms', 'dec_dms',
'SkyCoord', 'frame', 'header'])

directory3 = '/Volumes/GalagaResearch/Developer/masers/continuum/w51n.spw0thru19.14500.robust0.thr0.075mJy.mfs.I.startmod.selfcal7.image.tt0.pbcor.fits'
# directory3 = '/Volumes/GalagaResearch/Developer/masers/continuum/modified_VLA/W51_North_KbandAarray_cont_spws_continuum_cal_clean_2terms_robust0_selfcal3.image.tt0.pbcor_ICRS_CRVAL_CA3_crval_change.fits'

# d2
ca3_positions_3mm = fit_2dgaussian(table=True, cube=False, directory=directory3, 
                                          maser_table=ca3_positions_3mm,
                                          amplitude=5e-2, x=7628, y=7226, x_stddev=12, y_stddev=12, theta=11*np.pi/6, 
                                          size=10, freq='92.9823', source_name='d2')

# CA3
ca3_positions_3mm = fit_2dgaussian(table=True, cube=False, directory=directory3, 
                                          maser_table=ca3_positions_3mm,
                                          amplitude=1, x=7047, y=7342, x_stddev=9, y_stddev=9, theta=np.pi/4, 
                                          size=8, freq='92.9823', source_name='CA3')

ca3_positions_Kband = pd.DataFrame(columns=['source','freq', 'x_pix', 'y_pix', 'amplitude', 'amplitude_error', 
'x_mean', 'x_mean_error', 'y_mean', 'y_mean_error', 'x_mean_sky', 'x_mean_sky_error', 
'y_mean_sky', 'y_mean_sky_error', 
'ra_hms', 'dec_dms',
'SkyCoord', 'frame', 'header'])

directory3 = '/Volumes/GalagaResearch/Developer/masers/continuum/modified_VLA/W51_North_KbandAarray_cont_spws_continuum_cal_clean_2terms_robust0_selfcal3.image.tt0.pbcor_ICRS_CRVAL_CA3_crval_change.fits'
# directory3 = '/Volumes/GalagaResearch/Developer/masers/continuum/modified_VLA/W51_North_KbandAarray_cont_spws_continuum_cal_clean_2terms_robust0_selfcal3.image.tt0.pbcor_ICRS_CRVAL_CA3_crval_change.fits'

# d2
ca3_positions_Kband = fit_2dgaussian(table=True, cube=False, directory=directory3, 
                                          maser_table=ca3_positions_Kband,
                                          amplitude=1, x=5219, y=4959, x_stddev=3, y_stddev=3, theta=np.pi/4, 
                                          size=6, freq='22.4601', source_name='d2')

# CA3
ca3_positions_Kband = fit_2dgaussian(table=True, cube=False, directory=directory3, 
                                          maser_table=ca3_positions_Kband,
                                          amplitude=1, x=4948, y=5015, x_stddev=3, y_stddev=3, theta=np.pi/4, 
                                          size=4, freq='22.4601', source_name='CA3')


ca3_positions_Qband = pd.DataFrame(columns=['source','freq', 'x_pix', 'y_pix', 'amplitude', 'amplitude_error', 
'x_mean', 'x_mean_error', 'y_mean', 'y_mean_error', 'x_mean_sky', 'x_mean_sky_error', 
'y_mean_sky', 'y_mean_sky_error', 
'ra_hms', 'dec_dms',
'SkyCoord', 'frame', 'header'])

directory3 = '/Volumes/GalagaResearch/Developer/masers/continuum/W51_North_QbandAarray_cont_spws_continuum_cal_clean_2terms_robust0_selfcal4_final.image.tt0.pbcor_ICRS_CRVAL_CA3_crval_change.fits'

ca3_positions_Qband = fit_2dgaussian(table=True, cube=False, directory=directory3, 
                                          maser_table=ca3_positions_Qband,
                                          amplitude=1, x=8158, y=7381, x_stddev=2, y_stddev=3, theta=np.pi/3, 
                                          size=10, freq='46.99', source_name='d2')

ca3_positions_Qband = fit_2dgaussian(table=True, cube=False, directory=directory3, 
                                          maser_table=ca3_positions_Qband,
                                          amplitude=1, x=7345, y=7544, x_stddev=3, y_stddev=3, theta=np.pi/4, 
                                          size=12, freq='46.99', source_name='CA3')

ca3_positions_Kaband = pd.DataFrame(columns=['source','freq', 'x_pix', 'y_pix', 'amplitude', 'amplitude_error', 
'x_mean', 'x_mean_error', 'y_mean', 'y_mean_error', 'x_mean_sky', 'x_mean_sky_error', 
'y_mean_sky', 'y_mean_sky_error', 
'ra_hms', 'dec_dms',
'SkyCoord', 'frame', 'header'])

directory3 = '/Volumes/GalagaResearch/Developer/masers/continuum/modified_VLA/W51_North_KabandAarray_cont_spws_continuum_cal_clean_2terms_robust0_selfcal2.image.tt0.pbcor_ICRS_CRVAL_CA3_crval_change.fits'

ca3_positions_Kaband = fit_2dgaussian(table=True, cube=False, directory=directory3, 
                                          maser_table=ca3_positions_Kaband,
                                          amplitude=1, x=5219, y=4960, x_stddev=2, y_stddev=3, theta=np.pi/3, 
                                          size=5, freq='33.0368', source_name='d2')

ca3_positions_Kaband = fit_2dgaussian(table=True, cube=False, directory=directory3, 
                                          maser_table=ca3_positions_Kaband,
                                          amplitude=1, x=4949, y=5015, x_stddev=3, y_stddev=3, theta=np.pi/4, 
                                          size=5, freq='33.0368', source_name='CA3')

# %% cell 3
tables = [
    ca3_positions_1mm,
    ca3_positions_3mm,
    ca3_positions_Kband,
    ca3_positions_Qband,
    ca3_positions_Kaband,
]

ca3_all_d2 = pd.concat(
    [t[t["source"] == "d2"] for t in tables],   # list-comprehension of filtered dfs
    ignore_index=True                           # re-index from 0,1,2,…
)
ca3_all_d2['freq'] = ca3_all_d2['freq'].astype(float)
ca3_all_d2 = ca3_all_d2.sort_values(by='freq', ascending=False)



ca3_all_CA3 = pd.concat(
    [t[t["source"] == "CA3"] for t in tables],   # list-comprehension of filtered dfs
    ignore_index=True                           # re-index from 0,1,2,…
)
ca3_all_CA3['freq'] = ca3_all_CA3['freq'].astype(float)
ca3_all_CA3 = ca3_all_CA3.sort_values(by='freq', ascending=False)

# %% cell 4
print(ca3_all_CA3)

# %% cell 5
print(ca3_all_d2['amplitude'].astype(float) * 1000)

# %% cell 6
print(ca3_all_d2['amplitude_error'].astype(float) * 1000)

# %% cell 7
print(ca3_all_CA3)

# %% cell 8
print(ca3_all_CA3['freq'].astype(float))

# %% cell 9
freq_CA3 = ca3_all_CA3['freq'].astype(float)*u.GHz
amp_CA3 = ca3_all_CA3['amplitude'].astype(float)*1000*u.mJy
amp_error_CA3 = ca3_all_CA3['amplitude_error'].astype(float)*1000*u.mJy

freq_d2 = ca3_all_d2['freq'].astype(float)*u.GHz
amp_d2 = ca3_all_d2['amplitude'].astype(float)*1000*u.mJy
amp_error_d2 = ca3_all_d2['amplitude_error'].astype(float)*1000*u.mJy

# %% cell 10
print(freq_CA3)

# %% [markdown] cell 11
# # CA3

# %% cell 12
data = pd.read_csv('/Volumes/GalagaResearch/Developer/masers/W51/W51-Masers-and-Starforming-Regions/notebooks/W51CA3_Fitting/hiimodel/hiiregion_W51_CA3_pbcor_03-with_Updated_band6.csv')

# %% cell 13
print(data['mu_GHz'])

# %% cell 14
plt.rcParams['lines.linewidth'] = 3
plt.rcParams['lines.markersize'] = 50

fig , ax = plt.subplots(figsize=(10,10), dpi=300)

ca3 = HIIregion(freq_CA3, 
                amp_CA3, 
                fluxerr=amp_error_CA3,
                dist_kpc=5.1*u.kpc, EMguess=1e10)
# [A2] uncomment to refit CA3 at the adopted distance (replaces the fit above, so the plot shows the new size):
ca3 = HIIregion(freq_CA3, amp_CA3, fluxerr=amp_error_CA3, dist_kpc=5.41*u.kpc, EMguess=1e10)

ca3.loglogplot(numin=5*u.GHz, numax=250*u.GHz, annotation_xpos=0.5)
ax.text(0.1,.9,"CA3",fontsize = 20, weight='bold',transform=ax.transAxes)



plt.xlim(20,250)
plt.ylim(0.5, 10)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xlabel(xlabel='Frequency (GHz)',fontsize=20)
plt.ylabel(ylabel='Flux Density (mJy)',fontsize=20)
fig.savefig(fname='CA3_HII_region_fit.png', transparent=False)
print('CA3:', ca3.dist_kpc, ca3.physprops())
# [NOTE] HII_model.physprops() takes the size from the LOWEST-frequency flux (22.46 GHz here),
# [NOTE] not the 3 mm (92.98 GHz) flux that main.tex line 263 says was used.

# %% [markdown] cell 15
# # d2

# %% cell 16
data = pd.read_csv('/Volumes/GalagaResearch/Developer/masers/W51/W51-Masers-and-Starforming-Regions/notebooks/W51CA3_Fitting/hiimodel/sed_d2_pbcor.csv')

# %% cell 17
fig , ax = plt.subplots(figsize=(10,10), dpi=300)

# ca3 = HIIregion(data['mu_GHz']*u.GHz, data['amp_Jy']*1000*u.mJy, fluxerr=data['amp_Jy_error']*1000*u.mJy,
#                 dist_kpc=5.1*u.kpc, EMguess=1e10)

d2 = HIIregion(freq_d2, 
                amp_d2, 
                fluxerr=amp_error_d2,
                dist_kpc=5.1*u.kpc, EMguess=1e10)
# [A2] uncomment to refit d2 at the adopted distance (replaces the fit above, so the plot shows the new size):
d2 = HIIregion(freq_d2, amp_d2, fluxerr=amp_error_d2, dist_kpc=5.41*u.kpc, EMguess=1e10)

d2.loglogplot(numin=5*u.GHz, numax=250*u.GHz, annotation_xpos=0.5)
ax.text(0.1,.9,"d2",fontsize = 20, weight='bold',transform=ax.transAxes)

plt.xlim(5,250)

plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xlabel(xlabel='Frequency (GHz)',fontsize=20)
plt.ylabel(ylabel='Flux Density (mJy)',fontsize=20)
fig.savefig(fname='d2_HII_region_fit.png', transparent=False)
print('d2:', d2.dist_kpc, d2.physprops())

# %% cell 18

