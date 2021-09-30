import sys
import warnings
sys.path.append('/orange/adamginsburg/ALMA_IMF/reduction/analysis/')

import numpy as np
import pylab as pl
from astropy import units as u
from astropy.wcs import WCS
from astropy import visualization

from spectral_cube import SpectralCube

from quicklook_multicolorbar import quicklook

fn = '/orange/adamginsburg/w51/2017.1.00293.S/may2021_imaging/w51n.spw0thru19.14500.robust0.thr0.075mJy.mfs.I.startmod.selfcal7.image.tt0'
fh = SpectralCube.read(fn, format='casa_image')

with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    ww = WCS(fh[0].header)


def w51north_plot():
    fig = pl.figure(figsize=(9, 9))
    ax = pl.subplot(projection=ww)
    im = ax.imshow(fh[0].data, cmap='gray',
                   norm=visualization.simple_norm(fh[0].data, stretch='log', max_percent=100.00))
    cb = pl.colorbar(mappable=im)
    cb.set_label(f"$S_\\nu$ [{fh[0].header['BUNIT']}]", fontsize=16)
    # cb.set_ticks([np.nanmin(fh[0].data), -0.005, 0.00 ,0.005, 0.01, 0.015, 0.020, 0.025])
    ax.axis([6750, 7750, 6750, 7750])

    radesys = ww.wcs.radesys

    _ = ax.set_xlabel(f"Right Ascension {ww.wcs.radesys}")
    _ = ax.set_ylabel(f"Declination {ww.wcs.radesys}")

    tick_fontsize = 14
    fontsize = 16
    ra = ax.coords['ra']
    ra.set_major_formatter('hh:mm:ss.s')
    dec = ax.coords['dec']
    ra.set_axislabel(f"RA ({radesys})", fontsize=fontsize)
    dec.set_axislabel(f"Dec ({radesys})", fontsize=fontsize, minpad=0.0)
    ra.ticklabels.set_fontsize(tick_fontsize)
    ra.set_ticklabel(exclude_overlapping=True)
    dec.ticklabels.set_fontsize(tick_fontsize)
    dec.set_ticklabel(exclude_overlapping=True)

    return ax

#unknown region
cen = fh.shape[1]/2, fh.shape[2]/2

ax = w51north_plot()
ax.axis([cen[1]-215,cen[1]-190,cen[0]+80,cen[0]+105])

ax.imshow(fh[0].data, cmap='gray', norm=visualization.simple_norm(fh[0].data, stretch='log', max_percent=100.00))


#[cen[1]-215,cen[1]-190,cen[0]+80,cen[0]+105]

