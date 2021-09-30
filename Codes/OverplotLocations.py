import sys
import numpy as np
sys.path.append('/orange/adamginsburg/ALMA_IMF/reduction/analysis/')
from quicklook_multicolorbar import quicklook
from spectral_cube import SpectralCube
import pylab as pl
from astropy import units as u

fn = '/orange/adamginsburg/w51/2017.1.00293.S/may2021_imaging/w51n.spw0thru19.14500.robust0.thr0.075mJy.mfs.I.startmod.selfcal7.image.tt0'
img = SpectralCube.read(fn, format='casa_image')
cen = img.shape[1]/2, img.shape[2]/2
fig = pl.figure(figsize=(10,10))
fig, ax = quicklook(img[:,cen[0]-400:cen[0]+400, cen[1]-400:cen[1]+600], fig=fig, inner_stretch='linear')

ax.plot(290.9158956*u.deg, 14.51800069*u.deg, marker='+',  markersize=25, color='pink', transform=ax.get_transform('world'),label='(5,3)')
ax.plot(290.9159260*u.deg, 14.51800888*u.deg, marker='+',  markersize=25, color='pink', transform=ax.get_transform('world'))


ax.plot(290.9159003*u.deg, 14.51801042*u.deg, marker='+',  markersize=25, color='orange', transform=ax.get_transform('world'),label='(6,2)')


ax.plot(290.9158939*u.deg, 14.51800111*u.deg, marker='+',  markersize=25, color='white', transform=ax.get_transform('world'), label='(6,3)')
ax.plot(290.9159262*u.deg, 14.51800959*u.deg, marker='+',  markersize=25, color='white', transform=ax.get_transform('world'))

ax.plot(290.9158986*u.deg, 14.51800816*u.deg, marker='+',  markersize=25, color='yellow', transform=ax.get_transform('world'),label='(7,4)')

ax.plot(290.9158982*u.deg, 14.51801105*u.deg, marker='+',  markersize=25, color='red', transform=ax.get_transform('world'),label='(7,5)')

ax.plot(290.9170759*u.deg, 14.518231700*u.deg, marker='+',  markersize=25, color='blue', transform=ax.get_transform('world'),label='(7,6)')
ax.plot(290.9170768*u.deg, 14.51823206*u.deg, marker='+',  markersize=25, color='blue', transform=ax.get_transform('world'))
ax.plot(290.9170765*u.deg, 14.518231240*u.deg, marker='+',  markersize=25, color='blue', transform=ax.get_transform('world'))

ax.plot(290.9170797*u.deg, 14.51823187*u.deg, marker='+',  markersize=25, color='green', transform=ax.get_transform('world'),label='(7,7)')
ax.plot(290.9170811*u.deg, 14.51823057*u.deg, marker='+',  markersize=25, color='green', transform=ax.get_transform('world'))
ax.plot(290.9170758*u.deg, 14.51823158*u.deg, marker='+',  markersize=25, color='green', transform=ax.get_transform('world'))
ax.plot(290.9170760*u.deg, 14.5182316*u.deg, marker='+',  markersize=25, color='green', transform=ax.get_transform('world'))


ax.legend()
pl.show()