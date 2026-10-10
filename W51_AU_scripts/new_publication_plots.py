# Converted from notebooks/new_publication_plots.ipynb (exact code; one '# %%' block per notebook cell)
#
# Paper item this script produces: Figure 1 (Plots/w51north_overplot_apr25.png) and its scale bars
# (250 au, 250 au, 1400 au, 0.02 pc).
#
# Suggestions are commented out. Uncomment a whole tagged block to apply it:
#   [A2]    size the scale bars for the adopted distance, 5.41 kpc (Sato et al. 2010).
#           Uncomment the helper in cell 2 first; the bar lines in cell 11 use it.
#   [STYLE] write the unit as lowercase "au".
#   [OPTION] also save Figure 1 under the file name main.tex uses.
#
# The figure is saved to a file, never opened.

# Save figures to files instead of opening plot windows (must run before pyplot is imported):
import matplotlib
matplotlib.use('Agg')

# %% cell 0
import numpy as np
import pylab as pl
from astropy import units as u
from astropy.wcs import WCS
from astropy import visualization
from spectral_cube import SpectralCube
import warnings
import pandas as pd
import matplotlib.pyplot as plt

# %% cell 1
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset


# %% cell 2
fn = '/Volumes/GalagaResearch/Developer/masers/continuum/w51n.spw0thru19.14500.robust0.thr0.075mJy.mfs.I.startmod.selfcal7.image.tt0.pbcor'
fh = SpectralCube.read(fn, format='casa_image')
cen = fh.shape[1]/2, fh.shape[2]/2
fh = fh.to(u.K)
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    ww = WCS(fh[0].header)
# [A2] Scale-bar helper used by the [A2] lines in cell 11. Uncomment it first: remove the "#" (and the
# space after it) from every line from `from astropy...` down to `print(...)`, including the function body.
# The old bar widths were typed in by hand in pixels (width=7, 39.2, 140); this computes them
# from the image's pixel scale and the adopted distance: theta = length / D, then theta / pixel size.
# 
from astropy.wcs.utils import proj_plane_pixel_scales
pix = (proj_plane_pixel_scales(ww.celestial)[0] * u.deg).to(u.arcsec)   # arcsec per pixel

def bar_px(length, D=5.41*u.kpc):
    """Return the length of a scale bar in image pixels.

    length : physical size of the bar, with units (e.g. 250*u.au or 0.02*u.pc)
    D      : distance to the source; defaults to 5.41 kpc (Sato et al. 2010)
    """
    theta = (length / D).to(u.arcsec, u.dimensionless_angles())   # angle = size / distance, in arcsec
    n_pixels = (theta / pix).decompose()                           # arcsec / (arcsec per pixel) = pixels
    return n_pixels.value                                          # plain number, which matplotlib needs

print('pixel scale:', pix, '| bar widths (px): 250 au =', bar_px(250*u.au), ', 1400 au =', bar_px(1400*u.au), ', 0.02 pc =', bar_px(0.02*u.pc))

# %% cell 3
from astropy.io import fits
blue_outflow = fits.open('/Volumes/GalagaResearch/Developer/masers/LacyJet_SO65m54_blue30_50.fits', use_dask=True)
red_outflow = fits.open('/Volumes/GalagaResearch/Developer/masers/LacyJet_SO65m54_red70_95.fits', use_dask=True)
from astropy import visualization, wcs

# %% cell 4
#im = ax.imshow(fh[0].data, cmap='Greys', norm=visualization.simple_norm(fh[0].data, stretch='linear', max_percent=99.99))
img = fh[0].data

# %% cell 5
def w51north_plot():

    fig = pl.figure(figsize=(16, 9))
    ax = pl.subplot(projection=ww)
    im = ax.imshow(fh[0].data, cmap='Greys', norm=visualization.simple_norm(fh[0].data, stretch='linear', max_percent=99.99))
    
    #divider = make_axes_locatable(ax)
    #colorbar_axes = divider.append_axes("right",size="10%",pad=0.1)
    
#     cb = pl.colorbar(mappable=im,fraction=0.04, pad=0.04)#cax=colorbar_axes)
#     cb.set_label(f"$T_B$ [{fh[0].header['BUNIT']}]", fontsize = 16)
    

#     #cb.set_ticks([np.nanmin(fh[0].data), -0.005, 0.00 ,0.005, 0.01, 0.015, 0.020, 0.025])
#     cb.ax.tick_params(labelsize=14)
    
    
    # CHANGED: colorbar attached to the image with make_axes_locatable, so it is exactly as tall as the image.
    # (was: cb = pl.colorbar(mappable=im, ax=ax, fraction=0.025, pad=0.01))
    # axes_class=maxes.Axes is needed because `ax` is a WCSAxes; the colorbar must be an ordinary Axes.
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    import matplotlib.axes as maxes
    cax = make_axes_locatable(ax).append_axes("right", size="3%", pad=0.1, axes_class=maxes.Axes)
    cb = pl.colorbar(mappable=im, cax=cax)
    
    cb.set_label(f"$T_B$ [{fh[0].header['BUNIT']}]", fontsize = 16)


    #cb.set_ticks([np.nanmin(fh[0].data), -0.005, 0.00 ,0.005, 0.01, 0.015, 0.020, 0.025])
    cb.ax.tick_params(labelsize=14)
    
    
    
    ax.axis([6450,8050,6850,7750])

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
    ra.set_ticklabel(exclude_overlapping=True)
    dec.set_ticklabel(exclude_overlapping=True)
    ra.tick_params(labelsize=tick_fontsize)
    dec.tick_params(labelsize=tick_fontsize)
    return ax

# %% cell 6
def w51north_settings():

    ax.axis([6750,7750,6750,7750])

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
    ra.set_ticklabel(exclude_overlapping=True)
    dec.set_ticklabel(exclude_overlapping=True)

    ra.tick_params(labelsize=tick_fontsize)
    dec.tick_params(labelsize=tick_fontsize)
    
    return ax

# %% cell 7
def maser_points_ax():
    transform = ax.get_transform('world')
    
    maser_points = pd.read_csv('/Volumes/GalagaResearch/Developer/masers/W51/W51-Masers-and-Starforming-Regions/notebooks/w51_irs2_ammonia_maser_positions_CRVAL_CRPIX_change.csv')
    
    for i in range(len(maser_points)):
        ax.errorbar(maser_points['x_mean'][i], maser_points['y_mean'][i], xerr=maser_points['x_mean_error'][i], yerr=maser_points['y_mean_error'][i], color = maser_points['color'][i], marker = '+', transform=transform,)

    ax.legend(title='NH3 (J,K)', fontsize=16)
    return


# %% cell 8
def maser_points_axins():
    transform = axins.get_transform('world')
    
    # Read the CSV
    maser_points = pd.read_csv('/Volumes/GalagaResearch/Developer/masers/W51/W51-Masers-and-Starforming-Regions/notebooks/w51_irs2_ammonia_maser_positions_CRVAL_CRPIX_change.csv')
    
    # Filter rows where 'region' == 'CA3'
    maser_points_CA3 = maser_points.loc[maser_points['source'] == 'CA3']
    
    # Iterate through filtered rows using iterrows()
    for _, row in maser_points_CA3.iterrows():
        axins.errorbar(
            row['x_mean'],
            row['y_mean'],
            xerr=row['x_mean_error'],
            yerr=row['y_mean_error'],
            color=row['color'],  # Use color from the 'color' column
            marker='+',
            transform=transform,
            label=row['line'], 
        )
    return

# %% cell 9
def maser_points_axins2():
    transform = axins2.get_transform('world')
    
    # Read the CSV
    maser_points = pd.read_csv('/Volumes/GalagaResearch/Developer/masers/W51/W51-Masers-and-Starforming-Regions/notebooks/w51_irs2_ammonia_maser_positions_CRVAL_CRPIX_change.csv')
    
    # Filter rows where 'region' == 'CA3'
    maser_points_d2 = maser_points.loc[maser_points['source'] == 'd2']
    
    # Iterate through filtered rows using iterrows()
    for _, row in maser_points_d2.iterrows():
        axins2.errorbar(
            row['x_mean'],
            row['y_mean'],
            xerr=row['x_mean_error'],
            yerr=row['y_mean_error'],
            color=row['color'],  # Use color from the 'color' column
            marker='+',
            transform=transform,
            label = str(row['line'])
        )
    return

# %% cell 10
def maser_points_axins3():
    transform = axins3.get_transform('world')
    
    # Read the CSV
    maser_points = pd.read_csv('/Volumes/GalagaResearch/Developer/masers/W51/W51-Masers-and-Starforming-Regions/notebooks/w51_irs2_ammonia_maser_positions_CRVAL_CRPIX_change.csv')
    
    # Filter rows where 'region' == 'CA3'
    maser_points_W51North = maser_points.loc[maser_points['source'] == 'W51North']
    
    # Iterate through filtered rows using iterrows()
    for _, row in maser_points_W51North.iterrows():
        axins3.errorbar(
            row['x_mean'],
            row['y_mean'],
            xerr=row['x_mean_error'],
            yerr=row['y_mean_error'],
            color=row['color'],  # Use color from the 'color' column
            marker='+',
            transform=transform,
            label = row['line']
        )
    return

# %% cell 11
#### from mpl_toolkits.axes_grid1 import make_axes_locatable

fig = pl.figure(figsize=(20,12.5))
ax = pl.subplot(projection=ww)
ax = w51north_plot()


maser_points_ax()

axins = inset_axes(ax,
                   loc=2, width=2.4, height=2.4,   # CHANGED: was width=3, height=3 (smaller inset)
                   axes_class=visualization.wcsaxes.core.WCSAxes,
                   axes_kwargs=dict(wcs=ww))
axins.axis([cen[1]-215,cen[1]-190,cen[0]+80,cen[0]+105])

mark_inset(parent_axes=ax, inset_axes=axins,
                   fc="none", ec="b", lw=0.8, alpha=0.6, loc1=3, loc2=1)   # CHANGED: added lw=0.8, alpha=0.6 (thinner, lighter lines)



maser_points_axins()

axins.imshow(img, cmap='gray_r', norm=visualization.simple_norm(fh[0].data, stretch='linear', max_percent=99.999))



#
axins2 = inset_axes(ax,
                   loc=4, width=2.4, height=2.4,   # CHANGED: was width=3, height=3 (smaller inset)
                   axes_class=visualization.wcsaxes.core.WCSAxes,
                   axes_kwargs=dict(wcs=ww))
# axins2.axis([cen[0]-50,cen[0]-0,cen[1]+350,cen[1]+400])
axins2.axis([cen[0]+362,cen[0]+392,cen[1]-40,cen[1]-10])
mark_inset(parent_axes=ax, inset_axes=axins2,
                   fc="none", ec="b", lw=0.8, alpha=0.6, loc1=1, loc2=3)   # CHANGED: added lw=0.8, alpha=0.6 (thinner, lighter lines)

maser_points_axins2()
axins2.imshow(img, cmap='gray_r',norm=visualization.simple_norm(fh[0].data, stretch='linear', max_percent=99.998))






#W51North
axins3 = inset_axes(ax,
                   loc=3, width=2.4, height=2.4,   # CHANGED: was width=3, height=3 (smaller inset)
                   axes_class=visualization.wcsaxes.core.WCSAxes,
                   axes_kwargs=dict(wcs=ww))
# axins2.axis([cen[0]-50,cen[0]-0,cen[1]+350,cen[1]+400])
axins3.axis(axins3.axis([cen[0]-180,cen[0]-20,cen[1]-20,cen[1]+140]))
mark_inset(parent_axes=ax, inset_axes=axins3,
                   fc="none", ec="b", lw=0.8, alpha=0.6, loc1=4, loc2=2)   # CHANGED: added lw=0.8, alpha=0.6 (thinner, lighter lines)

maser_points_axins3()
axins3.imshow(img, cmap='gray_r',norm=visualization.simple_norm(fh[0].data, stretch='linear', max_percent=99.9998))


# divider = make_axes_locatable(axins)
# cax = divider.append_axes('right', size='5%', pad=0.05)
# pl.colorbar(ax2cb, cax)
# ax.arrow(290.9169891, 14.51815923, 0, 0.2,
#          head_width=0, head_length=0,
#          fc='red', ec='red', width=10,
#          transform=ax.get_transform(ww))
#each pixel is 0.007"

# cell size = 0.028 arcsec
axins.arrow(cen[0]-197, cen[1]+82.5, 0, 0.05,
         head_width=0, head_length=0,
         fc='red', ec='red', width=7, transform=axins.get_transform(ww))
# [A2] uncomment to resize the CA3 250 au bar above for 5.41 kpc (needs the [A2] helper in cell 2):
axins.patches[-1].set_data(width=bar_px(250*u.au))

plt.text(cen[0]-197, cen[1]+84, '250 AU',
         color='red', rotation=0,
         transform=axins.get_transform(ww), fontsize=20, ha='center')
# [STYLE] uncomment to write the unit as lowercase au (IAU/AAS style, as in the paper text):
# plt.gca().texts[-1].set_text('250 au')


ax.arrow(cen[0]+210, cen[1]-350, 0, 0.2,
         head_width=0, head_length=0,
         fc='red', ec='red', width=140,
         transform=ax.get_transform(ww))
# [A2] uncomment to resize the 0.02 pc bar above for 5.41 kpc (needs the [A2] helper in cell 2).
# At 0.007"/pixel, width=140 is 0.98" = 0.02 pc only at ~4.2 kpc, so this bar was off even at 5.1 kpc:
ax.patches[-1].set_data(width=bar_px(0.02*u.pc))

plt.text(cen[0]+210,cen[1]-330, '0.02 pc',
         color='red', rotation=0,
         transform=ax.get_transform(ww), fontsize=20, ha='center')

axins2.arrow(cen[0]+384, cen[1]-39, 0, 0.05,
         head_width=0, head_length=0,
         fc='red', ec='red', width=7, transform=axins2.get_transform(ww))
# [A2] uncomment to resize the d2 250 au bar above for 5.41 kpc (needs the [A2] helper in cell 2):
axins2.patches[-1].set_data(width=bar_px(250*u.au))
plt.text(cen[0]+384, cen[1]-37, '250 AU',
         color='red', rotation=0,
         transform=axins2.get_transform(ww), fontsize=20, ha='center')
# [STYLE] uncomment to write the unit as lowercase au:
# plt.gca().texts[-1].set_text('250 au')


axins3.arrow(cen[0]-70, cen[1]-10, 0, 0.05,
         head_width=0, head_length=0,
         fc='red', ec='red', width=39.2,#714 au
         transform=axins3.get_transform(ww))
# [A2] uncomment to resize the W51 North 1400 au bar above for 5.41 kpc (needs the [A2] helper in cell 2):
axins3.patches[-1].set_data(width=bar_px(1400*u.au))
plt.text(cen[0]-70, cen[1]+0, '1400 AU',
         color='red', rotation=0,
         transform=axins3.get_transform(ww), fontsize=20, ha='center')
# [STYLE] uncomment to write the unit as lowercase au:
# plt.gca().texts[-1].set_text('1400 au')

ax.contour(red_outflow[0].data, levels=[0.21, 0.42, 0.82, 1.30, 2.14], colors=[(1,0.3,0.3,0.7)]*6, transform=ax.get_transform(wcs.WCS(red_outflow[0].header)))
ax.contour(blue_outflow[0].data, levels=[0.45, 1.69, 3.21], colors=[(0.25,0.5,1,0.8)]*6, transform=ax.get_transform(wcs.WCS(blue_outflow[0].header)))

# axins3.contour(red_outflow[0].data, levels=[0.21, 0.42, 0.82, 1.30, 2.14],
#                colors=[(1,0.3,0.3,0.7)]*6,
#                transform=axins3.get_transform(wcs.WCS(red_outflow[0].header)))
# axins3.contour(blue_outflow[0].data, levels=[0.45, 1.69, 3.21],
#                colors=[(0.25,0.5,1,0.8)]*6,
#                transform=axins3.get_transform(wcs.WCS(blue_outflow[0].header)))

ra = axins2.coords['ra']
dec = axins2.coords['dec']

axins.xaxis.set_visible(False)
axins.yaxis.set_visible(False)
ra.set_ticks_visible(False)
dec.set_ticks_visible(False)
ra.set_axislabel('')
dec.set_axislabel('')
ra.set_ticklabel_visible(False)
dec.set_ticklabel_visible(False)


ra = axins.coords['ra']
dec = axins.coords['dec']

axins.xaxis.set_visible(False)
axins.yaxis.set_visible(False)
ra.set_ticks_visible(False)
dec.set_ticks_visible(False)
ra.set_axislabel('')
dec.set_axislabel('')
ra.set_ticklabel_visible(False)
dec.set_ticklabel_visible(False)

ra = axins3.coords['ra']
dec = axins3.coords['dec']

axins3.xaxis.set_visible(False)
axins3.yaxis.set_visible(False)
ra.set_ticks_visible(False)
dec.set_ticks_visible(False)
ra.set_axislabel('')
dec.set_axislabel('')
ra.set_ticklabel_visible(False)
dec.set_ticklabel_visible(False)


from matplotlib.lines import Line2D

# Define the colors, markers, and labels
colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple',
          'tab:cyan', 'tab:olive', 'tab:gray', 'tab:brown', 'tab:pink']
labels = ['(5,3)', '(6,2)', '(6,3)', '(7,4)', '(7,5)',
          '(7,6)', '(7,7)', '(8,5)', '(9,8)', '(11,9)']

# Create custom legend handles
handles = [Line2D([0], [0], color=color, marker='+', linestyle='', markersize=10, label=label)
           for color, label in zip(colors, labels)]

# Create a figure and add the legend
# CHANGED: added ncol=2 so the legend is half as tall and clears the d2 inset
legend = ax.legend(handles=handles, title="NH3 (J,K)", loc='upper right', fontsize=16, title_fontsize=16, ncol=2)

# plt.show()  # removed: figures are saved to files instead of opened

# pl.show()  # removed: figures are saved to files instead of opened
pl.rcParams['figure.dpi'] = 600
pl.rcParams['savefig.dpi'] = 600
# Changed from fig.savefig(...): `fig` is the figure made at the top of cell 11, but w51north_plot()
# draws everything on a NEW figure, so fig.savefig saved an empty plot. ax.figure is the one with the plot.
ax.figure.savefig('Overview_W51IRS2_6_5410pc', dpi=500, bbox_inches='tight', pad_inches=0.02)
# [OPTION] uncomment to also save it under the file name main.tex uses for Figure 1:
# ax.figure.savefig('w51north_overplot_apr25.png', dpi=800)

# %% cell 12