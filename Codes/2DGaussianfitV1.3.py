import warnings
import numpy as np
import matplotlib.pyplot as plt
from spectral_cube import SpectralCube, Projection
from astropy.io import fits
from astropy.modeling import models, fitting
from astropy.convolution import convolve_models
from astropy import units as u

from astropy.visualization import quantity_support

cube = SpectralCube.read("/orange/adamginsburg/w51/vla/19A-254/derod/sample18.fits")

cube.beam_threshold = 0.5
moment0 = cube.moment0(axis=0)
# maser_channel = 579

x, y = 3227, 5226
size = 1000

moment0_cutout = moment0[y - size:y + size, x - size:x + size]
moment0_cutout.quicklook(use_aplpy=True)

yy, xx = moment0_cutout.spatial_coordinate_map

x_mean = xx[size, size]
x_mean2 = x_mean.to(u.arcsec)
y_mean = yy[size, size]
y_mean2 = y_mean.to(u.arcsec)

#amplitude1 = u.pixel_scale(10000)

p_init_gauss2D = models.Gaussian2D(x_mean=x_mean2, y_mean=y_mean2, amplitude=1000,
                                   x_stddev=0.02 * u.arcsec, y_stddev=0.02 * u.arcsec)



fit_p = fitting.LevMarLSQFitter()

moment0_cutout_quant = moment0_cutout
moment0_cutout_quant[np.isnan(moment0_cutout_quant)] = 0.0

with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    p_gauss2D = fit_p(p_init_gauss2D, xx, yy, moment0_cutout_quant)

print(np.std(cube[:]))
print(p_gauss2D)
cord = cube.world[:, y, x]
print(cord)

quantity_support()

plt.figure(figsize=(18, 6))
plt.title("W51North NH3 (7,7)")
plt.subplot(1, 3, 1)
plt.title("Image")
plt.imshow(moment0_cutout.value, origin='lower', cmap='inferno')
plt.colorbar()
plt.subplot(1, 3, 2)
plt.title("Model")
plt.imshow(p_gauss2D(xx, yy).value, origin='lower', cmap='inferno')
plt.colorbar()
plt.subplot(1, 3, 3)
plt.title("Residual")
plt.imshow(moment0_cutout.value - p_gauss2D(xx, yy).value, origin='lower', cmap='inferno')
plt.colorbar()

plt.tight_layout()

plt.show()