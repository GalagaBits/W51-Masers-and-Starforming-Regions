'''Reference 
Fitting with spectral-cube and astropy

Eric Koch
'''

import warnings
import numpy as np
import matplotlib.pyplot as plt
from spectral_cube import SpectralCube, Projection
from astropy.io import fits
from astropy.modeling import models, fitting
from astropy.convolution import convolve_models
from astropy import units as u

cube = SpectralCube.read("sample18.fits")

print(np.std(cube[:]))

cube.beam_threshold=0.5
moment0 = cube.moment0(axis=0)

x, y = 3227, 5226
size = 1000

moment0_cutout = moment0[y-size:y+size, x-size:x+size]
moment0_cutout.quicklook(use_aplpy=True)

yy, xx = moment0_cutout.spatial_coordinate_map

p_init_gauss2D = models.Gaussian2D(x_mean=xx[size, size], y_mean=yy[size, size], amplitude=10000, x_stddev=0.02 * u.arcsec, y_stddev=0.02 * u.arcsec)

fit_p = fitting.LevMarLSQFitter()

moment0_cutout_quant = moment0_cutout.value
moment0_cutout_quant[np.isnan(moment0_cutout_quant)] = 0.0

with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    p_gauss2D = fit_p(p_init_gauss2D, xx, yy, moment0_cutout_quant)
    
print(p_gauss2D)

plt.figure(figsize=(18, 6))

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
