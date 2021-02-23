import warnings
import numpy as np
import matplotlib.pyplot as plt
from spectral_cube import SpectralCube, Projection
from astropy.io import fits
from astropy.modeling import models, fitting
from astropy.convolution import convolve_models
from astropy import units as u

cube = SpectralCube.read("/Users/dealderod/Desktop/W51North_spw_6.fits")

print(cube)

moment0 = cube.moment0(axis=0)

moment0.quicklook(use_aplpy=True)

y, x = 505, 485
size = 10

moment0_cutout = moment0[y-size:y+size, x-size:x+size]
moment0_cutout.quicklook(use_aplpy=True)

# Define the spatial grid for the fit centered at y, x = 145, 342
yy, xx = moment0_cutout.spatial_coordinate_map

# Define a single 2D Gaussian model.
p_init_gauss2D = models.Gaussian2D(x_mean=xx[size, size], y_mean=yy[size, size],
                                   amplitude=0.005,
                                   x_stddev=0.0001 * u.arcsec, y_stddev=0.0001 * u.arcsec)


# And fit with the Levenberg-Marquardt algorithm and least squares statistic.
fit_p = fitting.LevMarLSQFitter()

# NEEDS TO BE FIXED. should be able to use with_fill_value for projections
# fill value is NOT working for Projection
# mom0_sub.with_fill_value(fill_value=0.0).filled_data[:]

#moment0_cutout_quant = moment0_cutout.quantity
moment0_cutout_quant = moment0_cutout.value
moment0_cutout_quant[np.isnan(moment0_cutout_quant)] = 0.0

with warnings.catch_warnings():
    # Ignore model linearity warning from the fitter
    warnings.simplefilter('ignore')
    p_gauss2D = fit_p(p_init_gauss2D, xx, yy, moment0_cutout_quant)
    


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
plt.imshow(moment0_cutout.value - p_gauss2D(xx, yy).value, origin='lower', cmap='inferno', vmin=-0.9, vmax=0.9)
plt.colorbar()

plt.tight_layout()

