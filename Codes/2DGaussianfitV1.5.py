import warnings
import numpy as np
import matplotlib.pyplot as plt
from spectral_cube import SpectralCube, Projection
from astropy.io import fits
from astropy.modeling import models, fitting
from astropy.convolution import convolve_models
from astropy import units as u

from astropy.visualization import quantity_support

#Accessing Cube Data
cube = SpectralCube.read("/orange/adamginsburg/w51/vla/19A-254/derod/sample18.fits")
print(cube.max,"Cube Max)")
print(cube.unit,"Unit Flux")
cube.beam_threshold = 0.5

# maser_channel = 579

x, y = 3227, 5226
size = 1000


p_init_gauss2D = models.Gaussian2D(x_mean=290.9170756769094 * u.deg, y_mean=14.518231385948818 * u.deg, amplitude=10000 * (u.Jy/u.beam),
                                   x_stddev=0.02 * u.arcsec, y_stddev=0.02 * u.arcsec)
print(p_init_gauss2D,"Printed models.Gaussian2D")

yy, xx = cube.spatial_coordinate_map
print(yy,xx)

fit_p = fitting.LevMarLSQFitter()

with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    p_gauss2D = fit_p(p_init_gauss2D, xx, yy, cube[0,:,:])

print(p_gauss2D)


plt.figure(figsize=(18, 6))
plt.title("W51North NH3 (7,7)")

plt.subplot(1, 3, 2)
plt.title("Model")
plt.imshow(p_gauss2D(xx, yy).value, origin='lower', cmap='inferno')
plt.colorbar("Brightness (Jy)")
plt.xlabel("x Location (degrees)")
plt.ylabel("y Location (degrees)")

plt.tight_layout()

plt.show()
