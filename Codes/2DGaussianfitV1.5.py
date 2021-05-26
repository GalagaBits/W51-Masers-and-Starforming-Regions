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
size = 250

cube1 = cube[0,:,:]

max_2d_cube = cube.max(axis=0)
cube_cutout = max_2d_cube[y - size:y + size, x - size:x + size]
cube_cutout.quicklook()

p_init_gauss2D = models.Gaussian2D(x_mean=290.9170756769094 * u.deg, y_mean=14.518231385948818 * u.deg, amplitude=10000 * (u.Jy/u.beam),
                                   x_stddev=0.02 * u.arcsec, y_stddev=0.02 * u.arcsec)

print(p_init_gauss2D,"Printed models.Gaussian2D")

yy, xx = cube_cutout.spatial_coordinate_map

fit_p = fitting.LevMarLSQFitter()

cube_cutout_quant = cube_cutout
cube_cutout_quant[np.isnan(cube_cutout_quant)] = 0.0

with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    p_gauss2D = fit_p(p_init_gauss2D, xx, yy, cube_cutout_quant)

print(p_gauss2D.amplitude)
print(p_gauss2D.x_mean)
print(p_gauss2D.y_mean)
print(p_gauss2D.x_stddev)
print(p_gauss2D.y_stddev)

print(p_gauss2D)

plt.figure(figsize=(18, 6))
plt.suptitle("W51North NH3 (7,7)")
plt.subplot(1, 3, 1)
plt.title("Image")
plt.imshow(cube_cutout.value, origin='lower', cmap='inferno')
plt.colorbar()
plt.xlabel("x (degrees)")
plt.ylabel("y (degrees)")
plt.subplot(1, 3, 2)
plt.title("Model")
plt.imshow(p_gauss2D(xx, yy).value, origin='lower', cmap='inferno')
plt.colorbar()
plt.xlabel("x (degrees)")
plt.ylabel("y (degrees)")
plt.text(20, 360, "A(0) = "+str(p_gauss2D.amplitude[0])+" Jy\n"+"x(0) = "+str(p_gauss2D.x_mean[0])+
        " deg\n"+"y(0) = "+str(p_gauss2D.y_mean[0])+ " deg\n" + "σ_x(0) = "+str(p_gauss2D.x_stddev[0])+
        " deg\n"+"σ_y(0) = "+str(p_gauss2D.y_stddev[0])+" deg\n" +"θ(0) = "+str(p_gauss2D.theta[0])+" rad",
        style='italic',
        bbox={'facecolor': 'white', 'alpha': 0.9, 'pad': 10},
        horizontalalignment='left')
plt.subplot(1, 3, 3)
plt.title("Residual")
plt.imshow(cube_cutout.value - p_gauss2D(xx, yy).value, origin='lower', cmap='inferno')
plt.colorbar(label='S (Jy)')
plt.xlabel("x (degrees)")
plt.ylabel("y (degrees)")
plt.tight_layout()
plt.show()