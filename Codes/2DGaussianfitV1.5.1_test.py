import warnings
import numpy as np
import matplotlib.pyplot as plt
from spectral_cube import SpectralCube, Projection
from astropy.io import fits
from astropy.modeling import models, fitting
from astropy.convolution import convolve_models
from astropy import units as u

from astropy.visualization import quantity_support

import pylab as pl

# Saving methods
field = 'W51North'
spw = 8
directory1 = '/orange/adamginsburg/w51/vla/19A-254/derod/W51-Masers-and-Starforming-Regions/gaussianplots/Gaussian2D_plots/'


def saveplotfig_gaussianfit():
    answer = None
    while answer not in ("yes", "no"):
        answer = input("Save the gaussian fit (yes or no)?")
        if answer == "yes":
            pl.savefig(directory1 + field + '_spw_' + str(spw) + '_GaussianFit2D_V1.0.png')
            print('File saved.')


# Accessing Cube Data
cube = SpectralCube.read("/orange/adamginsburg/w51/vla/19A-254/derod/W51North_spw_8_corrected13_1.image",
                         format='casa_image')
# cube = cube[,:,:]
print(cube.max, "Cube Max)")
print(cube.unit, "Unit Flux")
cube.beam_threshold = 0.5

# maser_channel = 579

x, y = 637, 358
size = 20

cube = cube[84, :, :]
cube.quicklook()
# cube = cube.max(axis=0)
cube_cutout = cube[y - size:y + size, x - size:x + size]
cube_cutout.quicklook()

'''
p_init_gauss2D = models.Gaussian2D(x_mean=290.915 * u.deg, y_mean=14.517248 * u.deg, amplitude=1000 * (u.Jy / u.beam),
                                   x_stddev=1.11111e-5 * u.degree, y_stddev=1.11111e-5 * u.degree, theta=(np.pi / 2))

print(p_init_gauss2D, "Printed models.Gaussian2D")

yy, xx = cube_cutout.spatial_coordinate_map

fit_p = fitting.LevMarLSQFitter()

cube_cutout_quant = cube_cutout
cube_cutout_quant[np.isnan(cube_cutout_quant)] = 0.0

# p_init_gauss2D.theta.fixed = True

with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    p_gauss2D = fit_p(p_init_gauss2D, xx, yy, cube_cutout_quant)

# p_gauss2D_avg = (p_gauss2D.x_stddev[0], p_gauss2D.y_stddev[0])
# p_gauss2D.x_stddev[0], p_gauss2D.y_stddev[0] = p_gauss2D_avg

fitted_x_stddev = p_gauss2D.x_stddev
fitted_y_stddev = p_gauss2D.y_stddev

p_gauss2D_avg = (fitted_x_stddev + fitted_y_stddev) / 2
p_gauss2D.x_stddev = p_gauss2D_avg
p_gauss2D.y_stddev = p_gauss2D_avg

print("\n param_cov")

print(fit_p.fit_info)

cov = np.diag(fit_p.fit_info['param_cov'])
errors = np.sqrt(cov)

print("The errors are:\n", errors)

amplitude_error = np.format_float_scientific(errors[0], precision=9)
x_mean_error = np.format_float_scientific(errors[1], precision=9)
y_mean_error = np.format_float_scientific(errors[2], precision=9)
x_stddev_error = np.format_float_scientific(errors[3], precision=9)
y_stddev_error = np.format_float_scientific(errors[4], precision=9)
theta_error = np.format_float_scientific(errors[5], precision=9)

amplitude = np.format_float_scientific(p_gauss2D.amplitude[0], precision=9)
x_mean = np.format_float_scientific(p_gauss2D.x_mean[0], precision=9)
y_mean = np.format_float_scientific(p_gauss2D.y_mean[0], precision=9)
x_stddev = np.format_float_scientific(p_gauss2D.x_stddev[0], precision=9)
y_stddev = np.format_float_scientific(p_gauss2D.y_stddev[0], precision=9)
theta = np.format_float_scientific(p_gauss2D.theta[0], precision=9)

print(p_gauss2D.amplitude)
print(p_gauss2D.x_mean)
print(p_gauss2D.y_mean)
print(p_gauss2D.x_stddev)
print(p_gauss2D.y_stddev)

print(p_gauss2D)

plt.figure(figsize=(18, 6))
plt.suptitle("W51North NH3 (7,6) | Peak 3")
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

plt.subplot(1, 3, 3)
plt.title("Residual")
plt.imshow(cube_cutout.value - p_gauss2D(xx, yy).value, origin='lower', cmap='inferno')
plt.colorbar(label='S (Jy)')
plt.xlabel("x (degrees)")
plt.ylabel("y (degrees)")

plt.subplot(1, 3, 2)
plt.text(1, 2,
         "A(0) = " + str(amplitude) + " ± " + str(amplitude_error) + " Jy\n" + "x(0) = " + str(x_mean) + " ± " + str(
             x_mean_error) +
         " deg\n" + "y(0) = " + str(y_mean) + " ± " + str(y_mean_error) + " deg\n" + "σ_x(0) = " + str(
             x_stddev) + " ± " + str(x_stddev_error) +
         " deg\n" + "σ_y(0) = " + str(y_stddev) + " ± " + str(y_stddev_error) + " deg\n" + "θ(0) = " + str(
             theta) + " ± " + str(theta_error) + " rad",
         style='italic',
         bbox={'facecolor': 'black', 'alpha': 0.9, 'pad': 10}, color="white",
         horizontalalignment='left')

plt.tight_layout()

saveplotfig_gaussianfit()

plt.show()

'''

plt.show()