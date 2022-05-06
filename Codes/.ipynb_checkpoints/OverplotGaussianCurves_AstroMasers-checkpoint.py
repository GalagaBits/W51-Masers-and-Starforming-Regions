import pyspeckit
from spectral_cube import SpectralCube

#Import data
cube = SpectralCube.read('W51North_spw_28.image', format='casa_image')

import pylab as pl

cube = cube.with_spectral_unit(u.km / u.s, velocity_convention='radio', rest_value= 25.715141 * u.GHz)

#Pixel location (x,y)
x, y = 504, 484
sp = pyspeckit.Spectrum(xarr=cube.spectral_axis, data=cube[:, x, y])


import numpy as np
import matplotlib.pyplot as plt

a_77_1      = 0.20
FWHM_77_1    = 1.38
rest_77_1 = 47.29#LSR in km s-1

a_77_2     = 0.18
FWHM_77_2    = 1.23
rest_77_2 = 49.92

a_77_3    = 1.70
FWHM_77_3   = 1.46
rest_77_3 = 45.60#LSR in km s-1

a_77_4     = 2.48
FWHM_77_4    = 1.27
rest_77_4 = 46.08


def gaussian(xval, a, FWHM, rest_f):
    return (a*(np.exp(-np.log(2)*(((xval-rest_f)**2)/((0.5*FWHM)**2)))))

xval = cube.spectral_axis.value

fig = pl.figure(figsize=(8,6))

sp.plotter(color='black', xmin=40, xmax=60, ymax=3, figure=fig)

#sp.plotter(color='black')
pl.title('W51North Deal v. Henkel', fontname="Times New Roman")
pl.xlabel('Velocity (km s-1)', fontname="Times New Roman")
pl.ylabel('S (Jy)',fontname="Times New Roman")

pl.plot(cube.spectral_axis.value, np.array(cube[:, x, y].value), label="Deal")

pl.plot(xval, gaussian(xval, a_77_1, FWHM_77_1, rest_77_1),color='blue', label="Henkel, 2012")
pl.plot(xval, gaussian(xval, a_77_2, FWHM_77_2, rest_77_2),color='blue')

pl.plot(xval, gaussian(xval, a_77_3, FWHM_77_3, rest_77_3),color='limegreen', label="Henkel, 2008")
pl.plot(xval, gaussian(xval, a_77_4, FWHM_77_4, rest_77_4),color='limegreen')

pl.legend()
