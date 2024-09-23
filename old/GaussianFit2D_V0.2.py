#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 19 17:34:50 2021

@author: dealderod
"""

import astropy
import numpy as np
from spectral_cube import SpectralCube
from astropy import units as u
#Import data

cube = SpectralCube.read('W51North_spw_28.image', format='casa_image')

cube = cube.with_spectral_unit(u.km / u.s, velocity_convention='radio', rest_value= 25.715141 * u.GHz)
x, y = 504, 484
cube = cube[:,x,y]

x_mean = cube.mean(axis=-1)
y_mean = cube.mean(axis=0)

x_stddev = cube.mad_std(axis=-1)
y_stddev = cube.mad_std(axis=0)

#Pixel location (x,y)

print(cube.mean(axis=0))
print(cube.mean(axis=-1))

def gaussian2D(xvals,yvals,x,y,x_mean,y_mean,A):
    return A*np.exp(-(((xvals-x)^2)/(2*x_mean^2)) + (((yvals-y)^2)/(2*y_mean^2)))

import pylab as pl


#pl.plot(xval, gaussian(xval, a_77_2, FWHM_77_2, rest_77_2),color='blue')