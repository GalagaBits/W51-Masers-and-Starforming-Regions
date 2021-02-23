import numpy as np
import pyspeckit
import pyspeckit.spectrum.fitters
from spectral_cube import SpectralCube

# Save metadata
field='W51North'
spw=8
directory1 = '/orange/adamginsburg/w51/vla/19A-254/derod/'

maser_channel = 318

spwfile = field+'_spw_'+str(spw)+'.fits'
directory = directory1+spwfile

# Position of maser in pixels (use ds9).
x = 486
y = 485

# Cuts spectra to the appropriate length for gaussian fits.

maser_channel1 = (maser_channel) - 1
mc1 = maser_channel1 - 40
mc2 = maser_channel1 + 40

cube = SpectralCube.read(directory)
extracted_spectrum = cube[:,y,x]

sp = pyspeckit.Spectrum(data=extracted_spectrum[mc1:mc2],
                        xarr=extracted_spectrum.spectral_axis[mc1:mc2])
#sp.xarr.convert_to_unit('km/s')
sp.plotter()

sp.specfit(fittype='gaussian')
sp.specfit(fittype='gaussian')

print(directory1+field+'_'+str(spw)+'_gaussian.png')

cube[maser_channel, :, :].quicklook()

print(sp.specfit(fittype='gaussian'))

answer = None 
while answer not in ("yes", "no"): 
    answer = input("Enter yes or no: ") 
    if answer == "yes": 
         sp.plotter.savefig(directory1+field+'_spw_'+str(spw)+'_gaussian.png')
         print('File saved.')

