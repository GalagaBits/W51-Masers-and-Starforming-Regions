import pyspeckit
from spectral_cube import SpectralCube

#Import data
cube = SpectralCube.read('/orange/adamginsburg/w51/vla/19A-254/derod/W51North_spw_53_corrected2.image', format='casa_image')

import pylab as pl
pl.style.use('default')

#Pixel location (x,y)
x, y = 637, 358

sp = pyspeckit.Spectrum(xarr=cube.spectral_axis, data=cube[:, y, x])
sp.plotter(color='black')

#Saving methods
field='W51North'
spw=53
directory1 = '/orange/adamginsburg/w51/vla/19A-254/derod/gaussianplots/'

#Save gaussianfit
def saveplotfig_gaussianfit():
    answer = None 
    while answer not in ("yes", "no"): 
        answer = input("Save the gaussian fit (yes or no)?") 
        if answer == "yes": 
            sp.plotter.savefig(directory1+field+'_spw_'+str(spw)+'_GaussianFit1D_V1.3.png')
            print('File saved.')
        
#Save residual        
def saveplotfig_residual():
    answer = None 
    while answer not in ("yes", "no"): 
        answer = input("Save the residual (yes or no)?") 
        if answer == "yes": 
            pl.savefig(directory1+field+'_spw_'+str(spw)+'_GaussianFit1D_residual_V1.3.png')
            print('File saved.')
            
from astropy import units as u

#Plot the gassiuan fits in to each peak in the spectra
fig = pl.figure(figsize=(12,8))
sp.plotter(color='black', xmin=19.2145*u.GHz, xmax=19.2155*u.GHz, figure=fig)

sp.baseline(interactive=True, subtract=False)

sp.specfit(interactive=True)


saveplotfig_gaussianfit()


