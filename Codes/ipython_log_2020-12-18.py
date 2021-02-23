########################################################
# Started Logging At: 2020-12-18 16:46:10
########################################################
########################################################
# # Started Logging At: 2020-12-18 16:46:11
########################################################
import pyspeckit
from spectral_cube import SpectralCube
cube = SpectralCube.read('../W51North_spw_6.fits')
x = 486
y = 485
spectrum = cube[:,y,x]
sp = pyspeckit.Spectrum(data=extracted_spectrum,
                        xarr=extracted_spectrum.spectral_axis)
extracted_spectrum = cube[:,y,x]
sp = pyspeckit.Spectrum(data=extracted_spectrum,
                        xarr=extracted_spectrum.spectral_axis)
sp.plotter()
peak_intensity = cube.max(axis=0)
peak_intensity = cube.max(axis=0, allow_huge=True)
peak_intensity = cube.max(axis=0, allow_huge_operations=True)
cube
#[Out]# VaryingResolutionSpectralCube with shape=(512, 1000, 1000) and unit=Jy / beam:
#[Out]#  n_x:   1000  type_x: RA---SIN  unit_x: deg    range:   290.909691 deg:  290.924024 deg
#[Out]#  n_y:   1000  type_y: DEC--SIN  unit_y: deg    range:    14.511228 deg:   14.525103 deg
#[Out]#  n_s:    512  type_s: FREQ      unit_s: Hz     range: 22916326611.850 Hz:22924310656.192 Hz
cube = SpectralCube.read('../W51North_spw_6.fits', use_dask=True)
cube
#[Out]# DaskVaryingResolutionSpectralCube with shape=(512, 1000, 1000) and unit=Jy / beam and chunk size (256, 250, 250):
#[Out]#  n_x:   1000  type_x: RA---SIN  unit_x: deg    range:   290.909691 deg:  290.924024 deg
#[Out]#  n_y:   1000  type_y: DEC--SIN  unit_y: deg    range:    14.511228 deg:   14.525103 deg
#[Out]#  n_s:    512  type_s: FREQ      unit_s: Hz     range: 22916326611.850 Hz:22924310656.192 Hz
peak_intensity = cube.max(axis=0)
pl.imshow(peak_intensity)
import pylab as pl
pl.imshow(peak_intensity)
#[Out]# <matplotlib.image.AxesImage at 0x2ae6a4b835c0>
pl.imshow(peak_intensity.value)
#[Out]# <matplotlib.image.AxesImage at 0x2ae6a4d3c320>
pl.imshow(peak_intensity.value, norm=visualization.simple_norm(peak_intensity.value, min_percent=1, max_percent=99.5))
pl.imshow(peak_intensity.value, norm=visualization.simple_norm(peak_intensity.value, min_percent=1, max_percent=99.5))
import pylab as pl
from astropy import visualization
pl.imshow(peak_intensity.value, norm=visualization.simple_norm(peak_intensity.value, min_percent=1, max_percent=99.5))
#[Out]# <matplotlib.image.AxesImage at 0x2ae6a6212898>
pl.imshow(peak_intensity.value, norm=visualization.simple_norm(peak_intensity.value, min_percent=1, max_percent=99.5))
pl.axis([400,600,400,600])
#[Out]# [400, 600, 400, 600]
pl.imshow(peak_intensity.value, norm=visualization.simple_norm(peak_intensity.value, min_percent=1, max_percent=99.5))
pl.axis([525,575,450,500])
#[Out]# [525, 575, 450, 500]
pl.imshow(peak_intensity.value, norm=visualization.simple_norm(peak_intensity.value, min_percent=1, max_percent=99.5))
pl.axis([475,525,450,500])
#[Out]# [475, 525, 450, 500]
pl.imshow(peak_intensity.value, norm=visualization.simple_norm(peak_intensity.value, min_percent=1, max_percent=99.5))
pl.axis([400,605,450,500])
#[Out]# [400, 605, 450, 500]
pl.imshow(peak_intensity.value, norm=visualization.simple_norm(peak_intensity.value, min_percent=1, max_percent=99.5))
pl.axis([400,605,400,500])
#[Out]# [400, 605, 400, 500]
pl.imshow(peak_intensity.value, norm=visualization.simple_norm(peak_intensity.value, min_percent=1, max_percent=99.5))
pl.axis([400,605,500,600])
#[Out]# [400, 605, 500, 600]
pl.imshow(peak_intensity.value, norm=visualization.simple_norm(peak_intensity.value, min_percent=1, max_percent=99.5))
pl.axis([450,500,450,550])
#[Out]# [450, 500, 450, 550]
pl.imshow(peak_intensity.value, norm=visualization.simple_norm(peak_intensity.value, min_percent=1, max_percent=99.5))
pl.axis([450,550,450,550])
#[Out]# [450, 550, 450, 550]
pl.imshow(peak_intensity.value, norm=visualization.simple_norm(peak_intensity.value, min_percent=1, max_percent=99.5))
pl.axis([450,550,450,550])
pl.plot(x,y,color='w')
#[Out]# [<matplotlib.lines.Line2D at 0x2ae6a6646780>]
pl.imshow(peak_intensity.value, norm=visualization.simple_norm(peak_intensity.value, min_percent=1, max_percent=99.5))
pl.axis([450,550,450,550])
pl.plot(x,y,color='w',marker='x')
#[Out]# [<matplotlib.lines.Line2D at 0x2ae6a66fd048>]
pl.imshow(peak_intensity.value, norm=visualization.simple_norm(peak_intensity.value, min_percent=1, max_percent=99.5))
pl.axis([450,550,450,550])
y = 505
x = 482
pl.plot(x,y,color='w',marker='x')
#[Out]# [<matplotlib.lines.Line2D at 0x2ae6a672f048>]
pl.imshow(peak_intensity.value, norm=visualization.simple_norm(peak_intensity.value, min_percent=1, max_percent=99.5))
pl.axis([470,505,480,520])
y = 505
x = 482
pl.plot(x,y,color='w',marker='x')
#[Out]# [<matplotlib.lines.Line2D at 0x2ae6a694aba8>]
pl.imshow(peak_intensity.value, norm=visualization.simple_norm(peak_intensity.value, min_percent=1, max_percent=99.5))
pl.axis([470,505,480,520])
y = 503
x = 484
pl.plot(x,y,color='w',marker='x')
#[Out]# [<matplotlib.lines.Line2D at 0x2ae6a6a0e6d8>]
pl.imshow(peak_intensity.value, norm=visualization.simple_norm(peak_intensity.value, min_percent=1, max_percent=99.5))
pl.axis([470,505,480,520])
y = 504
x = 483
pl.plot(x,y,color='w',marker='x')
#[Out]# [<matplotlib.lines.Line2D at 0x2ae6a5fde160>]
pl.imshow(peak_intensity.value, norm=visualization.simple_norm(peak_intensity.value, min_percent=1, max_percent=99.5))
pl.axis([470,505,480,520])
y = 504
x = 483
pl.plot(x,y,color='k',marker='x')
#[Out]# [<matplotlib.lines.Line2D at 0x2ae6a607fc88>]
extracted_spectrum = cube[:,y,x]
sp = pyspeckit.Spectrum(data=extracted_spectrum,
                        xarr=extracted_spectrum.spectral_axis)
extracted_spectrum = cube[:,y,x]
sp = pyspeckit.Spectrum(data=extracted_spectrum,
                        xarr=extracted_spectrum.spectral_axis)
sp.plotter()
sp.plotter(xmin=22.92*u.GHz, xmax=22.922*u.GHz)
from astropy import units as u
sp.plotter(xmin=22.92*u.GHz, xmax=22.922*u.GHz)
sp.plotter(xmin=22.925*u.GHz, xmax=22.922*u.GHz)
sp.plotter(xmin=22.9215*u.GHz, xmax=22.922*u.GHz)
sp.plotter(xmin=22.9205*u.GHz, xmax=22.922*u.GHz)
sp.specfit()
sp.plotter(xmin=22.9205*u.GHz, xmax=22.922*u.GHz)
sp.specfit()
sp.plotter(xmin=22.9205*u.GHz, xmax=22.922*u.GHz)
sp.specfit(guesses=[0.05,22.9210e9,5e4, 0.05, 22.92115e9, 5e4, 0.05, 22.913e9, 5e4])
sp.plotter(xmin=22.9205*u.GHz, xmax=22.922*u.GHz)
sp.specfit(guesses=[0.05,22.9209e9, 2e4, 0.05, 22.92115e9, 2e4, 0.05, 22.913e9, 2e4])
get_ipython().run_line_magic('matplotlib', 'notebook')
sp.plotter(xmin=22.9205*u.GHz, xmax=22.922*u.GHz)
sp.specfit(guesses=[0.05,22.9209e9, 2e4, 0.05, 22.92115e9, 2e4, 0.05, 22.913e9, 2e4])
get_ipython().run_line_magic('matplotlib', 'widget')
sp.plotter(xmin=22.9205*u.GHz, xmax=22.922*u.GHz)
sp.specfit(guesses=[0.05,22.9209e9, 2e4, 0.05, 22.92115e9, 2e4, 0.05, 22.913e9, 2e4])
get_ipython().run_line_magic('matplotlib', 'inline')
sp.plotter(xmin=22.9205*u.GHz, xmax=22.922*u.GHz)
sp.specfit(guesses=[0.05,22.9209e9, 2e4, 0.05, 22.92115e9, 2e4, 0.05, 22.913e9, 2e4])
sp.plotter(xmin=22.9205*u.GHz, xmax=22.922*u.GHz)
sp.specfit(guesses=[0.05,22.92095e9, 2e4, 0.05, 22.92115e9, 2e4, 0.05, 22.913e9, 2e4])
sp.plotter(xmin=22.9205*u.GHz, xmax=22.922*u.GHz)
sp.specfit(guesses=[0.05,22.92095e9, 2e4, 0.05, 22.92115e9, 2e4, 0.05, 22.9135e9, 2e4])
sp.plotter(xmin=22.9205*u.GHz, xmax=22.922*u.GHz)
sp.specfit()
sp.plotter(xmin=22.9205*u.GHz, xmax=22.922*u.GHz)
sp.specfit(guesses=[0.05,22.9209e9, 2e4, 0.05, 22.92115e9, 2e4, 0.05, 22.913e9, 2e4])
get_ipython().run_line_magic('matplotlib', 'inline')
sp.plotter(xmin=22.9205*u.GHz, xmax=22.9215*u.GHz)
sp.specfit(guesses=[0.05,22.9209e9, 2e4, 0.05, 22.92115e9, 2e4, 0.05, 22.913e9, 2e4])
sp.plotter(xmin=22.9208*u.GHz, xmax=22.9215*u.GHz)
sp.specfit(guesses=[0.05,22.9209e9, 2e4, 0.05, 22.92115e9, 2e4, 0.05, 22.913e9, 2e4])
sp.plotter(xmin=22.9208*u.GHz, xmax=22.9215*u.GHz)
sp.specfit(guesses=[0.05, 22.9210e9, 2e4, 0.05, 22.92115e9, 2e4, 0.1, 22.9129e9, 2e4])
sp.plotter(xmin=22.9208*u.GHz, xmax=22.9215*u.GHz)
sp.specfit(guesses=[0.0566, 22.9212e9, 1.15e5, 0.05, 22.92115e9, 2e4, 0.1, 22.9129e9, 2e4])
sp.plotter(xmin=22.9208*u.GHz, xmax=22.9215*u.GHz)
sp.specfit(guesses=[0.0566, 22.9212e9, 1.15e5, 0.05, 22.92110e9, 2e4, 0.1, 22.9129e9, 2e4])
sp.xarr.convert_to_unit(u.km/u.s)
cube = SpectralCube.read('../W51North_spw_6.fits', use_dask=True)
cube
#[Out]# DaskVaryingResolutionSpectralCube with shape=(512, 1000, 1000) and unit=Jy / beam and chunk size (256, 250, 250):
#[Out]#  n_x:   1000  type_x: RA---SIN  unit_x: deg    range:   290.909691 deg:  290.924024 deg
#[Out]#  n_y:   1000  type_y: DEC--SIN  unit_y: deg    range:    14.511228 deg:   14.525103 deg
#[Out]#  n_s:    512  type_s: FREQ      unit_s: Hz     range: 22916326611.850 Hz:22924310656.192 Hz
sp.xarr.convert_to_unit(u.km/u.s, rest_value=22.92494000*u.GHz)
sp.xarr.convert_to_unit(u.km/u.s, refX=22.92494000*u.GHz)
get_ipython().run_line_magic('pinfo', 'sp.xarr.convert_to_unit')
sp.xarr.refX = 22.92494000*u.GHz
sp.xarr.convert_to_unit(u.km/u.s)
sp.xarr.velocity_convention
sp.xarr.velocity_convention = 'radio'
sp.xarr.convert_to_unit(u.km/u.s)
sp.plotter()
sp.plotter(xmin=40, xmax=60)
sp.plotter(xmin=45, xmax=53)
sp.plotter(xmin=45, xmax=53)
sp.specfit(guesses=[0.09, 47.8, 0.1, 0.06, 49.5, 0.1, 0.05, 51.5, 0.01])
sp.specfit.parinfo
#[Out]# [Param #0   AMPLITUDE0 =    0.0773156 +/-       0.0148047 ,
#[Out]#  Param #1       SHIFT0 =      47.6831 +/-       0.0591147 ,
#[Out]#  Param #2       WIDTH0 =     0.292679 +/-       0.0704842   Range:   [0,inf),
#[Out]#  Param #3   AMPLITUDE1 =    0.0453511 +/-       0.0057263 ,
#[Out]#  Param #4       SHIFT1 =      50.0001 +/-        0.360758 ,
#[Out]#  Param #5       WIDTH1 =      1.79169 +/-        0.435567   Range:   [0,inf),
#[Out]#  Param #6   AMPLITUDE2 =   -0.0164288 +/-       0.0159103 ,
#[Out]#  Param #7       SHIFT2 =      52.6956 +/-        0.298531 ,
#[Out]#  Param #8       WIDTH2 =     0.276786 +/-        0.375059   Range:   [0,inf)]
sp.specfit.parinfo[:3].fixed = True
for par in sp.specfit.parinfo[:3]:
    par.fixed=True
sp.plotter(xmin=45, xmax=53)
sp.specfit.refit()
sp.plotter(xmin=45, xmax=53)
sp.specfit(guesses=[0.09, 47.8, 0.1, 0.06, 49.5, 0.1, 0.05, 51.5, 0.1])
15.625 / 22e9
#[Out]# 7.102272727272728e-10
15.625 / 22e9 * 3e5
#[Out]# 0.00021306818181818183
15.625*1e3 / 22e9 * 3e5
#[Out]# 0.2130681818181818
