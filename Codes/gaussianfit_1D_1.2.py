import pyspeckit
from spectral_cube import SpectralCube

cube = SpectralCube.read('/orange/adamginsburg/w51/vla/19A-254/derod/', use_dask=True)

sp.plotter()
