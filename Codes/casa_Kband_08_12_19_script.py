import casatools

tclean(vis='/orange/adamginsburg/w51/vla/19A-254/derod/19A-254_2019_08_12_T06_24_59.084/19A-254.sb36820841.eb37052030.58707.15631115741.ms',
    field='3',spw='3~10',imagename='W51North_spw_03~10_2019_08_12_010',
    cell='0.1arcsec',niter=100000,specmode='cube',imsize=800,pblimit=0.1,interactive=False,
    robust=0,chanchunks=-1,savemodel='modelcolumn',threshold='0.0400435Jy',outframe='LSRK',
    weighting='briggs')