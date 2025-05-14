"""
============
UCHII Fitter
============

Fit a free-free spectrum to an SED.

.. moduleauthor:: Adam Ginsburg <adam.g.ginsburg@gmail.com>

"""
import numpy as np
try:
    from scipy import optimize
except ImportError:
    print("scipy not installed: UCHIIfitter may fail")
from pyspeckit import mpfit
from astropy import units as u
from astropy import constants
from astropy.modeling.models import custom_model
from astropy.modeling.fitting import LevMarLSQFitter
from astropy import log

import radio_beam
from dust_emissivity import dust

def val_with_unit(var, unit):
    if hasattr(var, 'unit'):
        return var.to(unit).value
    else:
        return var

def with_unit(var, unit):
    if hasattr(var, 'unit'):
        return var.to(unit)
    else:
        return var*unit


#kb = 1.38e-16
#c=3e10
#mu = 1.4
#mh = 1.67e-24
#msun = 1.9889e33
#pc = 3.08568e18      # cm
#au = 1.496e13        # cm
#msun = 1.99e33       # g

unitfactor={'mJy':1e-26,'Jy':1e-23,'cgs':1.0}
freqfactor={'GHz':1e9,'Hz':1.0}
muh = 2.8

default_te = 8500*u.K
default_freq = 1*u.GHz
# using the value from http://www.cv.nrao.edu/~sransom/web/Ch4.html
alpha_b = 3e-13*u.cm**3*u.s**-1
# emission measure unit
emu = u.cm**-6*u.pc

def tnu(Te, nu, EM):
    """
    Optical Depth for a given electron temperature, frequency, and emission
    measure

    Parameters
    ----------
    Te : K
        excitation temperature
    nu : GHz
        frequency in GHz
    EM : pc cm^-6
        Emission Measure

    Calculates optical depth as a function of temperature, frequency, and
    emission measure from Rohlfs and Wilson 2000's eqns 9.33 and 9.34.

    (edit: later edition, this is equations 10.35 and 10.36)

    """
    Te = with_unit(Te, u.K)
    nu = with_unit(nu, u.GHz)
    EM = with_unit(EM, emu)

#    nu0 = .3045 * Te**-.643 * EM**.476

    # nu0  is the frequency above which the Gaunt factor is exactly one
    gff_lownu = (np.log(4.955e-2 * (nu/u.GHz).decompose().value**-1) + 1.5 * np.log(Te/u.K))  # <gff> Gaunt factor for free-free


    nu0 = ((Te/u.K)**1.5 * u.MHz).to(u.GHz)
    # solve gff_lownu = 1
    nu0 = 4.955e-2 / np.exp(1-1.5 * np.log(Te/u.K)) * u.GHz


    answer_lownu = (nu < nu0) * 3.014e-2 * (Te/u.K)**-1.5 * (nu/u.GHz)**-2 * (EM/emu) * gff_lownu
    answer_highnu = (nu > nu0) * 3.014e-2 * (Te/u.K)**-1.5 * (nu/u.GHz)**-2 * (EM/emu)

    tau = answer_lownu+answer_highnu
    ## altenhoff version
    #tau = 8.235e-2 * Te**-1.35 * nu**-2.1 * EM
    return tau

def tau(nu, EM, Te=default_te):
    """
    Another equation for optical depth, less explicit.
    This is some approximation from Rohlfs & Wilson, but the above is probably
    better to use.
    """
    return (3.28e-7 * (Te/(1e4*u.K))**-1.35 * (nu/u.GHz)**-2.1 *
            (EM/(u.cm**-6*u.pc)))


def Inu(nu, tau, Te, I0=0):
    """
    Calculates flux for a given optical depth, frequency, and temperature
    assuming Rayleigh-Jeans

    nu - frequency in Hz
    tau - optical depth
    Te - excitation temperature (K)

    Parameters
    ----------
    I0 : float
        Scale factor.  If zero, will be determined.  If not, will still be
        determined.
    """
    nu = with_unit(nu, u.GHz)
    Te = with_unit(Te, u.K)

    if I0==0 and isinstance(nu,np.ndarray):
        whtau1 = np.argmin(np.abs(tau-1))
        nutau1 = nu[whtau1]
        taufactor = 1
    else:
        nutau1 = nu
        taufactor = tau
        """ assumes I0 is set"""
    I0 = 2 * constants.k_B * Te * nutau1**2 / constants.c**2 * taufactor
    thin = (tau < 1) * np.exp(1-tau) * I0
    thick = 2 * constants.k_B * Te * (nu * (tau > 1))**2 / constants.c**2
    #return I0 * np.exp(1-tau) # WRONG.
    return thin+thick

@custom_model
def inufit(nu, em=1e7, normfac=1e-10, Te=default_te.value, nu0=1.5):
    """
    Computes the expected intensity as a function of frequency
    for a given emission measure and normalization factor
    nu - array of frequencies (array)
    em - emission measure (float)
    normfac - normalization factor (float)
            - 1/solid angle of source.  1000 AU at 1 kpc = 206265.

    Units: mJy
    """

    nu = with_unit(nu, u.GHz)
    nu0 = with_unit(nu0, u.GHz)
    em = with_unit(em, emu)
    Te = with_unit(Te, u.K)

    I0 = 2 * constants.k_B * Te * nu0**2 / constants.c**2
    tau = tnu(Te, nu, em)  # tnu takes GHz
    model_intensity = Inu(nu, tau, Te, I0=I0)
    model_norm = normfac * model_intensity
    return model_norm.to(u.mJy).value


@custom_model
def inufit_dust(nu, em=1e7, normfac=1e-10, alpha=3.75, normfac2=1e2,
                Te=default_te.value, nu0=1.5):
    """
    inufit with dust added

    Dust spectral index is given by alpha
    """

    nu = with_unit(nu, u.GHz)
    nu0 = with_unit(nu0, u.GHz)
    em = with_unit(em, emu)
    Te = with_unit(Te, u.K)

    I0 = 2 * constants.k_B * Te * nu0**2 / constants.c**2
    model_intensity = Inu(nu, tnu(Te,nu,em), Te, I0=I0)

    dustem = normfac2 * u.mJy * (nu/nu0)**alpha

    model_norm = (normfac * model_intensity + dustem)

    return model_norm.to(u.mJy).value

@custom_model
def inufit_dustT(nu, em=1e7, normfac=1e-10, beta=1.75, normfac2=1e-5,
                 dustT=50, Te=default_te.value, nu0=1.5):
    """
    Parameters
    ----------
    nu : frequency
    """

    #print(f"em={em.item():0.2g} normfac={normfac.item():0.2g} beta={beta.item():0.2g} normfac2={normfac2.item():0.2g} dustT={dustT.item()}")

    nu = with_unit(nu, u.GHz)
    nu0 = with_unit(nu0, u.GHz)
    em = with_unit(em, emu)
    Te = with_unit(Te, u.K)
    dustT = with_unit(dustT, u.K)

    I0 = 2 * constants.k_B * Te * nu0**2 / constants.c**2
    model_intensity = Inu(nu,tnu(Te,nu,em),Te,I0=I0)
    #dustem = (2*constants.h*(nu)**(3+beta) / constants.c**2 *
    #          (np.exp(constants.h*nu*1e9/(constants.k_B*np.abs(dustT))) -
    #           1)**-1)
    dustem = dust.blackbody.modified_blackbody(nu=nu,
                                               temperature=dustT,
                                               column=u.Quantity(1e18, u.cm**-2),
                                               beta=beta) * u.sr
    #print(model_intensity.decompose())
    #print(dustem.decompose())
    model_norm = normfac * model_intensity + normfac2*dustem
    return model_norm.to(u.mJy).value


def mpfitfun(freq, flux, err=None, dust=False, dustT=False):
    """ wrapper around inufit to be passed into mpfit """

    flux = with_unit(flux, u.Jy)
    freq = with_unit(freq, u.GHz)
    if err is not None:
        err = with_unit(err, flux.unit)

    if dust:
        if err is None:
            def f(p,fjac=None): return [0,(flux-inufit_dust(*p)(freq).to(u.Jy)).decompose().value]
        else:
            def f(p,fjac=None): return [0,((flux-inufit_dust(*p)(freq).to(u.Jy))/err).decompose().value]
        return f
    elif dustT:
        if err is None:
            def f(p,fjac=None): return [0,(flux-inufit_dustT(*p)(freq).to(u.Jy)).decompose().value]
        else:
            def f(p,fjac=None): return [0,((flux-inufit_dustT(*p)(freq).to(u.Jy))/err).decompose().value]
        return f
    else:
        if err is None:
            def f(p,fjac=None): return [0,(flux-inufit(*p)(freq).to(u.Jy)).decompose().value]
        else:
            def f(p,fjac=None): return [0,((flux-inufit(*p)(freq).to(u.Jy))/err).decompose().value]
        return f

def emtau(freq, flux, err=None, EMguess=1e7*emu, Te=default_te, normfac=5e-6,
          quiet=1, dust=False, dustT=False, alpha=3.5, normfac2=1e-6,
          beta=1.75, use_mpfit=False, maxiter=500, **kwargs):
    """
    Returns emission measure & optical depth given radio continuum data points
    at frequency freq with flux density flux.

    return bestEM,nu(tau=1),chi^2
    """
    EMguess = with_unit(EMguess, emu)
    Te = with_unit(Te, u.K)
    flux = with_unit(flux, u.Jy)
    err = with_unit(err, u.Jy)

    ok = np.isfinite(freq) & np.isfinite(flux) & np.isfinite(err)

    guesses = [(EMguess/emu).decompose().value, normfac]
    if dust:
        guesses += [alpha, normfac2]
    elif dustT:
        guesses += [beta, normfac2, dustT.to(u.K).value]

    if use_mpfit:
        mp = mpfit.mpfit(mpfitfun(freq[ok], flux[ok], err[ok], dust=dust,
                                  dustT=bool(dustT)),
                         xall=guesses,
                         maxiter=maxiter,
                         quiet=quiet)
        mpp = mp.params
        mpperr = mp.perror
        chi2 = mp.fnorm
        bestEM = mpp[0]
        normfac = mpp[1]
    else:
        fitter = LevMarLSQFitter()
        model = (inufit_dustT if bool(dustT) else inufit_dust if dust else inufit)
        m_init = model(*guesses)

        m_init.Te.fixed = True
        m_init.nu0.fixed = True

        m_init.em.min = 1 # cannot be less than 1 cm^-6 pc
        m_init.normfac.min = 0
        if hasattr(m_init, 'beta'):
            m_init.beta.min = 1
            m_init.beta.max = 10
            m_init.normfac2.min = 0
        elif hasattr(m_init, 'alpha'):
            m_init.alpha.min = 0
            m_init.alpha.max = 12

        fitted = fitter(m_init,
                        freq[ok].to(u.GHz).value,
                        flux[ok].to(u.mJy).value,
                        weights=1/err[ok].to(u.mJy).value,
                        maxiter=maxiter
                       )
        mp = fitter
        mp.params = fitted.parameters
        mpp = fitted.parameters
        try:
            assert 'cov_x' in fitter.fit_info and fitter.fit_info['cov_x'] is not None
            mpperr = fitter.fit_info['cov_x'].diagonal()**0.5
        except AssertionError:
            mpperr = [np.nan] * len(guesses)
        chi2 = (((u.Quantity(fitted(freq[ok]), u.mJy) - flux[ok])/err[ok])**2).sum() / (ok.sum() - len(guesses))
        bestEM = mpp[0]
        normfac = mpp[1]
        log.info(fitter.fit_info['message'])

    # bestEM is unitless w/ emu units
    nu_tau = (((Te/u.K)**1.35 / (bestEM) / 8.235e-2)**(-1/2.1)).decompose()

    return with_unit(bestEM, emu), nu_tau, normfac, chi2, mp

class HIIregion(object):
    """
    An HII region has properties frequency, flux, and error, which must be
    numpy ndarrays of the same length
    """

    def __init__(self, nu, flux, fluxerr, fluxunit='mJy', frequnit='GHz',
                 beamsize_as2=0.25*u.arcsec**2, dist_kpc=1.0*u.kpc, resolved=False,
                 Te=default_te, **kwargs):
        order = np.argsort(np.asarray(nu))
        self.frequnit     = u.Unit(frequnit)
        self.fluxunit     = u.Unit(fluxunit)
        self.nu           = with_unit(np.asarray(nu)[order], self.frequnit)
        self.flux         = with_unit(np.asarray(flux)[order], self.fluxunit)
        self.fluxerr      = with_unit(np.asarray(fluxerr)[order], self.fluxunit)
        self.beamsize_as2 = beamsize_as2
        self.dist_kpc = dist_kpc
        self.resolved = resolved
        self.Te = Te

        self.fit(**kwargs)

    def fit(self, **kwargs):
        """  """
        self.em, self.nutau, self.normfac, self.chi2, self.mp = emtau(freq=self.nu,
                                                                      flux=self.flux,
                                                                      err=self.fluxerr,
                                                                      Te=self.Te,
                                                                      **kwargs)
        self.params = self.mp.params
        if 'dustT' in kwargs and kwargs['dustT']:
            self.beta = self.params[2]
            self.dustT = self.params[4]
            self.normfac2 = self.params[3]
        elif 'dust' in kwargs and kwargs['dust']:
            self.alpha = self.params[2]
            self.normfac2 = self.params[3]


    def loglogplot(self, numin=1.0*u.GHz, numax=10.0*u.GHz, plottitle='',
                   params=None, do_annotations=True, dust=False, dustT=False,
                   annotation_xpos=0.7,
                   **kwargs):
        import pylab as pl
        x = np.logspace(np.log10(numin.to(u.GHz).value),
                        np.log10(numax.to(u.GHz).value),
                        500)

        if params is None:
            params = self.params

        if dust:
            y = inufit_dust(*params)(x)
        elif dustT:
            y = inufit_dustT(*params)(x)
        else:
            y = inufit(*params)(x)

        y = u.Quantity(y, u.mJy)

        pl.loglog(x, y, **kwargs)
        pl.xlabel('Frequency (GHz)')
        pl.ylabel('Flux Density (mJy)')
        pl.title(plottitle)

        pl.errorbar(self.nu.value, self.flux.to(u.mJy).value,
                    yerr=self.fluxerr.to(u.mJy).value,
                    zorder=10,
                    fmt='s', markersize=6, alpha=0.75, **kwargs)

        self.physprops()
        if do_annotations:
            #pl.annotate("size (as): {0:0.2g}".format(self.srcsize), [annotation_xpos, .35],textcoords='axes fraction',xycoords='axes fraction')
            if hasattr(self, 'beta'):
                pl.annotate("$\\beta$: {0:0.3g}".format(self.beta), [annotation_xpos, .4],textcoords='axes fraction',xycoords='axes fraction',fontsize=20)
                pl.annotate("$T_{{dust}}$: {0:0.2g} K".format(self.dustT), [annotation_xpos, .35],textcoords='axes fraction',xycoords='axes fraction',fontsize=20)
            elif hasattr(self, 'alpha'):
                pl.annotate("$\\alpha$: {0:0.3g}".format(self.alpha), [annotation_xpos, .35],textcoords='axes fraction',xycoords='axes fraction',fontsize=20)

            pl.annotate("size (au): {0.value:0.2g}{0.unit:latex}".format(self.srcsize), [annotation_xpos, .3],textcoords='axes fraction',xycoords='axes fraction',fontsize=20)
            pl.annotate("mass (msun): {0.value:0.2g}{0.unit:latex}".format(self.mass), [annotation_xpos, .25],textcoords='axes fraction',xycoords='axes fraction',fontsize=20)
            pl.annotate("EM: {0.value:0.2g}{0.unit:latex}".format(self.em), [annotation_xpos, .2],textcoords='axes fraction',xycoords='axes fraction',fontsize=20)
            pl.annotate("Nu(Tau=1): {0:0.2g}".format(self.nutau), [annotation_xpos, .15],textcoords='axes fraction',xycoords='axes fraction',fontsize=20)
            pl.annotate("N(lyc): {0.value:0.2g}{0.unit:latex}".format(self.Nlyc), [annotation_xpos,.1],textcoords='axes fraction',xycoords='axes fraction',fontsize=20)
            pl.annotate("dens: {0.value:0.2g}{0.unit:latex}".format(self.dens), [annotation_xpos,.05],textcoords='axes fraction',xycoords='axes fraction',fontsize=20)

    def physprops(self):
        """
        Get the source size (au), density (cm^-3),
        mass (msun), and Nlyc of the UCHII

        Also return EM and nutau

        ERROR IN CURRENT VERSION
        """
        if self.resolved:
            self.srcsize = ((self.beamsize_as2 *
                             (self.dist_kpc.to(u.au))**2)**0.5).to(u.au,
                                                                   u.dimensionless_angles())
        else:

            flux_to_use = np.argmax(np.isfinite(self.flux))

            self.srcsize = np.sqrt(self.flux[flux_to_use] /
                                   (2*constants.k_B*self.Te)
                                   * (constants.c /
                                      (self.nu[flux_to_use]))**2
                                   * (self.dist_kpc)**2 / np.pi).to(u.au)
        self.dens = np.sqrt(self.em/(self.srcsize)).to(u.cm**-3)
        self.mass = (self.dens * 4.0/3.0 * np.pi * self.srcsize**3 * muh *
                     constants.m_p).to(u.M_sun)

        ## I think this is a solution to the Stromgren equation with... some....
        ## questionable values?
        #U = self.dens**(2/3.) * self.srcsize/u.pc
        ## From table 14.1 of Wilson+ 2009
        #U_O7 = 68 * u.pc * u.cm**(2/3.)

        #self.Nlyc = (8.04e46*(self.Te/u.K)**-.85 * (U/u.cm**-3)**3).decompose()

        # number rate of lyman continuum photons is just the recombination rate
        # assuming uniform density
        self.Nlyc = (4/3 * np.pi * self.dens**2 * alpha_b * self.srcsize**3).to(u.s**-1)
        


        return {'srcsize': self.srcsize,
                'density': self.dens,
                'mass': self.mass,
                'Nlyc': self.Nlyc,
                'EM': self.em,
                'nutau': self.nutau}



    # Cara test data:
    # nu = array([1.4,5,8.33]); flux=array([4.7,9.2,9.1]); err=array([.52,.24,.07])
    # em,nutau,normfac,chi2 = UCHIIfitter.emtau(nu,flux,err)

def dens(Qlyc=1e45*u.s**-1, R=0.1*u.pc, alpha_b=alpha_b):
    return (((3 * Qlyc)/(4 * np.pi * R**3 * alpha_b))**0.5).to(u.cm**-3)

def EM(Qlyc=1e45*u.s**-1, R=0.1*u.pc, alpha_b=alpha_b):
    return (R * (((3 * Qlyc)/(4 * np.pi * R**3 *
                              alpha_b))**0.5)**2).to(u.cm**-6*u.pc)

def Tb(Te=default_te, nu=95*u.GHz, EM=EM()):
    return Te * (1-np.exp(-tau(nu=nu, EM=EM, Te=Te)))
    #return (8.235e-2 * (Te/(u.K))**-0.35 * (nu/u.GHz)**-2.1 * (EM/u.cm**-6/u.pc)*u.K).to(u.K)

def Tb_beamdiluted(Te=default_te, nu=95*u.GHz, R=0.1*u.pc, Qlyc=1e45*u.s**-1,
                   beam=4000*u.au):
    tb = Tb(Te=Te, nu=nu, EM=EM(R=R, Qlyc=Qlyc))
    if beam < R:
        return tb
    else:
        return (tb * (R/beam)**2).to(u.K)

def Snu(Te=default_te, nu=95*u.GHz, R=0.1*u.pc, Qlyc=1e45*u.s**-1, beam=4000*u.au,
        angular_beam=0.5*u.arcsec):
    tb = Tb(Te=Te, nu=nu, EM=EM(R=R, Qlyc=Qlyc))
    if beam < R:
        return tb.to(u.mJy,
                     u.brightness_temperature(radio_beam.Beam(angular_beam),
                                              nu))
    else:
        return (tb * (R/beam)**2).to(u.mJy,
                                     u.brightness_temperature(radio_beam.Beam(angular_beam),
                                                              nu))

def snu_dust(density=1e4*u.cm**-3, Td=40*u.K, radius=4000*u.au,
             distance=8.4*u.kpc, cfreq=95*u.GHz):
    mass = (density * 2.8 * u.Da * 4/3. * radius**3).to(u.M_sun)
    #print(mass)
    beam = radio_beam.Beam((radius/distance).to(u.arcsec,
                                                u.dimensionless_angles()))
    flux = dust.snuofmass(nu=cfreq, mass=mass, beamomega=beam, temperature=Td,
                          distance=distance)
    return flux

def EM_of_T(TB, Te=default_te, nu=default_freq):
    " eqn 4.61 of Condon & Ransom inverted "
    return (-3.05e6 * (Te/(1e4*u.K))**1.35 * (nu/u.GHz)**2.1 * np.log(1-TB/Te)
            * u.cm**-6 * u.pc)

def qlyc_of_tb(TB, Te=default_te, nu=default_freq, radius=1*u.pc):
    EM = EM_of_T(TB, Te=Te, nu=nu)
    result = (4/3. * np.pi * radius**3 * alpha_b * EM / radius)
    return result.to(u.s**-1)
    return (-4/3. * np.pi * radius**3 * alpha_b * (3.28e-7)**-1 *
            (Te/(1e4*u.K))**1.35 * (nu/u.GHz)**2.1 * np.log(1-TB/Te) * u.cm**-6 *
            u.pc).to(u.s**-1)

__all__ = [tnu,Inu,unitfactor,freqfactor,inufit,emtau,mpfitfun,HIIregion]
