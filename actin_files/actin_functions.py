#!/usr/bin/env python

from pylab import *
import astropy.io.fits as pyfits
import numpy as np


# call this function get_win and redo calc_ind, etc.
def frac_pixels(wave, flux, wmin, wmax):
    """
    TESTED
    """

    pix_size = np.diff(wave)
    # reduce wave of first element to have same dim as pix_size
    wave = wave[1:]
    flux = flux[1:]

    cond = (wave >= wmin) & (wave <= wmax)

    frac_left = 1 - (wave[cond][0] - wmin)/pix_size[(wave < wmin)][-1]
    frac_right = (wmax - wave[cond][-1])/pix_size[(wave < wmax)][-1]

    # Diaz method:
    # fraction_left = 0.5 - (wmin - wave[cond])/pix_size[cond]
    # fraction_right = 0.5 - (wave[cond] - wmax)/pix_size[cond]
    # fraction_left = np.where(fraction_left > 1, 1.0, fraction_left)
    # fraction_right = np.where(fraction_right > 1, 1.0, fraction_right)
    # frac_left = fraction_left[0]
    # frac_right = fraction_right[-1]
    # frac_flux_left = flux[cond][0] * frac_left
    # frac_flux_right = flux[cond][-1] * frac_right
    # flux_win = flux[cond]
    # flux_win[0] = frac_flux_left
    # flux_win[-1] = frac_flux_right

    pixels_win = len(wave[cond]) + frac_left + frac_right

    frac_flux_left = flux[(wave < wmin)][-1] * frac_left
    frac_flux_right = flux[(wave > wmax)][0] * frac_right

    flux_win = flux[cond]
    flux_win = np.insert(flux_win, 0, frac_flux_left)
    flux_win = np.append(flux_win, frac_flux_right)

    #frac_flux = frac_flux_left + frac_flux_right

    return flux_win#, frac_flux


def compute_flux(wave, flux, blaze, ln_ctr, ln_win, ln_c, bandtype, weight=None):

    if weight == 'blaze':
        weight = blaze

    # no weight, normalized by number of pixels
    elif weight == None:
        weight = np.ones(len(flux))

    else:
        print "*** ERROR: compute_flux weight option must be 'blaze' or None."
        return

    wmin = ln_ctr - ln_win/2.
    wmax = ln_ctr + ln_win/2.

    flux_win = frac_pixels(wave, flux, wmin, wmax)
    blaze_win = frac_pixels(wave, blaze, wmin, wmax)

    print "MEAN FLUX_WIN = %s" % np.mean(flux_win)
    print "MEAN BLAZE_WIN = %s" % np.mean(blaze_win)
    print "SUM FLUX/BLAZE WIN = %s" % sum(flux_win/blaze_win)
    print "N FLUX_WIN = %s" % len(flux_win)

    flux_win_norm = flux_win/blaze_win

    if bandtype == 'sq':
        bandfunc = None

    if bandtype == 'tri':
        # Add triangular weight function for line:
        bandfunc = np.where(wave >= ln_ctr, -(wave - ln_ctr)/ln_win*2. + 1, (wave - ln_ctr)/ln_win*2. + 1)

        bandfunc = np.where(bandfunc > 0, bandfunc, bandfunc*0.0)

        weight = weight * bandfunc

    weight_win = frac_pixels(wave, weight, wmin, wmax)

    f_sum = ln_c * sum(flux_win_norm * weight_win)/sum(weight_win)


    # err = sqrt(flux), var = err**2
    variance = flux_win/blaze_win**2

    f_sum_var = sum(variance * weight_win**2)/sum(weight_win)**2
    f_sum_err = ln_c * np.sqrt(f_sum_var)

    return f_sum, f_sum_err, bandfunc


def get_target(fits_file):

    hdu = pyfits.open(fits_file)

    try: obj = hdu[0].header['OBJECT']
    except:
        try: obj = hdu[0].header['ESO OBS TARG NAME']
        except:
            try: obj = hdu[0].header['TNG OBS TARG NAME']
            except:
                print "Cannot identify object"
                return

    return obj



# not used yet, not tested yet
def calc_ind(table):
    """
    table with keys: ind_var, flux, flux_err, flg, frac_neg, snr
    """

    var = np.asarray(table['ind_var'])
    flux = np.asarray(table['flux'])
    flux_err = np.asarray(table['flux_err'])

    num = flux['L' in var]
    num_err = flux_err['L' in var]
    denom = flux['R' in var]
    denom_err = flux_err['R' in var]

    ind = sum(num) / sum(denom)

    # error using propagation of errors for lines and ref lines
    ind_err = np.sqrt(sum(num_err**2) + ind**2 * sum(denom_err**2))/sum(denom)

    return ind, ind_err





def flag_negflux(flux):
    """
    Tests if flux as negative values and returns flag 'flg' as 'negFlux' if they are present, None otherwise, and the fraction of pixels with negative values, 'frac_n_neg'.
    """

    flag_array = np.where(flux < 0.0, 'negFlux', None)

    n_neg = len(flag_array[flag_array == 'negFlux'])
    frac_n_neg = float(n_neg)/len(flux)

    if 'negFlux' in flag_array: flg = 'negFlux'
    else: flg = None

    return flg, frac_n_neg


def plot_params(width=6, height=3.5):
    """
    Parameters for plots.
    """

    rcdefaults()
    rc('axes',linewidth = 1.2)
    rc('lines',markeredgewidth=1)
    rc('xtick.major',size=7)
    rc('ytick.major',size=7)
    rc('xtick.minor',size=4)
    rc('ytick.minor',size=4)
    rc("font", family="sans-serif")
    rc("font", size=10)

    rc("figure.subplot", left=(0.15))
    rc("figure.subplot", right=(0.93))
    rc("figure.subplot", bottom=(0.12))
    rc("figure.subplot", top=(0.95))

    return width, height


def check_targ(fits_file, targets):
    """
    Checks if a fits file belongs to a target in a list of targets.
    """

    print "Targets = %s" % targets

    print fits_file

    try: fits = pyfits.open(fits_file)
    except:
        print "*** ERROR: Cannot read %s." % fits_file
        return False

    # to chose if using HARPS or HARPS-N
    tel = fits[0].header['TELESCOP']

    try: obj = fits[0].header['OBJECT']
    except:
        try: obj = fits[0].header['%s OBS TARG NAME' % tel[:3]]
        except:
            print "*** ERROR: Cannot get object identification."
            return False

    print "Object = %s" % obj

    if obj in targets: return True
    else: return False
