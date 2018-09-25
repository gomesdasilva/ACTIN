#!/usr/bin/env python

# compatibility with python 2/3:
from __future__ import print_function
from __future__ import division

import os
import string
from pylab import *
import astropy.io.fits as pyfits
import numpy as np


def frac_pixels(wave, flux, wmin, wmax):
    """
    For given 'wave' and 'flux' (wavelength and flux) calculates the fraction of pixels inside a bandpass with limits 'wmin' and 'wmax' and returns the total flux inside the bandpass taking the fraction of pixels into account.
    """

    pix_size = np.diff(wave)
    # reduce wave of first element to have same dim as pix_size
    wave = wave[1:]
    flux = flux[1:]

    cond = (wave >= wmin) & (wave <= wmax)

    frac_left = 1 - (wave[cond][0] - wmin)/pix_size[(wave < wmin)][-1]
    frac_right = (wmax - wave[cond][-1])/pix_size[(wave < wmax)][-1]

    pixels_win = len(wave[cond]) + frac_left + frac_right

    frac_flux_left = flux[(wave < wmin)][-1] * frac_left
    frac_flux_right = flux[(wave > wmax)][0] * frac_right

    flux_win = flux[cond]
    flux_win = np.insert(flux_win, 0, frac_flux_left)
    flux_win = np.append(flux_win, frac_flux_right)

    return flux_win, pixels_win


def compute_flux(wave, flux, blaze, ln_ctr, ln_win, ln_c, bandtype, weight=None, norm='npixels'):
    """
    Calculates the sum of the flux inside a bandpass, using different bandpass functions and weights and taking into account fractions of pixels in the bandpass.

    Parameters:
    -----------
    wave, flux, blaze : lists
        Wavelength [angstroms], flux and blaze function.
    ln_ctr, ln_win : float
        Centre and bandpass of the spectral line to compute flux [angstroms].
    ln_c : float
        Constant multiplied to the integrated flux.
    bandtype : string
        Function to be used in the integration of flux. If 'sq', uses a square function with limits 'ln_win', if 'tri' a triangular function with full-width-at-half-maximum given by 'ln_win'.
    weight : string
        Function to weight the integrated flux. If 'blaze' the flux is multiplied by the blaze function, if None (default) the flux is not weighted (default).
    norm : string
        Normalization of the flux: if 'band' the sum is normalized by the bandpass wavelength value in angstroms, if 'npixels' by the number of pixels in the bandpass (default), if 'weight' by the sum of the weight function inside the bandpass, if None the integrated flux is not normalized.

    Returns:
    --------
    f_sum, f_sum_err : float
        Integrated flux and respective error (considering only photon noise)
    bandpass : np.float
        Bandpass function. If None is a square function with limits ln_ctr - ln_win/2 and ln_ctr + ln_win/2.
    """

    if weight == 'blaze': weight = blaze
    elif weight == None: weight = np.ones(len(flux))
    else:
        print("*** ERROR: 'weight' option must be 'blaze' or None, '%s' was given" % weight)
        quit()

    wmin = ln_ctr - ln_win/2.
    wmax = ln_ctr + ln_win/2.

    flux_win = frac_pixels(wave, flux, wmin, wmax)[0]
    blaze_win = frac_pixels(wave, blaze, wmin, wmax)[0]

    flux_win_norm = flux_win/blaze_win

    if bandtype == 'sq':
        bandfunc = None

    elif bandtype == 'tri':
        # Add triangular weight function for line:
        bandfunc = np.where(wave >= ln_ctr, -(wave - ln_ctr)/ln_win*2. + 1, (wave - ln_ctr)/ln_win*2. + 1)

        bandfunc = np.where(bandfunc > 0, bandfunc, bandfunc*0.0)

        weight = weight * bandfunc
    else:
        print("*** ERROR: 'bandtype' (in config file) must be either 'sq' or 'tri' but '%s' was given. The config file path can be found by calling 'actin -cfg True'." % bandtype)
        quit()

    weight_win = frac_pixels(wave, weight, wmin, wmax)[0]

    # err = sqrt(flux), var = err**2
    variance = ln_c * flux_win/blaze_win**2

    if norm == 'band':
        f_sum = ln_c * sum(flux_win_norm * weight_win)/ln_win
        f_sum_var = sum(variance * weight_win**2)/ln_win**2
    elif norm == 'npixels':
        f_sum = ln_c * sum(flux_win_norm * weight_win)/len(weight_win)
        f_sum_var = sum(variance * weight_win**2)/len(weight_win)**2
    elif norm == 'weight':
        f_sum = ln_c * sum(flux_win_norm * weight_win)/sum(weight_win)
        f_sum_var = sum(variance * weight_win**2)/sum(weight_win)**2
    elif norm == None:
        f_sum = ln_c * sum(flux_win_norm * weight_win)
        f_sum_var = sum(variance * weight_win**2)
    else:
        print("*** ERROR: 'norm' option must be either 'band', 'npixels', 'weight', or None, but '%s' was given" % norm)
        quit()

    f_sum_err = np.sqrt(f_sum_var)

    return f_sum, f_sum_err, bandfunc


def flag_negflux(flux):
    """
    Tests if flux has negative values and returns flag 'flg' as 'negFlux' if detected, None otherwise, and the fraction of pixels with negative values of flux, 'frac_neg'.
    """
    negflux_array = np.where(flux < 0.0, flux, 0.0)
    #negflux = sum(negflux_array)

    negflux_only = [negflux_array[x] for x in range(len(negflux_array)) if negflux_array[x] < 0.0]

    # fraction of pixels with negative flux
    frac_neg = len(negflux_only)/len(flux)

    #frac_neg = sum(abs(negflux))/sum(flux)

    flag_array = np.where(flux < 0.0, 'negFlux', None)

    if 'negFlux' in flag_array: flg = 'negFlux'
    else: flg = None

    return flg, frac_neg


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


def remove_output(files, save_output, targ_list=None):
    """
    Removes output directories for given fits files 'files', 'save_output' directory and list of targets 'targ_list' (if available).
    """
    err_msg = None
    try:
        if targ_list is None:
            obj = [get_target(files[k]) for k in range(len(files))]
            obj = list(set(obj)) # remove duplicates
        elif targ_list is not None:
            obj = targ_list

        if 'ADP' in files[0]: ft = 'ADP'
        if 'e2ds' in files[0]: ft = 'e2ds'
        if 's1d' in files[0]: ft = 's1d'
        if 'rdb' in files[0]: ft = 'rdb'

        for k in range(len(obj)):
            ff = "%s/%s/%s_%s_actin.rdb" % (save_output, obj[k], obj[k], ft)
            os.remove(ff)
            print("Output file removed: %s" % ff)

    except: err_msg = "*** ERROR: Could not delete output file(s)"
    return err_msg


def get_target(fits_file):
    """
    Returns the object targeted in the fits file 'fits_file'.
    """

    hdu = pyfits.open(fits_file)

    try:
        obj = hdu[0].header['OBJECT']
        hdu.close()
    except:
        try:
            obj = hdu[0].header['ESO OBS TARG NAME']
            hdu.close()
        except:
            try:
                obj = hdu[0].header['TNG OBS TARG NAME']
                hdu.close()
            except:
                print("Cannot identify object")
                return

    return obj


def check_targ(fits_file, targets):
    """
    Checks if a fits file belongs to a target in a list of targets.
    """

    print("Targets = %s" % targets)

    print(fits_file)

    try: fits = pyfits.open(fits_file)
    except:
        print("*** ERROR: Cannot read %s." % fits_file)
        return False

    # to chose if using HARPS or HARPS-N
    tel = fits[0].header['TELESCOP']

    try: obj = fits[0].header['OBJECT']
    except:
        try: obj = fits[0].header['%s OBS TARG NAME' % tel[:3]]
        except:
            print("*** ERROR: Cannot get object identification.")
            return False

    print("Object = %s" % obj)

    if obj in targets: return True
    else: return False


def read_rdb(filename):
    """
    Reads an .rdb file and returns a dictionary with the headers as keys and data as lists ('output'), and also a list of headers ('keys').
    """
    # use: table = pyrdb.read_rdb(file)[0] for data
    # use: table = pyrdb.read_rdb(file)[1] to get the keys by order

    f = open(filename, 'r')
    data = f.readlines()
    f.close()

    key = string.split(data[0][:-1],'\t')
    output = {}
    for i in range(len(key)): output[key[i]] = []

    for line in data[2:]:
        qq = string.split(line[:-1],'\t')
        for i in range(len(key)):
            try: value = float(qq[i])
            except ValueError: value = qq[i]
            output[key[i]].append(value)

    return output, key


def save_rdb(dic,keys,file):
    """
    From a disctionary 'dic' saves the columns related to the specified 'keys' into an .rdb file named 'file'.
    """

    out = open(file,'w')
    n_keys = len(keys)

    for k in range(n_keys):
        if k != n_keys-1: out.write('%s\t' % (keys[k]))
        else: out.write('%s' % (keys[k]))
    out.write('\n')
    for k in range(n_keys):
        if k != n_keys-1: out.write('%s\t' % ('-'*len(keys[k])))
        else: out.write('%s' % ('-'*len(keys[k])))
    out.write('\n')

    for i in range(len(dic[keys[0]])):
        for k in range(n_keys):
            if k != n_keys-1: out.write('%s\t' % (dic[keys[k]][i]))
            else: out.write('%s' % (dic[keys[k]][i]))
        out.write('\n')
    out.close()


def add_rdb(dic,keys,file_name):
    """
    Adds data to an existing .rdb file. The 'keys' must match the headers already present in the file.
    """

    out = open(file_name,'a')
    n_keys = len(keys)
    for i in range(len(dic[keys[0]])):
        for k in range(n_keys):
            if k != n_keys-1: out.write('%s\t' % (dic[keys[k]][i]))
            else: out.write('%s\t' % (dic[keys[k]][i]))
        out.write('\n')
    out.close()
