#!/usr/bin/env python

# compatibility with python 2/3:
from __future__ import print_function
from __future__ import division

import os, sys
import string
from pylab import *
import astropy.io.fits as pyfits
import numpy as np
import glob


from matplotlib import pylab as plt

import ac_settings as ac_set


def compute_flux(wave, flux, blaze, noise, ln_ctr, ln_win, bandtype, frac=True, test_plot=False):
    """
    Calculates the flux inside a bandpass. Interpolates between flux between pixel edges and bandpass limits.
    """
    step = 1e-4

    def interpolate_band_lims(array, wave, wmin, wmax, step):
        from scipy.interpolate import interp1d

        res_ratio = np.mean(np.diff(wave))/step

        mask = (wave >= wmin) & (wave <= wmax)
        wave_int_low = (wave[wave < wmin][-1], wave[wave >= wmin][0])
        wave_int_high = (wave[wave <= wmax][-1], wave[wave > wmax][0])

        interv_low = (wave >= wave_int_low[0]) & (wave <= wave_int_low[1])
        interv_high = (wave >= wave_int_high[0]) & (wave <= wave_int_high[1])

        wave_low = wave[interv_low]
        wave_high = wave[interv_high]

        array_low = array[interv_low]
        array_high = array[interv_high]

        interp_low = interp1d(wave_low, array_low, kind='linear')
        interp_high= interp1d(wave_high, array_high, kind='linear')

        wave_i_low = np.arange(np.min(wave_low), np.max(wave_low), step)
        array_i_low = interp_low(wave_i_low)

        wave_i_high = np.arange(np.min(wave_high), np.max(wave_high), step)
        array_i_high = interp_high(wave_i_high)

        wave_i = np.r_[wave_i_low, wave[mask], wave_i_high]

        array_i = np.r_[array_i_low/(array_i_low.size*res_ratio), array[mask], array_i_high/(array_i_high.size*res_ratio)]

        return wave_i, array_i

    ctr = ln_ctr
    win = ln_win

    print("Executing compute_flux")

    if blaze is None: blaze = np.ones(len(flux))

    # make all important values the size of px_size
    px_size = np.diff(wave)
    wave  = wave[1:]
    flux  = flux[1:]
    blaze = blaze[1:]

    # BANDPASS TYPES
    if bandtype == 'tri':
        wmin = ln_ctr - ln_win
        wmax = ln_ctr + ln_win
        # used for frac = false:
        cond = (wave > wmin) & (wave < wmax)
        bandfunc = -np.abs(wave-ln_ctr)/ln_win + 1.
        bandfunc = np.where(bandfunc > 0, bandfunc, bandfunc*0.0)
    if bandtype == 'sq':
        wmin = ln_ctr - ln_win/2.
        wmax = ln_ctr + ln_win/2.
        # used for frac = false:
        cond = (wave > wmin) & (wave < wmax)
        bandfunc = np.ones(len(wave))
        bandfunc = np.where(cond, 1., 0.)


    # HARPS DRS METHOD
    if frac == False:
        flux_i  = flux[cond]
        wave_i  = wave[cond]
        blaze_i = blaze[cond]
        dflux_i = flux[cond]/blaze[cond]
        bp_i    = bandfunc[cond]
        dflux_i = flux_i/blaze_i
    # ACTIN METHOD
    if frac == True:
        wave_i, flux_i = interpolate_band_lims(flux, wave, wmin, wmax, step)
        _, dflux_i = interpolate_band_lims(flux/blaze, wave, wmin, wmax, step)

        if bandtype == 'tri':
            bp_i = 1 - np.abs(wave_i - ctr)/win
            bp_i = np.where(bp_i > 0, bp_i, bp_i*0.0)
        elif bandtype == 'sq':
            bp_mask = (wave_i >= ctr - win/2) & (wave_i <= ctr + win/2)
            bp_i = np.where(bp_mask, 1, 0.0)


    r_neg_ln, _, _, flg_negflux = check_negflux(dflux_i, verb=False)

    # Flux sum and variance for line:
    if not noise: noise = 0.0
    f_sum = sum(dflux_i * bp_i)/win
    f_sum_var = sum((flux_i + noise**2) * bp_i**2)/win**2

    return f_sum, f_sum_var, bandfunc, flg_negflux, r_neg_ln


def check_negflux(flux, verb=True):
    flux = np.asarray(flux)

    # Total absolute flux in given spectral line:
    tot_flux = np.sum(abs(flux))

    # Number of pixels with negative flux
    neg_pixs = flux[flux < 0].size

    # Negative flux in given spectral line
    neg_flux = np.sum(flux[flux < 0])

    # Positive flux in given spectral line
    pos_flux = np.sum(flux[flux > 0])

    # Negative flux ratio for a given spectral line:
    r_neg_ln = abs(neg_flux)/tot_flux

    flg_negflux = "OK"

    if neg_pixs > 0:
        if verb:
            print(f"Negative flux detected")
        flg_negflux = "negFlux"

    return r_neg_ln, neg_flux, tot_flux, flg_negflux


def remove_output2(star_name, instr, file_type, save_output):

    file_rmv = "{}_{}_{}_data.rdb".format(star_name, instr, file_type)
    files = glob.glob(os.path.join(save_output, star_name, file_rmv))
    if files:
        for file in files:
            if os.path.isfile(file):
                print(f"Removing file {file}")
                os.remove(file)
    #else:
    #    print("There are no files to remove.")


def files_by_star_and_ftype(files):
    """
    Organize files by file type.
    """
    # different file types in files
    ft = []
    for k in range(len(files)):
        ft.append(get_file_type(files[k]))
    ft = list(set(ft))

    # different stars in files
    sp = []
    for k in range(len(files)):
        sp.append(os.path.split(files[k])[0])
    sp = list(set(sp))

    files_list = []
    for k in range(len(sp)):
        lists_sp = []
        for i in range(len(ft)):
            lists_ft = []
            for j in range(len(files)):
                if sp[k] in files[j] and ft[i] in files[j]:
                    lists_ft.append(files[j])
            if lists_ft:
                lists_sp.append(lists_ft)
        files_list.append(lists_sp)

    return files_list


def override_obj(obj, obj_name):
    """
    Override object name with name given in obj_name option.
    """
    if type(obj_name) is list and len(obj_name) == 1:
        return obj_name[0]
    elif type(obj_name) is list and len(obj_name) > 1:
        sys.exit("*** ERROR: obj_name requires only one name, more than one given.")
    else: return obj_name


def check_files(files):
    for k in range(len(files)):
        if not os.path.isfile(files[k]):
            sys.exit("*** ERROR: File {} not found.".format(files[k]))


def get_file_type(file):
    """
    Checks file name for known file types and returns it.
    """
    for k in range(len(ac_set.ftypes['all'])):
        if ac_set.ftypes['all'][k] in file:
            return ac_set.ftypes['all'][k]


def get_instr(fits_file):
    if ".rdb" not in fits_file:
        hdu = pyfits.open(fits_file)

        tel = hdu[0].header['TELESCOP']
        instr = hdu[0].header['INSTRUME']

        hdu.close()

        return tel, instr
    if ".rdb" in fits_file: return False, False


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
                print("*** ERROR: Cannot identify object.")
                return

    return obj


def check_targ(fits_file, targets):
    """
    Checks if a fits file belongs to a target in a list of targets.
    """
    print("Executing: check_targ")
    print("Targets = {}".format(targets))

    obj = get_target(fits_file)

    print("Object = {}".format(obj))
    if obj in targets: return True
    else:
        print("*** INFO: {} is not in the target list.".format(obj))
        print("*** INFO: file {}".format(fits_file))
        return False


def test_actin(test, path, calc_index):
    if not calc_index:
        calc_index = ("I_CaII", "I_Ha06")
    if test == "S1D":
        files = os.path.join(path, "test_files", "HD41248_1_1_S1D_A.fits")
    elif test == "S2D":
        files = os.path.join(path, "test_files", "HD41248_1_1_S2D_A.fits")
    elif test == "e2ds":
        files = os.path.join(path, "test_files", "HARPS.2003-12-13T06:19:48.371_e2ds_A.fits")
    elif test == "s1d":
        files = os.path.join(path, "test_files", "HARPS.2010-09-18T23:42:36.178_s1d_A.fits")
    elif test == "adp":
        files = os.path.join(path, "test_files", "ADP.2014-09-16T11:04:45.123.fits")
    elif test == 'rdb':
        files = os.path.join(path, "test_files", "2010-09-18T23:42:36.178_spec.rdb")
    else:
        print("*** ERROR:")
        print("*** Test can only be 'S1D', 'S2D', 'e2ds', 's1d', 'adp', or 'rdb'")
        return None, None
    return calc_index, files


def read_rdb(filename):
    """
    Reads an .rdb file and returns a dictionary with the headers as keys and data as lists ('output'), and also a list of headers ('keys').

    use: table = pyrdb.read_rdb(file)[0] for data
    use: table = pyrdb.read_rdb(file)[1] to get the keys by order
    """
    f = open(filename, 'r')
    data = f.readlines()
    f.close()

    key = str.split(data[0][:-1],'\t')
    output = {}
    for i in range(len(key)): output[key[i]] = []

    for line in data[2:]:
        qq = str.split(line[:-1],'\t')
        for i in range(len(key)):
            try: value = float(qq[i])
            except ValueError: value = qq[i]
            output[key[i]].append(value)

    return output, key


def save_rdb(dic, keys,file):
    """
    From a disctionary 'dic' saves the columns related to the specified 'keys' into an .rdb file named 'file'.
    """
    out = open(file,'w')
    n_keys = len(keys)

    for k in range(n_keys):
        if k != n_keys-1: out.write('{}\t'.format(keys[k]))
        else: out.write('{}'.format(keys[k]))
    out.write('\n')
    for k in range(n_keys):
        if k != n_keys-1: out.write('{}\t'.format('-'*len(keys[k])))
        else: out.write('{}'.format('-'*len(keys[k])))
    out.write('\n')

    for i in range(len(dic[keys[0]])):
        for k in range(n_keys):
            if k != n_keys-1: out.write('{}\t'.format(dic[keys[k]][i]))
            else: out.write('{}'.format(dic[keys[k]][i]))
        out.write('\n')
    out.close()


def add_rdb(dic,keys, file_name):
    """
    Adds data to an existing .rdb file. The 'keys' must match the headers already present in the file.
    """
    out = open(file_name,'a')
    n_keys = len(keys)
    for i in range(len(dic[keys[0]])):
        for k in range(n_keys):
            if k != n_keys-1: out.write('{}\t'.format(dic[keys[k]][i]))
            else: out.write('{}\t'.format(dic[keys[k]][i]))
        out.write('\n')
    out.close()


def plot_params(width=6, height=3.5):
        """
        Parameters for plots.
        """
        rcdefaults()
        rc('axes', linewidth=1)
        rc('lines', markeredgewidth=0.5)

        rc('errorbar', capsize=2)

        rc('ytick', direction='in')
        rc('xtick', direction='in')

        rc('ytick', right='True')
        rc('xtick', top='True')

        rc('xtick.minor', visible='True')


        rc('xtick.major', size=7)
        rc('ytick.major', size=7)
        rc('xtick.minor', size=3)
        rc('ytick.minor', size=3)
        rc("font", family="sans-serif")
        rc("font", size=10)

        rc("figure.subplot", left=(0.15))
        rc("figure.subplot", right=(0.93))
        rc("figure.subplot", bottom=(0.12))
        rc("figure.subplot", top=(0.95))

        return width, height
