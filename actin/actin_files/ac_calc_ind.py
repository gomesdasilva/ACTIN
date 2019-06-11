#!/usr/bin/env python

# compatibility with python 2/3:
from __future__ import print_function
from __future__ import division

import os, sys ### test
import numpy as np

# ACTIN FILES
import ac_settings as ac_set
import ac_get_win


def check_lines(wave, sel_lines):
    """
    Tests if the selected lines from config file fit inside the spectrograph
    wavelength range and fit inside any spectral orders for the case of 2d
    spectrum.
    """

    print("\nCHECKING LINES FOR WAVELENGTH RANGE AND SP. ORDERS")
    print("--------------------------------------------------")

    if type(wave[0]) is np.ndarray: # 2d spec
        min_spec_wave = wave[0][0]
        max_spec_wave = wave[-1][-1]
        spec_type = '2d'
    if type(wave[0]) in [np.float, np.float64]: # 1d spec
        min_spec_wave = wave[0]
        max_spec_wave = wave[-1]
        spec_type = '1d'

    # For each row (sp. line) in the config table calculate the min and max values of bandwidth
    rows = len(sel_lines['ln_id'])
    for k in range(rows):
        ln_id = sel_lines['ln_id'][k]
        ln_ctr = sel_lines['ln_ctr'][k]
        ln_win = sel_lines['ln_win'][k]

        if ln_win <= 0:
            sys.exit("*** ERROR: line {} bandwidth is not positive.".format(ln_id))

        min_wave = ln_ctr - ln_win/2.
        max_wave = ln_ctr + ln_win/2.

        # Check if line fits inside spectral range
        if min_wave < min_spec_wave or max_wave > max_spec_wave:
            sys.exit("*** ERROR: Line {} bandwidth outside spectral range.".format(ln_id))
        else: print("Line {} inside spectral range".format(ln_id))

        # If wave is 2d check if line fits inside sp. order
        if spec_type == '2d':
            order = None
            ln_ctr_orders = []
            order = []
            for i in range(len(wave)):
                if min_wave > wave[i][0] and max_wave < wave[i][-1]:
                    order.append(i)
                # used to show potencial orders with wavelength range
                if ln_ctr > wave[i][0] and ln_ctr < wave[i][-1]:
                    ln_ctr_orders.append(i)

            if order is None:
                print("*** ERROR: Could not determine sp. order for {}".format(ln_id))
                print("\tmin_wave = {:.2f}".format(min_wave))
                print("\tmax_wave = {:.2f}".format(max_wave))
                print("\tThe closest orders are:")
                for k in range(len(ln_ctr_orders)):
                    closest_ord = ln_ctr_orders[k]
                    print("\tOrder {}: {:.2f}-{:.2f}".format(closest_ord, wave[closest_ord][0], wave[closest_ord][-1]))
                sys.exit()
            else:
                for i in range(len(order)):
                    print("Line {} inside spectral order {}".format(ln_id, order[i]))
    return


def calc_flux_lines(data, sel_lines, ln_plts=False, frac=True):
    """
    Calculates the sum of the flux and associated errors for all spectral lines required to calculate the selected indices.

    Parameters:
    -----------
    data : dict
        Dictionary with data returned from fits files.

        The required keys are:

        ==========  ========================================================
        keys        Description
        ----------  --------------------------------------------------------
        flux        list of lists : Flux per pixel per order.
        wave        list of lists : Wavelength calibrated for BERV and RV
                    (at rest frame) per pixel per order [angstroms].
        blaze       list of lists : Blaze function.
        snr         list : SNR at each spectral order.
        obj         str : Object identification.
        date        str : Date of observation in the fits file format.
        ==========  ========================================================

    sel_lines : dict
        Dictionary containing the identification of the indices selected and
        the parameters of the spectral lines required for each index.

        The required keys are:

        ==========  ========================================================
        keys        Description
        ----------  --------------------------------------------------------
        ind_id      str : Index identification.
        ind_var     str : Variable assigned to a given line to be used in
                    the index equation. Ex: 'L1', 'L2', etc, for the core
                    lines, and 'R1', 'R2', etc, for reference lines.
        ln_id       str : Spectral line identification.
        ln_ctr      float : Wavelength of the line centre [angstroms].
        ln_win      float : Bandpass around the line centre to be used in
                    the flux integration [angstroms].
        bandtype    str : Type of bandpass used to integrate flux.
        ==========  ========================================================

    ln_plts : {str, False} (optional)
        Path for the saved line plots. If False, the plots are not saved
        (default).
    frac : bool (optional)
        Use fractional pixels if 'True' (default), use integral pixels if 'False'.

    Returns:
    --------
    sel_lines : dict
        Dictionary containing the identification of the indices selected and
        the parameters of the spectral lines required for each index.

        Included in the output are also the keys used as input (see above).

        The returned keys are:

        ==========  ========================================================
        keys        Description
        ----------  --------------------------------------------------------
        flux        list : Integrated flux for each line.
        error       list : Errors on integrated flux for each line.
        snr         list : SNR for each line.
        flg         list : Flag for each line, 'negFlux' if negative flux
                    found inside any bandpass.
        frac_neg    list : Fraction of pixels with negative values of flux
                    for each line.
        npixels     float : Number of pixels inside the bandpass.
        ==========  ========================================================
    """

    print()
    print("CALCULATING FLUX IN SELECTED LINES")
    print("----------------------------------")

    wave = data['wave']
    flux = data['flux']
    blaze = data['blaze']

    obj = data['obj']
    obs_date = data['obs_date']
    snr = data['snr']

    sel_lines['flux'] = []
    sel_lines['error'] = []
    sel_lines['snr'] = []
    sel_lines['flg'] = []
    sel_lines['frac_neg'] = []
    sel_lines['npixels'] = []

    ln_ctr = sel_lines['ln_ctr']
    ln_win = sel_lines['ln_win']
    ln_id = sel_lines['ln_id']
    ind_var = sel_lines['ind_var']
    bandtype = sel_lines['bandtype']

    if 'error_pixel' in list(data):
        # case where pixel errors are given in rdb file as "error_pixel"
        print("Using pixel errors from input rdb file.")
        err = data['error_pixel']
    else: err = None

    for k in range(len(ln_id)):
        print()
        print("Computing flux in line {}".format(ln_id[k]))
        print("-----------------------{}".format('-'*len(ln_id[k])))

        print("Using {} bandpass".format(bandtype[k]))

        win = ac_get_win.get_win(data, ln_id[k], ln_ctr[k], ln_win[k], bandtype[k], blaze=blaze, err=err, frac=frac, ln_plts=ln_plts)

        sel_lines['flux'].append(win['sum'])
        sel_lines['error'].append(win['sum_err'])
        sel_lines['snr'].append(win['snr'])
        sel_lines['flg'].append(win['flg'])
        sel_lines['frac_neg'].append(win['frac_neg'])
        sel_lines['npixels'].append(win['npixels'])

    return sel_lines


def calc_ind(sel_lines):
    """
    Calculates the indices identified in sel_lines as 'ind_id'.

    sel_lines : dict
        Dictionary containing the identification of the indices selected and
        the parameters of the spectral lines required for each index.

        The required keys are:

        ==========  ========================================================
        keys        Description
        ----------  --------------------------------------------------------
        ind_id      str : Index identification.
        ind_var     str : Variable assigned to a given line to be used in
                    the index equation. Ex: 'L1', 'L2', etc, for the core
                    lines, and 'R1', 'R2', etc, for reference lines.
        ln_id       str : Spectral line identification.
        flux        list : Integrated flux for each line.
        error       list : Errors on integrated flux for each line.
        snr         list : SNR for each line.
        flg         list : Flag for each line, 'negFlux' if negative flux
                    found inside any bandpass.
        frac_neg    list : Fraction of pixels with negative values of flux
                    for each line.
        npixels     float : Number of pixels inside the bandpass.
        ==========  ========================================================

    Returns:
    --------
    index : dict
        Dictionary containing the index values, errors and related
        information.

        The returned keys are:

        ==========  ========================================================
        keys        Description
        ----------  --------------------------------------------------------
        index       str : Index identification.
        value       float : Index value.
        error       float : Index error.
        flg         {str, None} : Index flag, 'negFlux' if negative flux
                    found inside any bandpass.
        mfrac_neg   float : Maximum fraction of pixels with negative values
                    of flux taking into account all lines used to compute
                    index.
        snr         float : Mean SNR of all lines used to compute index.
        npixels     float : Number of pixels inside the bandpass.
        ==========  ========================================================
    """

    print()
    print("CALCULATING INDICES")
    print("-------------------")

    # remove duplicates of ind_id and gives a list of selected indices
    sel_ind = list(set(sel_lines['ind_id']))
    sel_ind = np.asarray(sel_ind)

    index = {}
    index['index'] = []
    index['value'] = []
    index['error'] = []
    index['flg'] = []
    index['mfrac_neg'] = []
    index['snr'] = []
    npixels_list = []

    print("index\tvalue\terror\t\tsnr\tflag\tmfrac_neg")
    print("-----\t-----\t-----\t\t---\t----\t---------")

    ind_ids = np.asarray(sel_lines['ind_id'])
    rows = len(sel_lines['ln_id'])
    for i in range(len(sel_ind)): # each index

        var = [sel_lines['ind_var'][k] for k in range(rows) \
                                                if ind_ids[k] == sel_ind[i]]
        flux = [sel_lines['flux'][k] for k in range(rows) \
                                                if ind_ids[k] == sel_ind[i]]
        err = [sel_lines['error'][k] for k in range(rows) \
                                                if ind_ids[k] == sel_ind[i]]
        flg = [sel_lines['flg'][k] for k in range(rows) \
                                                if ind_ids[k] == sel_ind[i]]
        frac_neg = [sel_lines['frac_neg'][k] for k in range(rows) \
                                                if ind_ids[k] == sel_ind[i]]
        snr = [sel_lines['snr'][k] for k in range(rows) \
                                                if ind_ids[k] == sel_ind[i]]
        ln_c = [sel_lines['ln_c'][k] for k in range(rows) \
                                                if ind_ids[k] == sel_ind[i]] ##
        npixels = [sel_lines['npixels'][k] for k in range(rows) \
                                                if ind_ids[k] == sel_ind[i]]

        # Maximum fraction of flux with negative values of all lines
        mfrac_neg = max(frac_neg)

        if "negFlux" in flg: flg_ind = 'negFlux'
        else: flg_ind = None

        # Mean snr of lines:
        try: snr_ind = np.mean([float(snr[k]) for k in range(len(snr))])
        except: snr_ind = None

        for k in range(len(var)):
            if 'L' not in var[k] and 'R' not in var[k]:
                msg="*** ERROR: 'ind_var' variable (in config file config_lines.txt) must start with either an 'L' for core line or 'R' for reference line. Value given was '{}'".format(var[k])
                sys.exit(msg)

        # Add line variables for numerator or denominator:
        num = [ln_c[k]*flux[k] for k in range(len(var)) if 'L' in var[k]]
        num_err = [ln_c[k]*err[k] for k in range(len(var)) if 'L' in var[k]]
        denom = [ln_c[k]*flux[k] for k in range(len(var)) if 'R' in var[k]]
        denom_err = [ln_c[k]*err[k] for k in range(len(var)) if 'R' in var[k]]

        num = np.asarray(num)
        denom = np.asarray(denom)
        num_err = np.asarray(num_err)
        denom_err = np.asarray(denom_err)

        ind = sum(num) / sum(denom)

        # Error using propagation of errors for lines and ref lines
        ind_err = np.sqrt(sum(num_err**2) + ind**2 * sum(denom_err**2)) /sum(denom)

        if snr_ind is not None: snr_ind = round(snr_ind, 2)

        index['index'].append(sel_ind[i])
        index['value'].append(ind)
        index['error'].append(ind_err)
        index['flg'].append(flg_ind)
        index['mfrac_neg'].append(mfrac_neg)
        index['snr'].append(snr_ind)
        npixels_list.append(npixels)

        print("{}\t{:.4f}\t{:.6f}\t{}\t{}\t{:.4f}".format(index['index'][i], index['value'][i], index['error'][i], index['snr'][i], index['flg'][i], index['mfrac_neg'][i]))

    # npixels data
    for k in range(len(sel_lines['ln_id'])):
        index['{}_npixels'.format(sel_lines['ln_id'][k])] = sel_lines['npixels'][k]

    return index
