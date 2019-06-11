#!/usr/bin/env python

# compatibility with python 2/3:
from __future__ import print_function
from __future__ import division

import sys, os
import numpy as np
from matplotlib import pylab as plt

import ac_tools
import ac_settings as ac_set


def sel_order(wave_2d, ln_ctr, ln_win, file_type):
    """
    Selects a spectral order for a line based and its bandpass.

    Parameters:
    -----------
    wave_2d : list of lists
        Wavelength in each spectral order [angstroms].
    ln_ctr : float
        Centre of the spectral line [angstroms].
    ln_win : float
        Bandwidth around the line where the flux is to be integrated
        [angstroms].

    Returns:
    --------
    order : {int, None}
        Spectral order for the given line parameters, None if no order is
        found.
    """
    min_wave = ln_ctr-ln_win/2.
    max_wave = ln_ctr+ln_win/2.

    ord = None
    orders = []
    for k in range(len(wave_2d)):
        if min_wave > wave_2d[k][0] and max_wave < wave_2d[k][-1]:
            ord = k
            orders.append(ord)

    print("Orders available =", orders)

    if file_type == "S2D":
        ord_type = "espresso"
    else:
        ord_type = "harps" #####

    print("ord_type =", ord_type)

    if ord_type == "harps":
        # Using two orders for the CaIIK line (as in Harps pipeline)
        if ln_ctr != 3933.664:
            orders = [orders[-1]]

    if ord_type == "espresso":
        # first two for CaIIK
        if ln_ctr == 3933.664:
            rders = orders[:-2]
        # last two for CaIIH
        if ln_ctr == 3968.47:
            orders = orders[2:]

    if ord_type == "last": # First actin version
        orders = [orders[-1]]

    if ord_type == "all": # All possible orders
        orders = orders

    print("Selected order =", orders)

    if ord is None:
        print("\n*** ERROR: Could not determine spectral order for:")
        print("*** min_wave = {:.2f}".format(min_wave))
        print("*** max_wave = {:.2f}".format(max_wave))

    return orders, ord_type


def line_plot(data, ord, ln_id, ln_ctr, ln_win, bandfunc, bandtype, out_dir):
    """
    Saves spectral lines plots.

    Parameters:
    -----------
    wave, flux : lists
        Wavelength [angstroms] and flux.
    obj : str
        Object identification.
    date : str
        Date of observations in fits file format.
    ln_id : str
        Spectral line identification.
    ln_ctr, ln_win : float
        Spectral line centre and bandpass [angstroms].
    bandfunc : list
        Bandpass function.
    file_type : str
        Identification of the type of file used. Options are 'S2D', 'S1D', 'e2ds', 's1d', 'ADP', and 'rdb'.
    out_dir : str
        Output directory to save the plots.

    Returns:
    --------
    Plots of the regions around the spectral lines used to calculate the
    indices. The plots are saved in the directory: 'out_dir'/'obj'.
    """

    print("Executing line_plot")

    obj       = data['obj']
    obs_date  = data['obs_date']
    file_type = data['file_type']
    instr     = data['instr']

    # 2d spectrum
    if isinstance(ord, int):
        wave = data['wave'][ord]
        flux = data['flux'][ord]/data['blaze'][ord]
        snr  = data['snr'][ord]
    # 1d spectrum
    elif ord is None:
        wave       = data['wave']
        flux       = data['flux']
        median_snr = data['median_snr']
    else:
        sys.exit("***ERROR: 'ord' must be either 'None' or instance of 'int'.")

    if bandtype == 'sq':
        wmin = ln_ctr - ln_win/2
        wmax = ln_ctr + ln_win/2
    if bandtype == 'tri':
        wmin = ln_ctr - ln_win
        wmax = ln_ctr + ln_win

    # Computing data with fractional pixels
    wave = wave[1:]
    flux = flux[1:]

    cond_reg = (wave > wmin - 4) & (wave < wmax + 4)
    flux_reg = flux[cond_reg]
    wave_reg = wave[cond_reg]

    cond = (wave > wmin) & (wave < wmax)

    dwave_l = wave[cond][0]-wmin
    dwave_r = wmax-wave[cond][-1]

    wave_win = np.insert(wave[cond], 0, wave[cond][0]-dwave_l)
    wave_win = np.append(wave_win, wave[cond][-1]+dwave_r)

    cond_l = (wave < wave[cond][0])
    cond_r = (wave > wave[cond][-1])

    flux_win = np.insert(flux[cond], 0, flux[cond_l][-1])
    flux_win = np.append(flux_win, flux[cond_r][0])

    flux_norm = (flux_reg-min(flux_reg))/(max(flux_reg)-min(flux_reg))

    wave_tot = np.r_[wave[cond_l], wave_win, wave[cond_r]]
    flux_tot = np.r_[flux[cond_l], flux_win, flux[cond_r]]

    if bandtype == 'sq':
        bandfunc_win = np.insert(bandfunc[cond], 0, 1.)
        bandfunc_win = np.append(bandfunc_win, 1.)
    if bandtype == 'tri':
        bandfunc_win = np.insert(bandfunc[cond], 0, 0.)
        bandfunc_win = np.append(bandfunc_win, 0.)

    band_l = np.zeros(len(wave[cond_l]))
    band_r = np.zeros(len(wave[cond_r]))

    width, height = ac_tools.plot_params()
    plt.figure(figsize=(width, height))

    plt.plot(wave_reg, flux_norm, color='k', ls='-', marker='', linewidth=1)
    #plt.axhline(0.0, c='k',ls=':')

    if bandtype == 'sq':
        bandfunc_tot = np.r_[band_l, 0., bandfunc_win[:-1], band_r]
        plt.plot(wave_tot, bandfunc_tot, color='r', ls='-', marker='', drawstyle='steps-pre')
    if bandtype == 'tri':
        bandfunc_tot = np.r_[band_l, bandfunc_win, band_r]
        plt.plot(wave_tot, bandfunc_tot, color='r', ls='-', marker='')

    plt.axvline(ln_ctr, color='r', ls=':', linewidth=1)
    #plt.axvline(wmin, color='r', ls='-', linewidth=1)
    #plt.axvline(wmax, color='r', ls='-', linewidth=1)

    plt.xlim(wave_reg[0], wave_reg[-1])
    plt.ylim(-0.1, 1.3)

    plt.ylabel('Normalized {} flux'.format(ln_id))
    plt.xlabel('Wavelength [Ang]')

    if file_type in ac_set.ftypes['1d']:
        plt.annotate("sp. median S/N = {:.2f}".format(median_snr),
                        xy=(0.05,0.9),
                        xycoords='axes fraction',
                        textcoords='axes fraction')

    if file_type in ac_set.ftypes['2d']:
        plt.annotate("S/N = {:.2f}".format(snr),
                        xy=(0.05,0.9),
                        xycoords='axes fraction',
                        textcoords='axes fraction')

    # For case of reading from rdb, obj and date will be huge repeating lists
    if type(obj) == list:
        obj = obj[0]
        obs_date = obs_date[0]

    # Save plots
    if out_dir and out_dir != 'show':
        print("Saving plot of line {}".format(ln_id))
        if ord is not None:
            save_name = '{}_{}_{}_{}_{}_ord{}_{}'.format(obj, obs_date, instr, file_type, ln_id, ord, ac_set.fnames['ln_plt'])
        else:
            save_name = '{}_{}_{}_{}_{}_{}'.format(obj, obs_date, instr, file_type, ln_id, ac_set.fnames['ln_plt'])

        # Create directories
        if not os.path.isdir(out_dir):
            os.makedirs(out_dir)
        if not os.path.isdir(os.path.join(out_dir, obj)):
            os.mkdir(os.path.join(out_dir, obj))

        plt.savefig(os.path.join(out_dir, obj, save_name))

    # Show plots
    if out_dir == 'show': plt.show()

    plt.close()

    return


def get_win(data, ln_id, ln_ctr, ln_win, bandtype, blaze=None, err=None, frac=True, ln_plts=False):
    """
    Calculate the flux for each line.
    """

    file_type  = data['file_type']
    wave       = data['wave']
    flux       = data['flux']
    obs_date   = data['obs_date']
    obj        = data['obj']
    snr        = data['snr']
    median_snr = data['median_snr']
    blaze      = data['blaze']
    noise      = data['noise']

    # 2d spectra
    if type(wave[0]) is np.ndarray:

        # NOTE: sel_ord is using ln_win/2 (as in sq bandpass)
        orders, ord_type = sel_order(wave, ln_ctr, ln_win, file_type)

        num     = []
        denum   = []
        snr_ord = []
        for k in range(len(orders)):
            ord = orders[k]
            print('Order = {}'.format(ord))
            print('SNR = {:.2f}'.format(snr[ord]))

            wave_ord  = wave[ord]
            flux_ord  = flux[ord]
            blaze_ord = blaze[ord]

            snr_ord.append(snr[ord])

            # Compute flux for line parameters
            f_sum, f_sum_var, bandfunc, npixels, flg, frac_neg = ac_tools.compute_flux(wave_ord, flux_ord, blaze_ord, noise, ln_ctr, ln_win, bandtype=bandtype, frac=frac)

            # Plot line regions
            if ln_plts:
                line_plot(data, ord, ln_id, ln_ctr, ln_win, bandfunc, bandtype, out_dir=ln_plts)

            num.append(f_sum/f_sum_var)
            denum.append(1.0/f_sum_var)

        # Weighted sum of fluxes in the orders
        f_sum = sum(num)/sum(denum)

        # Compute error
        f_sum_var = 1.0/sum(denum)
        if f_sum_var >= 0.0:
            f_sum_err = np.sqrt(f_sum_var)
        else:
            print("*** ERROR: Variance of flux in the line is negative. Cannot compute error.")
            f_sum_err = 0.0

        # mMdian snr of the used orders
        snr = np.median(snr_ord)

    # 1d spectra
    elif type(wave[0]) in [np.float, np.float64]:
        ord      = None
        ord_type = None
        wave     = wave
        flux     = flux

        # Compute flux for line parameters
        f_sum, f_sum_var, bandfunc, npixels, flg, frac_neg = ac_tools.compute_flux(wave, flux, blaze, noise, ln_ctr, ln_win, bandtype=bandtype, frac=frac)

        # Compute error
        if f_sum_var >= 0.0:
            f_sum_err = np.sqrt(f_sum_var)
        else:
            print("*** ERROR: Variance of flux in the line is negative. Cannot compute error.")
            f_sum_err = 0.0

        # Plot line regions
        if ln_plts:
            line_plot(data, ord, ln_id, ln_ctr, ln_win, bandfunc, bandtype, out_dir=ln_plts)

    win = {}
    win['sum']      = f_sum
    win['sum_err']  = f_sum_err
    win['flg']      = flg
    win['frac_neg'] = frac_neg
    win['snr']      = snr
    win['order']    = ord
    win['bandfunc'] = bandfunc
    win['ord_type'] = ord_type
    win['npixels']  = npixels

    return win
