#!/usr/bin/env python

import numpy as np
import pylab as plt

### added to compare with malavolta
from astropy.io import fits

# SPECHA FILES:
import actin_functions as func


# THIS FILE USES R. DIAZ METHOD TO CALCULATE INDICES


def integrate(x, y, weight, xmin, xmax, vary=None):
    """
    Return integrated flux of array y (already integrated in
    each pixel) with respecto to x, with variance vary and
    weight function weigth between limits xmin and xmax.
    """

    # Keep only elements of x, y, and weight within the interval,
    # include fractional pixels
    deltax = np.diff(x)
    # Reduce x of first element for compatibility with deltax
    x = x[1:]

    cond = (x > xmin - deltax/2.0) & (x <= xmax + deltax/2.0)

    # Compute fraction of pixel within interval
    fraction_left = 0.5 - (xmin - x[cond])/deltax[cond]
    fraction_right = 0.5 - (x[cond] - xmax)/deltax[cond]

    fraction_left = np.where(fraction_left > 1, 1.0, fraction_left)
    fraction_right = np.where(fraction_right > 1, 1.0, fraction_right)
    #print "fraction_right = %s" % fraction_right

    fraction = np.minimum(fraction_left, fraction_right)
    print "N = %s" % (len(y[1:][cond]))

    # Sum contributions of pixels inside interval, considering fractions
    # when necessary
    summed_y = np.sum(y[1:][cond] * weight[1:][cond] * fraction)
    summed_weight = np.sum(weight[1:][cond] * fraction)

    integral = np.divide(summed_y, summed_weight)

    if vary is not None:
        # Also compute error
        summed_var = np.sum(vary[1:][cond] * (weight[1:][cond] * fraction)**2)
        var_integral = np.divide(summed_var, summed_weight**2)

    else:
        var_integral = None

    # Flag values with negative flux, and gives fraction of pixels w/ neg values
    flg, frac_n_neg = func.flag_negflux(y[cond])

    return integral, var_integral, flg, frac_n_neg


def compute_flux(e2ds, blaze, weight, ww, wwmin, wwmax, noise=0.0):

    # Find the spectral order(s) where the relevant window is located
    orders = np.unique(np.where(np.logical_and(ww > wwmin, ww < wwmax))[0])
    print "Found %i orders for this line: %s" % (len(orders),orders)

    # ADDED TO COMPARE WITH MY CODE
    #import actin_get_win as get_win
    #orders = get_win.sel_order(ww, (wwmin + (wwmax - wwmin)/2.), (wwmax - wwmin))
    #print orders
    #orders = [orders[-1]] ####
    #print orders
    #######

    # Initialize arrays
    flux = np.zeros(len(orders))
    vflux = np.zeros(len(orders))

    for i, order in enumerate(orders):
        #
        e2ds_order = e2ds[order]
        blaze_order = blaze[order]
        ww_order = ww[order]
        weight_order = weight[order]

        # Normalize
        e2dsn = e2ds_order/blaze_order

        # Compute variance in e2dsn (error is sqrt(e2ds))
        var_e2dsn = (e2ds_order + noise**2)/blaze_order**2

        flux[i], vflux[i], flg, frac_n_neg = integrate(ww_order, e2dsn, weight_order, wwmin, wwmax, vary=var_e2dsn)


    # Obtain weighted average of fluxes, if present in more than one order
    ff = np.sum(flux/vflux)/np.sum([1.0/vflux])
    vv = 1.0/np.sum([1.0/vflux])

    return ff, vv, flux, orders, flg, frac_n_neg


def get_win(flux_2d, blaze_2d, wave_2d, snr_2d, ln_ctr, ln_win, ln_c, bandtype):

    win = {}

    min_wave = ln_ctr-ln_win/2.
    max_wave = ln_ctr+ln_win/2.

    e2ds = flux_2d
    blaze = blaze_2d
    #weight = np.ones([len(blaze_2d),len(blaze_2d[0])])
    weight = blaze_2d
    ww = wave_2d
    wwmin = min_wave
    wwmax = max_wave

    if bandtype == 'sq':
        f,f_var,flux,orders,flg, frac_n_neg  = compute_flux(e2ds, blaze, weight, ww, wwmin, wwmax, noise=0.0)

    if bandtype == 'tri':
        # Add triangular weight function for line:
        triang1 = np.where(ww >= ln_ctr,
                            -(ww - ln_ctr)/ln_win + 1,
                            (ww - ln_ctr)/ln_win + 1)

        triang1 = np.where(triang1 > 0, triang1, triang1*0.0)

        f,f_var,flux,orders,flg,frac_n_neg  = compute_flux(e2ds, blaze, weight*triang1, ww, wwmin, wwmax, noise=0.0)

    win['sum'] = f
    win['sum_err'] = np.sqrt(f_var)

    win['error_pixel'] = f_var

    win['flg'] = flg
    win['frac_neg'] = frac_n_neg
    win['flux'] = flux
    #win['flux_deb'] =
    #win['wave']
    win['snr'] = snr_2d[orders[0]]
    win['order'] = orders[0]


    return win
