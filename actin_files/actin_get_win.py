#!/usr/bin/env python

import numpy as np

### added to compare with malavolta
from astropy.io import fits

import pylab as plt
import actin_functions as func


def sel_order(wave_2d, ln_ctr, ln_win):
	"""
	Selects a spectral order for a line based and its selected bandwidth.

	Parameters:
	-----------
	wave_2d : list
		Wavelength in each spectral order [angstrom].
	ln_ctr : float
		Centre of the spectral line [angstrom].
	ln_win : float
		Bandwidth around the line where the flux is to be integrated
		[angstrom].

	Returns:
	--------
	order : {int, None}
		Spectral order for the given line parameters, None if no order is
		found.
	"""

	ord = None
	min_wave = ln_ctr-ln_win/2.
	max_wave = ln_ctr+ln_win/2.
	order = []
	for k in range(len(wave_2d)):
		if min_wave > wave_2d[k][0] and max_wave < wave_2d[k][-1]:
			ord = k
			order.append(ord)

	if ord == None:
		print ""
		print "ERROR: Could not determine spectral order for:"
		print "       * min_wave = %f.2" % min_wave
		print "       * max_wave = %f.2" % max_wave

	return order


def get_win(wave, flux, ln_ctr, ln_win, ln_c, bandtype, blaze=None, snr=None, err=None, weight=None):

	win = {}

	if blaze != None:
		order = sel_order(wave, ln_ctr, ln_win)
		print order
		order = order[-1] # Higher orders have higher SNR
		print order

		wave = np.asarray(wave[order])
		flux = np.asarray(flux[order])
		blaze = np.asarray(blaze[order])
		snr = np.asarray(snr[order])

	elif blaze == None:
		order = None
		wave = np.asarray(wave)
		flux = np.asarray(flux)
		blaze = np.ones(len(flux))
		snr = None
		#print "*** WARNING: actin_get_win: Blaze file not available."

	cond = (wave >= ln_ctr - ln_win/2.) & (wave <= ln_ctr + ln_win/2.)

	# Function to flag negative flux
	flg, frac_neg = func.flag_negflux(flux[cond])

	if flg == 'negFlux':
		print "*** WARNING: Negative flux detected"
		print "Fraction of pixels with negative flux = %s" % frac_neg

	# Computing flux for line parameters
	f_sum, f_sum_err, bandfunc = func.compute_flux(wave, flux, blaze, ln_ctr, ln_win, ln_c, bandtype=bandtype, weight=weight)

	win['sum'] = f_sum
	win['sum_err'] = f_sum_err

	win['flg'] = flg
	win['frac_neg'] = frac_neg

	win['snr'] = snr
	win['order'] = order

	win['bandfunc'] = bandfunc

	return win
