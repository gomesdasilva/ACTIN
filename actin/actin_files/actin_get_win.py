#!/usr/bin/env python

# compatibility with python 2/3:
from __future__ import print_function
from __future__ import division

import numpy as np

import actin_functions as func


def sel_order(wave_2d, ln_ctr, ln_win):
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

	ord = None
	min_wave = ln_ctr-ln_win/2.
	max_wave = ln_ctr+ln_win/2.
	order = []
	for k in range(len(wave_2d)):
		if min_wave > wave_2d[k][0] and max_wave < wave_2d[k][-1]:
			ord = k
			order.append(ord)

	if ord is None:
		print("\nERROR: Could not determine spectral order for:")
		print("       * min_wave = %f.2" % min_wave)
		print("       * max_wave = %f.2" % max_wave)

	return order


def get_win(wave, flux, ln_ctr, ln_win, ln_c, bandtype, blaze=None, snr=None, err=None, weight=None, norm='npixels'):
	"""
	Calculates the sum of the flux and associated error for a given spectral line.

		- Flags regions where negative flux is detected.
		- For .rdb file: detects if errors in flux per pixel are given,
		  otherwise calculates them.

	Parameters:
	-----------
	wave : list of list of lists
		Wavelength per pixel (1d spectrum) or wavelength per pixel per order
		(2d spectrum) [angstroms].
	flux : list or list of lists
		Flux per pixel (1d spectrum) of Flux per pixel per order (2d spectrum).
		wave : list of list of lists
	ln_ctr : float
		Centre of the line [angstroms]
	ln_win : float
		Bandpass around the line where the flux is to be integrated
		[angstroms].
	ln_c : float
		Constant multiplied to the integrated flux.
	bandtype : string
		Function to be used in the integration of flux. If 'sq', uses a square function with limits 'ln_win', if 'tri' a triangular function with full-width-at-half-maximum given by 'ln_win'.
	blaze : list (optional)
		Blaze function. None is default.
	snr : {list, None} (optional)
		SNR per spectral order. None is default.
	err : {list, None} (optional)
		List of errors on the flux per pixel, None if not available (default).
	weight : {str, None} (optional)
		Function to weight the integrated flux. If 'blaze' the flux is multiplied by the blaze function, if None the flux is not weighted (default).
	norm : str (optional)
		Normalisation of the flux: if 'band' the sum is normalised by the bandpass wavelength value in angstroms, if 'npixels' by the number of pixels in the bandpass (default), if 'weight' by the sum of the weight function inside the bandpass, if None the integrated flux is not normalised.

	Returns:
	--------
	win : dict
		Dictionary of parameters returned from computing the flux in the
		regions around the selected lines.

		The returned keys are:

		==========  ========================================================
		keys		Description
		----------  --------------------------------------------------------
		sum 		float : Integrated flux inside the bandwidth.
		sum_err 	float : Error on the integrated flux, either calculated
					using the errors given in the input option or using
					Poisson noise.
					NOTE: When using normalized flux (as in the
					case of s1d files) the errors are only indicative and
					do not represent real estimates on the error. Errors on
					the flux should always use real flux.


		flg 		{str, None} : Flags the presence of negative fluxes
					inside the bandwidth as 'negFlux', None otherwise.
		frac_neg 	float : Fraction of pixels with negative values of flux.
		snr			float : SNR in the spectral order used.
		order		{int, None} : Spectral order used. None if 1d
					spectrum.
		bandfunc	{list, None} : Bandpass function, if triangular, None if
					square.
		norm		str : Normalization option used.
		==========  ========================================================
	"""

	win = {}

	if blaze is not None:
		order = sel_order(wave, ln_ctr, ln_win)
		order = order[-1]
		print("Using order %i" % order)

		wave = np.asarray(wave[order])
		flux = np.asarray(flux[order])
		blaze = np.asarray(blaze[order])
		snr = np.asarray(snr[order])

	elif blaze is None:
		order = None
		wave = np.asarray(wave)
		flux = np.asarray(flux)
		blaze = np.ones(len(flux))
		snr = None

	wmin = ln_ctr - ln_win/2.
	wmax = ln_ctr + ln_win/2.

	# Calculate the flux inside bandpass taking into account fraction of pixels
	flux_win, pixels_win = func.frac_pixels(wave, flux, wmin, wmax)

	print("Nr. of pixels in bandpass = %.2f" % pixels_win)

	# Function to flag negative flux inside
	flg, frac_neg = func.flag_negflux(flux_win)

	if flg == 'negFlux':
		print("*** WARNING: Negative flux detected")
		print("Fraction of flux with negative values = %s" % frac_neg)

	# Computing flux for line parameters
	f_sum, f_sum_err, bandfunc = func.compute_flux(wave, flux, blaze, ln_ctr, ln_win, ln_c, bandtype=bandtype, weight=weight,norm=norm)

	win['sum'] = f_sum
	win['sum_err'] = f_sum_err

	win['flg'] = flg
	win['frac_neg'] = frac_neg

	win['snr'] = snr
	win['order'] = order

	win['bandfunc'] = bandfunc
	win['norm'] = norm

	return win
