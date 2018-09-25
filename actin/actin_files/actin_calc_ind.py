#!/usr/bin/env python

# compatibility with python 2/3:
from __future__ import print_function
from __future__ import division

import numpy as np

# ACTIN FILES
import actin_get_win as get_win
import actin_save_data as save_data


def check_lines(wave, sel_lines):
	"""
	Tests if the selected lines from config file fit inside the spectrograph
	wavelength range and fit inside any spectral orders for the case of 2d
	spectrum.

	Parameters:
	-----------
	wave : list of floats (1d) or list of lists of floats (2d)
		Wavelength in 1d or 2d format [angstrom].
	sel_lines : dict
		Dictionary containing the identification of the indices selected and
		the parameters of the spectral lines required for each index.

		Each key entry is a list of parameters where the list indices form the
		rows related to the same spectral line identified with the key 'ln_id'
		which is related to the spectral index identified by 'ind_id'.

		The used keys are:

		==========  ========================================================
		keys		Description
		----------  --------------------------------------------------------
		ln_id		str : Spectral line identification.
		ln_ctr 		float : Wavelength of the line centre [angstrom].
		ln_win 		float : Bandwidth around the line centre to be used in
					the flux integration [angstrom].
		==========  ========================================================

	Returns
	-------
	exit() if test fails.
	"""

	print("\nCHECKING LINES FOR WAVELENGTH RANGE AND SP. ORDERS")
	print("--------------------------------------------------")

	if wave is None or sel_lines is None: return

	if type(wave[0]) is np.ndarray: # 2d spec
		min_spec_wave = wave[0][0]
		max_spec_wave = wave[-1][-1]
		spec_type = '2d'
	elif type(wave[0]) is np.float64: # 1d spec
		min_spec_wave = wave[0]
		max_spec_wave = wave[-1]
		spec_type = '1d'

	# For each row (sp. line) in the config table calculate the min and max values of bandwidth
	rows = len(sel_lines['ln_id'])
	for k in range(rows):
		ln_id = sel_lines['ln_id'][k]
		ln_ctr = sel_lines['ln_ctr'][k]
		ln_win = sel_lines['ln_win'][k]

		if ln_win == 0:
			print("*** ERROR: line %s bandwidth is zero" % ln_id)
			quit()

		min_wave = ln_ctr - ln_win/2.
		max_wave = ln_ctr + ln_win/2.

		# Check if line fits inside spectral range
		if min_wave < min_spec_wave or max_wave > max_spec_wave:
			print("*** ERROR: Line %s bandwidth outside spectral range" % ln_id)
			quit()
		else: print("Line %s inside spectral range" % ln_id)

		# If wave is 2d check if line fits inside sp. order
		if spec_type == '2d':
			order = None
			ln_ctr_orders = []
			for i in range(len(wave)):
				if min_wave > wave[i][0] and max_wave < wave[i][-1]:
					order = i
				# used to show potencial orders with wavelength range
				if ln_ctr > wave[i][0] and ln_ctr < wave[i][-1]:
					ln_ctr_orders.append(i)

			if order is None:
				print("*** ERROR: Could not determine sp. order for %s" % ln_id)
				print("\tmin_wave = %.2f" % min_wave)
				print("\tmax_wave = %.2f" % max_wave)
				print("\tThe closest orders are:")
				for k in range(len(ln_ctr_orders)):
					closest_ord = ln_ctr_orders[k]
					print("\tOrder %i: %.2f-%.2f" % (closest_ord,wave[closest_ord][0],wave[closest_ord][-1]))
				quit()
			else: print("Line %s inside spectral order %i" % (ln_id,order))
	return


def calc_flux_lines(data, sel_lines, save_plots=False, weight=None, norm='npixels'):
	"""
	Calculates the sum of the flux and associated errors for all spectral lines required to calculate the selected indices.

	Parameters:
	-----------
	data : dict
		Dictionary with data returned from fits files.

		The required keys are:

		==========  ========================================================
		keys		Description
		----------  --------------------------------------------------------
		flux 		list of lists : Flux per pixel per order.
		wave 		list of lists : Wavelength calibrated for BERV and RV
					(at rest frame) per pixel per order [angstroms].
		blaze		list of lists : Blaze function.
		snr 		list : SNR at each spectral order.
		obj 		str : Object identification.
		date 		str : Date of observation in the fits file format.
		==========  ========================================================

	sel_lines : dict
		Dictionary containing the identification of the indices selected and
		the parameters of the spectral lines required for each index.

		The required keys are:

		==========  ========================================================
		keys		Description
		----------  --------------------------------------------------------
		ind_id		str : Index identification.
		ind_var		str : Variable assigned to a given line to be used in
					the index equation. Ex: 'L1', 'L2', etc, for the core
					lines, and 'R1', 'R2', etc, for reference lines.
		ln_id		str : Spectral line identification.
		ln_c		float : Constant to be multilied to the flux of the line.
		ln_ctr 		float : Wavelength of the line centre [angstroms].
		ln_win 		float : Bandpass around the line centre to be used in
					the flux integration [angstroms].
		bandtype    str : Bandpass type to integrate flux.
		==========  ========================================================

	save_plot : {str, False} (optional)
		Path for the saved line plots. If False, the plots are not saved
		(default).
	weight : {str, None} (optional)
		Function to weight the integrated flux. If 'blaze' the flux is multiplied by the blaze function, if None (default) the flux is not weighted (default).
	norm : str (optional)
		Normalisation of the flux: if 'band' the sum is normalised by the bandpass wavelength value in angstroms, if 'npixels' by the number of pixels in the bandpass (default), if 'weight' by the sum of the weight function inside the bandpass, if None the integrated flux is not normalised.

	Returns:
	--------
	sel_lines : dict
		Dictionary containing the identification of the indices selected and
		the parameters of the spectral lines required for each index.

		Included in the output are also the keys used as input (see above).

		The returned keys are:

		==========  ========================================================
		keys		Description
		----------  --------------------------------------------------------
		flux		list : Integrated flux for each line.
		error		list : Errors on integrated flux for each line.
		snr			list : SNR for each line.
		flg			list : Flag for each line.
		frac_neg	list : Fraction of pixels with negative values of flux
					for each line.
		==========  ========================================================
	"""

	print("\nCALCULATING FLUX IN SELECTED LINES")
	print("----------------------------------")

	if data is None:
		print("*** ERROR: data dictionary is empty")
		return
	if sel_lines is None:
		print("*** ERROR: sel_lines dictionary is empty")
		return

	wave = np.asarray(data['wave'])
	flux = np.asarray(data['flux'])
	if 'blaze' in data.keys(): blaze = np.asarray(data['blaze'])
	else: blaze = None
	obj = data['obj']
	date = data['date']

	if 'snr' in data.keys(): snr = data['snr']
	else: snr = None

	sel_lines['flux'] = []
	sel_lines['error'] = []
	sel_lines['snr'] = []
	sel_lines['flg'] = []
	sel_lines['frac_neg'] = []

	ln_ctr = sel_lines['ln_ctr']
	ln_win = sel_lines['ln_win']
	ln_id = sel_lines['ln_id']
	ln_c = sel_lines['ln_c']
	ind_var = sel_lines['ind_var']
	bandtype = sel_lines['bandtype']

	if 'error_pixel' in data.keys():
		# case where pixel errors are given in rdb file as "error_pixel"
		print("Using pixel errors from input rdb file.")
		err = data['error_pixel']
	else: err = None

	for k in range(len(ln_id)):
		print("\nComputing flux in line %s" % ln_id[k])
		print("-----------------------%s" % ('-'*len(ln_id[k])))

		print("Using %s bandpass" % bandtype[k])

		win = get_win.get_win(wave, flux, ln_ctr[k], ln_win[k], ln_c[k], bandtype[k], blaze=blaze, snr=snr, err=err, weight=weight, norm=norm)

		sel_lines['flux'].append(win['sum'])
		sel_lines['error'].append(win['sum_err'])
		sel_lines['snr'].append(win['snr'])
		sel_lines['flg'].append(win['flg'])
		sel_lines['frac_neg'].append(win['frac_neg'])

		# Save plots in the line regions:
		if save_plots:
			print("Saving plot of line %s" % ln_id[k])

			bandfunc = win['bandfunc']

			if win['order'] is not None:
				file_type = 'e2ds'
				wave_plt = np.asarray(wave[win['order']])
				flux_plt = np.asarray(flux[win['order']])
				blaze_plt = np.asarray(blaze[win['order']])
				if win['bandfunc'] is not None:
					bandfunc_plt = np.asarray(bandfunc[win['order']])
				else: bandfunc_plt = None

				flux_plt = flux_plt/blaze_plt
			elif win['order'] is None:
				file_type = 's1d'
				wave_plt = wave
				flux_plt = flux

			save_data.line_plot(wave_plt, flux_plt, obj, date, ln_id[k], ln_ctr[k], ln_win[k], bandfunc, file_type, out_dir=save_plots)

	return sel_lines


def calc_ind(sel_lines):
	"""
	Calculates the indices identified in sel_lines as 'ind_id'.

	sel_lines : dict
		Dictionary containing the identification of the indices selected and
		the parameters of the spectral lines required for each index.

		The required keys are:

		==========  ========================================================
		keys		Description
		----------  --------------------------------------------------------
		ind_id		str : Index identification.
		ind_var		str : Variable assigned to a given line to be used in
					the index equation. Ex: 'L1', 'L2', etc, for the core
					lines, and 'R1', 'R2', etc, for reference lines.
		ln_id		str : Spectral line identification.
		flux		list : Integrated flux for each line.
		error		list : Errors on integrated flux for each line.
		snr			list : SNR for each line.
		flg			list : Flag for each line.
		frac_neg	list : Fraction of pixels with negative values of flux
					for each line.
		==========  ========================================================

	Returns:
	--------
	index : dict
		Dictionary containing the index values, errors and related
		information.

		The returned keys are:

		==========  ========================================================
		keys		Description
		----------  --------------------------------------------------------
		index		str : Index identification.
		value		float : Index value.
		error		float : Index error.
		flg			{str, None} : Index flag.
		mfrac_neg	float : Maximum fraction of pixels with negative values
					of flux taking into account all lines used to compute
					index.
		snr			float : Mean SNR of all lines used to compute index.
		==========  ========================================================
	"""

	print("\nCALCULATING INDICES")
	print("-------------------")

	if sel_lines is None:
		print("*** ERROR: sel_lines dictionary is empty.")
		return

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

		# Maximum fraction of flux with negative values of all lines
		mfrac_neg = max(frac_neg)

		if "negFlux" in flg: flg_ind = 'negFlux'
		else: flg_ind = None

		# Mean snr of lines:
		try: snr_ind = np.mean([float(snr[k]) for k in range(len(snr))])
		except: snr_ind = None

		for k in range(len(var)):
			if 'L' not in var[k] and 'R' not in var[k]:
				print("*** ERROR: 'ind_var' variable (in config file config_lines.txt) must start with either an 'L' for core line or 'R' for reference line. Value given was '%s'" % var[k])
				quit()

		# Add line variables for numerator or denominator:
		num = [flux[k] for k in range(len(var)) if 'L' in var[k]]
		num_err = [err[k] for k in range(len(var)) if 'L' in var[k]]
		denom = [flux[k] for k in range(len(var)) if 'R' in var[k]]
		denom_err = [err[k] for k in range(len(var)) if 'R' in var[k]]

		num = np.asarray(num)
		denom = np.asarray(denom)
		num_err = np.asarray(num_err)
		denom_err = np.asarray(denom_err)

		ind = sum(num) / sum(denom)

		# Error using propagation of errors for lines and ref lines
		ind_err = np.sqrt(sum(num_err**2) + ind**2 * sum(denom_err**2))/sum(denom)

		index['index'].append(sel_ind[i])
		index['value'].append(ind)
		index['error'].append(ind_err)
		index['flg'].append(flg_ind)
		index['mfrac_neg'].append(mfrac_neg)
		index['snr'].append(snr_ind)

		print("%s\t%.4f\t%.7f\t%s\t%s\t%.6f" % (index['index'][i], index['value'][i], index['error'][i],index['snr'][i],index['flg'][i], index['mfrac_neg'][i]))

	return index
