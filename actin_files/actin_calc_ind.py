#!/usr/bin/env python

import numpy as np


# ACTIN FILES
import actin_get_win as get_win
import actin_save_data as save_data

#import get_win_diaz

# add option to show flux in each line

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

	print "\nCHECKING LINES FOR WAVELENGTH RANGE AND SP. ORDERS"
	print "--------------------------------------------------"

	if wave == None or sel_lines == None: return

	# check if wave is 1d or 2d list
	try:
		if len(wave[0]): # 2d spec
			min_spec_wave = wave[0][0]
			max_spec_wave = wave[-1][-1]
			spec_type = '2d'
	except: # 1d spec
		min_spec_wave = wave[0]
		max_spec_wave = wave[-1]
		spec_type = '1d'

	# For each row (sp. line) in the config table calculate the min and max values of bandwidth
	rows = len(sel_lines['ln_id'])
	for k in range(rows):
		ln_id = sel_lines['ln_id'][k]
		ln_ctr = sel_lines['ln_ctr'][k]
		ln_win = sel_lines['ln_win'][k]

		min_wave = ln_ctr - ln_win/2.
		max_wave = ln_ctr + ln_win/2.

		# Check if line fits inside spectral range
		if min_wave < min_spec_wave or max_wave > max_spec_wave:
			print "*** ERROR: Line %s bandwidth outside spectral range" % ln_id
			exit()
		else: print "Line %s inside spectral range" % ln_id

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

			if order == None:
				print "*** ERROR: Could not determine sp. order for %s" % ln_id
				print "\tmin_wave = %.2f" % min_wave
				print "\tmax_wave = %.2f" % max_wave
				print "\tThe closest orders are:"
				for k in range(len(ln_ctr_orders)):
					closest_ord = ln_ctr_orders[k]
					print "\tOrder %i: %.2f-%.2f" % (closest_ord,wave[closest_ord][0],wave[closest_ord][-1])
				exit()
			else: print "Line %s inside spectral order %i" % (ln_id,order)
	return


def calc_flux_lines(data, sel_lines, save_plots=False, weight=None):

	print "\nCALCULATING FLUX IN SELECTED LINES"
	print "----------------------------------"

	if data == None:
		print "*** ERROR: actin_calc_ind: 'data' dictionary is None."
		return
	if sel_lines == None:
		print "*** ERROR: actin_calc_ind: 'sel_lines' dictionary is None."
		return

	wave = np.asarray(data['wave'])
	flux = np.asarray(data['flux'])
	try: blaze = np.asarray(data['blaze'])
	except: blaze = None
	obj = data['obj']
	date = data['date']

	snr = data['snr']

	sel_lines['flux'] = []
	sel_lines['error'] = []
	sel_lines['snr'] = []
	sel_lines['flg'] = []
	sel_lines['frac_neg'] = []

	ln_ctr = sel_lines['ln_ctr']
	ln_win = sel_lines['ln_win']
	ln_id = sel_lines['ln_id']
	ln_c = sel_lines['ind_c']
	ind_var = sel_lines['ind_var']
	bandtype = sel_lines['bandtype']

	if 'error_pixel' in data.keys():
		# case where pixel errors are given in rdb file as "error_pixel"
		print "Using pixel errors from input rdb file."
		err = data['error_pixel']
	else: err = None

	for k in range(len(ln_id)):
		print "\nComputing flux in line %s" % ln_id[k]
		print "-----------------------%s" % ('-'*len(ln_id[k]))

		print "Using %s bandwidth" % bandtype[k]

		win = get_win.get_win(wave, flux, ln_ctr[k], ln_win[k], ln_c[k], bandtype[k], blaze=blaze, snr=snr, err=err, weight=weight)

		# Diaz method:
		#win = get_win_diaz.get_win(flux, blaze, wave, snr, ln_ctr[k], ln_win[k], ln_c[k], bandtype[k])


		sel_lines['flux'].append(win['sum'])
		sel_lines['error'].append(win['sum_err'])
		sel_lines['snr'].append(win['snr'])
		sel_lines['flg'].append(win['flg'])
		sel_lines['frac_neg'].append(win['frac_neg'])

		# Save plots in the line regions:
		if save_plots != False:
			print "Saving plot of line %s" % ln_id[k]

			bandfunc = win['bandfunc']

			if win['order'] != None:
				file_type = 'e2ds'
				wave_plt = np.asarray(wave[win['order']])
				flux_plt = np.asarray(flux[win['order']])
				blaze_plt = np.asarray(blaze[win['order']])
				if win['bandfunc'] != None:
					bandfunc_plt = np.asarray(bandfunc[win['order']])
				else: bandfunc_plt = None

				flux_plt = flux_plt/blaze_plt
			elif win['order'] == None:
				file_type = 's1d'
				wave_plt = wave
				flux_plt = flux

			save_data.line_plot(wave_plt, flux_plt, obj, date, ln_id[k], ln_ctr[k], ln_win[k], bandfunc, file_type, out_dir=save_plots)

	return sel_lines


def calc_ind(sel_lines):

	print "\nCALCULATING INDICES"
	print "-------------------"

	if sel_lines == None:
		print "*** ERROR: actin_calc_ind: sel_lines is None."
		return

	# remove duplicates of ind_id and gives a list of selected indices
	sel_ind = list(set(sel_lines['ind_id']))

	index = {}
	index['index'] = []
	index['value'] = []
	index['error'] = []
	index['flg'] = []
	index['mfrac_neg'] = []
	index['snr'] = []

	print "index\tvalue\terror\t\tsnr\tflag\tmfrac_neg"
	print "-----\t-----\t-----\t\t---\t----\t---------"

	ind_ids = sel_lines['ind_id']
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

		# sum of percentage of pixels with neg flux in each line
		mfrac_neg = max(frac_neg)

		if "negFlux" in flg: flg_ind = 'negFlux'
		else: flg_ind = None

		# calc mean snr of lines:
		try: snr_ind = np.mean([float(snr[k]) for k in range(len(snr))])
		except: snr_ind = None

		# add line variables for numerator or denominator:
		num = [flux[k] for k in range(len(var)) if 'L' in var[k]]
		num_err = [err[k] for k in range(len(var)) if 'L' in var[k]]
		denom = [flux[k] for k in range(len(var)) if 'R' in var[k]]
		denom_err = [err[k] for k in range(len(var)) if 'R' in var[k]]

		num = np.asarray(num)
		denom = np.asarray(denom)
		num_err = np.asarray(num_err)
		denom_err = np.asarray(denom_err)

		ind = sum(num) / sum(denom)

		# error using propagation of errors for lines and ref lines
		ind_err = np.sqrt(sum(num_err**2) + ind**2 * sum(denom_err**2))/sum(denom)

		index['index'].append(sel_ind[i])
		index['value'].append(ind)
		index['error'].append(ind_err)
		index['flg'].append(flg_ind)
		index['mfrac_neg'].append(mfrac_neg)
		index['snr'].append(snr_ind)

		print "%s\t%.4f\t%.7f\t%s\t%s\t%.6f" % (index['index'][i], index['value'][i], index['error'][i],index['snr'][i],index['flg'][i], index['mfrac_neg'][i])

	# 	if sel_ind[i] == 'I_CaII':
	# 		# Calibration of I_CaII to S_MW scale using Baliunas stars
	# 		a = 0.5506850 # blaze, 0.6
	# 		b = 0.0857112 # blaze, 0.6
	# 		s_mw = a*ind + b
	# 		s_mw_err = a*ind_err
	# 		s_mw_flg = flg_ind
	# 		s_mw_mfrac_neg = mfrac_neg
	# 		s_mw_snr = snr_ind
	# 	else: pass
    #
	# if 'I_CaII' in sel_ind:
	# 	index['index'].append('s_mw')
	# 	index['value'].append(s_mw)
	# 	index['error'].append(s_mw_err)
	# 	index['flg'].append(s_mw_flg)
	# 	index['mfrac_neg'].append(s_mw_mfrac_neg)
	# 	index['snr'].append(s_mw_snr)
	# 	print "%s\t%.4f\t%.7f\t%s\t%s\t%.6f" % ('s_mw', s_mw, s_mw_err, s_mw_snr, s_mw_flg, s_mw_mfrac_neg)

	return index
