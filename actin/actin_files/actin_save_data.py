#!/usr/bin/env python

# compatibility with python 2/3:
from __future__ import print_function
from __future__ import division

import sys, os
import pylab as plt
import numpy as np

# ACTIN FILES:
import actin_functions as func


def check_duplicate(obj, date, file_type, out_dir):
	"""
	Check if measurement is a duplicate in the output file.

	Parameters:
	-----------
	obj : str
		Object identification.
	date : str
		Date of observation in the fits file format.
	file_type : str
		Type of file used: 'e2ds', 's1d', 'ADP', or 'rdb'.
	output_dir : str
		Directory of output file.

	Returns:
	--------
	bool
		True if duplicate, False otherwise.
	"""

	print("\nCHECKING DUPLICATE")
	print("------------------")

	if obj is None or date is None or out_dir is None: return

	file_name = "%s_%s_actin.rdb" % (obj, file_type)
	file_pname = os.path.join(out_dir, obj, file_name)
	if os.path.isfile(file_pname):

		rdb_data = func.read_rdb(file_pname)[0]

		if date in rdb_data['date']:
			print("Date %s already saved in %s" % (date,file_pname))
			print("*** ACTION: Ignoring measurement")
			return True
		else:
			print("Date %s not present in %s" % (date,file_pname))

	return False


def line_plot(wave, flux, obj, date, ln_id, ln_ctr, ln_win, bandfunc, file_type, out_dir):
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
	bandfunc : {list, None}
		Bandpass function, None if using a square function, list if using triangular.
	file_type : str
		Identification of the type of file used. Options are 'e2ds', 's1d',
		'ADP', and 'rdb'.
	out_dir : str
		Output directory to save the plots.

	Returns:
	--------
	Plots of the regions around the spectral lines used to calculate the
	indices. The plots are saved in the directory: 'out_dir'/'obj'.
	"""

	width, height = func.plot_params(6, 3.5)

	win_min = ln_ctr - ln_win/2
	win_max = ln_ctr + ln_win/2

	cond = (wave >= win_min) & (wave <= win_max)
	wave_win = wave[cond]
	flux_win = flux[cond]

	wave_min = win_min - 4.
	wave_max = win_max + 4.

	flux_reg = flux[(wave > wave_min + 4) & (wave < wave_max + 4)]

	plt.figure(figsize=(width, height))
	plt.plot(wave, flux/max(flux_reg),'k-')

	if bandfunc is not None:
		plt.plot(wave[(bandfunc > 0)], bandfunc[(bandfunc > 0)],'b-', linewidth=1.5)

	plt.xlim(wave_min, wave_max)

	plt.ylim(0., 1.3)

	if bandfunc is None:
		plt.axvline(ln_ctr, color='b', ls='--', linewidth=1.5)
		plt.axvline(win_min, color='b', ls='-', linewidth=1.5)
		plt.axvline(win_max, color='b', ls='-', linewidth=1.5)

	plt.ylabel('%s flux' % ln_id)
	plt.xlabel('Wavelength [Ang]')

	# For case of reading from rdb, obj and date will be huge repeating lists
	if type(obj) == list:
		obj = obj[0]
		date = date[0]

	save_name = '%s_%s_%s_%s_win.pdf' % (obj, date, file_type, ln_id)

	if not os.path.isdir(out_dir):
		os.makedirs(out_dir)
	if not os.path.isdir(os.path.join(out_dir, obj)):
		os.mkdir(os.path.join(out_dir, obj))

	plt.savefig(os.path.join(out_dir, obj, save_name))
	plt.close()

	return


def save_data(data, index, out_dir):
	"""
	Saves output data to rdb file.

	Parameters:
	-----------
	data : dict
		Dictionary with data returned from fits files.

		Each key is a list with data related to a given measurement date
		given by the key 'date'.

		The used keys are:

		==========  ========================================================
		keys		Description
		----------  --------------------------------------------------------
		obj 		str : Object (target) identification.
		median_snr  float : Median SNR of spectrum.
		date 		str : Date of observation in the fits file format.
		bjd 		float : Barycentric Julian date of observation [days].
		rv			float : Radial velocity [m/s] (if CCF file available).
		rv_err		float : Error on radial velocity (photon noise) [m/s]
					(if CCF file available).
		fwhm		float : Full-Width-at-Half-Maximum of the CCF line
                    profile [m/s] (if CCF file available).
        cont		float : Contrast of the CCF line profile [%] (if CCF
                    file available).
        bis			float : Bisector Inverse Span of the CCF line profile
                    [m/s] (if BIS file available).
        noise		float : CCF noise [m/s] (if CCF file available).
		data_flg 	str : Flag with value 'noDeblazed' when the blaze file
					was not found (and flux_deb is real flux), None
					otherwise.
		==========  ========================================================

	index : dict
		Dictionary containing the parameters related to the calculated spectral
		indices.

		Each key entry is a list of parameters where the list indices form the
		rows related to a given spectral index identified by the key 'index'.

		The returned keys are:

		==========  ========================================================
		keys		Description
		----------  --------------------------------------------------------
		index 		str : Identification of the spectral index as given in
					the configuration file.
		value 		float : Spectral index value.
		error 		float : Error of the index calculated by error propa-
					gation.
		snr 		{float, None} : Mean of the SNR at the lines spectral
					order if the SNR per order was given as input, median
					SNR of spectrum if median SNR given as input, 'None' if
					no SNR values given as input.
		flg 		{str, None} : Flags associated with the index: 'negFlux'
					if negative flux detected in any of the lines used to
					compute the index, None otherwise.
		mfrac_neg	list : Maximum fraction of flux with negative values
					when taking into account all lines for a given index.
		==========  ========================================================

	output_dir : str
		Directory of the saved rdb file.

	Returns:
	--------
	Saves output to .rdb file with name and path given by 'save_name'.

	save_name : str
		Output filename with path.
	"""

	print("\nSAVING DATA")
	print("-----------")

	if data is None or out_dir is None: return

	# For case of reading from rdb, obj and date will be huge repeating lists
	# reading from rdb file
	if type(data['obj']) == list:
		data['obj'] = data['obj'][0]
		data['date'] = data['date'][0]
		data['bjd'] = data['bjd'][0]
		data_keys = ['obj', 'date', 'bjd']
	else:
		data_keys = ['obj', 'instr', 'date', 'bjd', 'rv', 'rv_err', 'fwhm', 'cont', 'bis', 'noise', 'median_snr', 'data_flg']

	if not os.path.isdir(out_dir):
		os.makedirs(out_dir)
	if not os.path.isdir(os.path.join(out_dir, data['obj'])):
		os.mkdir(os.path.join(out_dir, data['obj']))

	# Include index information
	if index != None:
		# convert index data from table list to dictionary
		indices = [index['index'][k] for k in range(len(index['index']))]
		ind = {}
		for k in range(len(indices)):
			ind[indices[k]] = index['value'][k]
			ind['%s_err' % indices[k]] = index['error'][k]
			ind['%s_snr' % indices[k]] = index['snr'][k]
			ind['%s_flg' % indices[k]] = index['flg'][k]
			ind['%s_mfrac_neg' % indices[k]] = index['mfrac_neg'][k]

		index_keys = ind.keys()
		index_keys.sort()

		keys = data_keys + index_keys

		# merging data and index dictionaries to use below
		all_data = dict(data)
		all_data.update(ind)
	else:
		keys = data_keys
		all_data = data

	# for save_rdb to work need to have the dic with lists
	data_save = {}
	for k in range(len(keys)):
		data_save[keys[k]] = [all_data[keys[k]]]

	save_file = "%s_%s_actin.rdb" % (all_data['obj'], data['file_type'])
	save_name = os.path.join(out_dir, all_data['obj'], save_file)

	if not os.path.isfile(save_name):
		func.save_rdb(data_save, keys, save_name)
		print("%s | %s | Data saved to %s" % (all_data['obj'], all_data['date'], save_name))
	else:
		func.add_rdb(data_save, keys, save_name)
		print("%s | %s | Data added to %s" % (all_data['obj'], all_data['date'], save_name))

	return save_name
