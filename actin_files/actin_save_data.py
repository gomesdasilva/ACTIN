#!/usr/bin/env python

import sys,os
import pylab as plt
import numpy as np

# ACTIN FILES:
import pyrdb
import actin_functions as func


def check_duplicate(obj, date, file_type, out_dir):
	"""
	Check if measurement is a duplicate in the output file.

	Parameters:
	-----------
	obj : str
		Object (target) identification.
	date : str
		Date of observation in the fits file format.
	output_dir : str, optional
		Directory of output file, 'output' is default.

	Returns:
	--------
	bool
		True if duplicate, False otherwise.
	"""

	print "\nCHECKING DUPLICATE"
	print "------------------"

	if obj == None or date == None or out_dir == None: return

	file_pname = '%s/%s/%s_%s_actin.rdb' % (out_dir, obj, obj, file_type)
	if os.path.isfile(file_pname):
		
		rdb_data = func.read_rdb(file_pname)[0]

		if date in rdb_data['date']:
			print "Date %s already saved in %s" % (date,file_pname)
			print "*** ACTION: Ignoring measurement"
			return True
		else:
			print "Date %s not present in %s" % (date,file_pname)

	return False


def line_plot(wave, flux, obj, date, ln_id, ln_ctr, ln_win, bandfunc, file_type, out_dir):

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

	if bandfunc != None:
		#plt.plot(wave[(bandfunc > 0)], flux[(bandfunc > 0)]/max(flux),'r-')
		plt.plot(wave[(bandfunc > 0)], bandfunc[(bandfunc > 0)],'r-', linewidth=1.5)

	plt.xlim(wave_min, wave_max)

	plt.ylim(0., 1.3)#max(flux_reg) + max(flux_reg) / 4.)
	#plt.ylim(-0.3,1.5)

	if bandfunc == None:
		#plt.plot(wave_win, flux_win/max(flux),'r-')
		plt.axvline(ln_ctr, color='r', ls='--', linewidth=1.5)
		plt.axvline(win_min, color='r', ls='-', linewidth=1.5)
		plt.axvline(win_max, color='r', ls='-', linewidth=1.5)

	plt.ylabel('%s flux' % ln_id)
	plt.xlabel('Wavelength [Ang]')

	# For case of reading from rdb, obj and date will be huge repeating lists
	if type(obj) == list:
		obj = obj[0]
		date = date[0]

	save_name = '%s_%s_%s_win.pdf' % (date, file_type, ln_id)

	if not os.path.isdir(out_dir):
		os.makedirs(out_dir)
	if not os.path.isdir('%s/%s' % (out_dir, obj)):
		os.mkdir('%s/%s' % (out_dir, obj))

	path = '%s/%s' % (out_dir, obj)
	plt.savefig('%s/%s' % (path, save_name))
	plt.close()

	return


def save_data(data, index, out_dir):
	"""
	Saves data from specha to rdb.

		- Automatically detects if reading from 1d or 2d fits files.
		- Automatically detects number of selected indices.

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
		value 		float : Computed value for the spectral index.
		error 		float : Error of the index calculated by error propa-
					gation.
		snr 		float, None : Mean of the SNR at the lines spectral
					order if the SNR per order was given as input, median
					SNR of spectrum if median SNR given as input, 'None' if
					no SNR values given as input.
		flg 		str, None : Flags associated with the index: 'negFlux'
					if negative flux detected in any of the lines used to
					compute the index, None otherwise.
		==========  ========================================================
	output_dir : str, optional
		Directory of the saved rdb file, 'output' is default.
	"""

	print "\nSAVING DATA"
	print "-----------"

	if data == None or out_dir == None: return

	# For case of reading from rdb, obj and date will be huge repeating lists
	# reading from rdb file
	if type(data['obj']) == list:
		data['obj'] = data['obj'][0]
		data['date'] = data['date'][0]
		data['bjd'] = data['bjd'][0]
		data_keys = ['obj','date','bjd']
	else:
		data_keys = ['obj', 'date', 'bjd', 'median_snr', 'data_flg']

	if not os.path.isdir(out_dir):
		os.makedirs(out_dir)
	if not os.path.isdir('%s/%s' % (out_dir, data['obj'])):
		os.mkdir('%s/%s' % (out_dir, data['obj']))

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

	save_name = "%s/%s/%s_%s_actin.rdb" % (out_dir, all_data['obj'], all_data['obj'], data['file_type'])

	if not os.path.isfile(save_name):
		func.save_rdb(data_save, keys, save_name)
		print "%s | %s | Data saved to %s" % (all_data['obj'], all_data['date'], save_name)
	else:
		func.add_rdb(data_save, keys, save_name)
		print "%s | %s | Data added to %s" % (all_data['obj'], all_data['date'], save_name)

	return save_name
