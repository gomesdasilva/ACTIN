#!/usr/bin/env python

# compatibility with python 2/3:
from __future__ import print_function
from __future__ import division

import os
import glob
import numpy as np
import astropy.io.fits as pyfits

# Constants
lightspeed = 299792458.0 # [m/s]


def read_file_adp(adp_pfile, obj_name=None):
	"""
	Reads ADP fits file and recovers all necessary files to compute
	spectrum (including CCF). To be used with load_data_adp function.

		Functionalities:
		- Recognises if instrument is HARPS or HARPS-N.
		- Ignores file if not science (ThAr flux).

	Parameters:
	-----------
	adp_pfile : str
		ADP fits filename with path.
	obj_name : str (optional)
		Name given to object that overrides the name from fits file. None is default.

	Returns:
	--------
	files : dict
		Dictionary with the file names required to compute the spectrum.

		The returned keys are:

		==========  ========================================================
		keys		Description
		----------  --------------------------------------------------------
		adp 		str : ADP fits filename with path.
		ccf 		str : CCF fits filename with path.
		bis			{str, None} : BIS fits filename with path if found, None
					otherwise.
		instr		str : Instrument identification.
		obs			str : Code related to instrument to be used in fits
					headers.
		obj 		str : Object identification.
		date 		str : Date of observation in the fits file format.
		==========  ========================================================
	"""

	print("\nREADING FILES FROM ADP FITS FILE")
	print("--------------------------------")

	if not adp_pfile:
		print("*** ERROR: ADP file with path required")
		return

	try: adp = pyfits.open(adp_pfile)
	except:
		print("*** ERROR: Cannot read %s" % adp_pfile)
		return

	folder, adp_file = os.path.split(adp_pfile)

	if not folder: folder = os.getcwd()

	adp_file_info = adp_file.split('_')[0]

	print("WORKING FOLDER:\t%s%s" % (folder, os.path.sep))
	print("READING FILE:\t%s" % adp_file)

	tel = adp[0].header['TELESCOP'] # ESO-3P6 for HARPS, TNG for HARPS-N
	instr = adp[0].header['INSTRUME'] # instrument used

	print("TELESCOPE:\t%s" % tel)
	print("INSTRUMENT:\t%s" % instr)

	if instr == 'HARPS': obs = 'ESO'
	elif instr == 'HARPN': obs = 'TNG'
	else:
		print("*** ERROR: Instrument not recognized")

	try: obj = adp[0].header['OBJECT']
	except:
		try: obj = adp[0].header['%s OBS TARG NAME' % obs]
		except:
			print("*** ERROR: Cannot identify object")
			return

	print("OBJECT:\t\t%s" % obj)

	# Check if target is science or calibration file
	if obj in ('WAVE,WAVE,THAR1','WAVE,WAVE,THAR2'):
		print('*** ERROR: File is ThAr flux')
		return

	# Override object name with name given in obj_name option
	if obj_name:
		if type(obj_name) is list and len(obj_name) == 1:
			obj = obj_name[0]
		elif type(obj_name) is list and len(obj_name) > 1:
			print("*** ERROR: obj_name requires only one name, more than one given")
			return
		else: obj = obj_name

	# This is the observation date
	date = adp[0].header['DATE-OBS']

	adp.close()

	# Test CCF file
	try:
		filename = os.path.join(folder, "%s.%s*_ccf_*_A.fits" % (instr, date[:-2]))
		ccf_pfile = glob.glob(filename)[0]
		ccf_fits = pyfits.open(ccf_pfile)
		ccf_fits.close()
		print("CCF FILE:\t%s" % os.path.split(ccf_pfile)[-1])
	except:
		print("*** WARNING: Cannot read CCF file")
		ccf_pfile = None

	# Test BIS file
	try:
		filename = os.path.join(folder, "%s.%s*_bis_*_A.fits" % (instr, date[:-2]))
		bis_pfile = glob.glob(filename)[0]
		bis_fits = pyfits.open(bis_pfile)
		bis_fits.close()
		print("BIS FILE:\t%s" % os.path.split(bis_pfile)[-1])
	except:
		print("*** WARNING: Cannot read BIS file")
		bis_pfile = None

	files = {}
	files['adp'] = adp_pfile
	files['ccf'] = ccf_pfile
	files['bis'] = bis_pfile
	files['instr'] = instr
	files['obs'] = obs
	files['date'] = date
	files['obj'] = obj

	return files


def load_data_adp(files):
	"""
	Loads data from the files returned from read_file_adp.

		Functionalities:
		- Recognises if instrument is HARPS or HARPS-N.
		- Corrects wavelength for RV (rest frame).

	Parameters:
	-----------
	files : dict
		Dictionary with the file names required to compute the spectrum.

		The required keys are:

		==========  ========================================================
		keys		Description
		----------  --------------------------------------------------------
		adp 		str : ADP fits filename with path.
		ccf 		str : CCF fits filename with path.
		bis			{str, None} : BIS fits filename with path, None if not
					found.
		instr		str : Instrument identification.
		obs			str : Code related to instrument to be used in fits
					headers.
		==========  ========================================================

	Returns:
	--------
	data : dict
		Dictionary with data returned from fits files.

		Each key is a list with data related to a given measurement date
		given by the key 'date'.

		The returned keys are:

		==========  ========================================================
		keys		Description
		----------  --------------------------------------------------------
		wave 		list : Wavelength calibrated for BERV and RV (at rest
					frame) per pixel [angstroms].
		flux 		list : Deblazed flux per pixel.
		median_snr 	float : Median SNR of spectrum.
		obj 		str : Object identification.
		date 		str : Date of observation in the fits file format.
		bjd 		float : Barycentric Julian date of observation [days].
		rv			float : Radial velocity [m/s] (if CCF file available).
        rv_err		float : Error on radial velocity (photon noise) [m/s]
                    (if CCF file available).
        fwhm		float : Full-Width-at-Half-Maximum of the CCF line
                    profile [m/s] (if CCF file available).
        cont		float : Contrast of the CCF line profile [%] (if CCF
                    file available).
        bis			{float, None} : Bisector Inverse Span of the CCF line
                    [m/s] (if BIS file available).
        noise		float : CCF noise [m/s] (if CCF file available).
		instr		str : Instrument identification.
		data_flg 	None : Flag, not implemented.
		==========  ========================================================
	"""

	print("\nREADING DATA FROM ADP FITS FILE")
	print("-------------------------------")

	if not files:
		print("*** ERROR: files dictionary is empty")
		return

	flg = None

	obs = files['obs']

	adp = pyfits.open(files['adp'])

	flux = adp[1].data[0][1]
	print("Flux data read success")

	wave_orig = adp[1].data[0][0]
	print("Wave data read success")

	header = adp[0].header

	obj = files['obj']
	date = files['date']

	bjd = header['HIERARCH %s DRS BJD' % obs] # Barycentric Julian Day

	# median snr in all orders
	median_snr = header['SNR'] # Median SNR

	# SNR in all orders for HARPS
	if obs == 'ESO': orders = 72
	if obs == 'TNG': orders = 69
	snr = [float("%0.1f" % header['HIERARCH %s DRS SPE EXT SN%s' % (obs,k)]) for k in range(orders)]

	adp.close()

	if files['ccf']:
		ccf_fits = pyfits.open(files['ccf'])

		rv = ccf_fits[0].header['HIERARCH %s DRS CCF RVC' % obs] # [km/s], drift corrected
		rv_err = ccf_fits[0].header['HIERARCH %s DRS DVRMS' % obs] # [m/s]

		fwhm = ccf_fits[0].header['HIERARCH %s DRS CCF FWHM' % obs] # [km/s]
		cont = ccf_fits[0].header['HIERARCH %s DRS CCF CONTRAST' % obs] # [%]
		ccf_noise = ccf_fits[0].header['HIERARCH %s DRS CCF NOISE' % obs] # [km/s]

		ccf_fits.close()

		print("CCF data read success")

		# RV already corrected for BERV (from CCF file)
		rv = rv * 1000 # convert to m/s
		fwhm = fwhm * 1000 # convert to m/s
		ccf_noise = ccf_noise * 1000 # convert to m/s

		# Wavelength Doppler shift correction (BERV already corrected)
		delta_wave = rv * wave_orig / lightspeed
		wave = wave_orig - delta_wave
		print("Wavelength Doppler shift (to rest frame) corrected")

	elif not files['ccf']:
		rv = None; rv_err = None; fwhm = None; cont = None; ccf_noise = None
		print("*** ERROR: No CCF data available but required to calibrate wavelength")
		return

	if files['bis']:
		bis_fits = pyfits.open(files['bis'])
		bis = bis_fits[0].header['HIERARCH %s DRS BIS SPAN' % obs] # [km/s]
		bis = bis * 1000 # convert to m/s
		bis_fits.close()
		print("BIS data read success")
	elif not files['bis']:
		print("*** WARNING: No BIS data available")
		bis = None


	data = {}
	data['flux'] = flux
	data['wave'] = wave
	data['date'] = date
	data['obj'] = obj
	data['bjd'] = bjd
	data['median_snr'] = median_snr
	data['snr'] = snr
	data['rv'] = rv # [m/s]
	data['rv_err'] = rv_err # [m/s]
	data['fwhm'] = fwhm # [m/s]
	data['cont'] = cont
	data['bis'] = bis # [m/s]
	data['noise'] = ccf_noise # [m/s]
	data['instr'] = files['instr']
	data['data_flg'] = flg

	return data
