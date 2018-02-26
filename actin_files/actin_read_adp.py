#!/usr/bin/env python

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
		instr		str : Instrument identification.
		obs			str : Code related to instrument to be used in fits
					headers.
		obj 		str : Object identification.
		date 		str : Date of observation in the fits file format.
		==========  ========================================================
	"""

	print "\nREADING FILES FROM ADP FITS FILE"
	print "--------------------------------"

	if adp_pfile is None:
		print "*** ERROR: ADP file with path required"
		return

	try: adp = pyfits.open(adp_pfile)
	except:
		print "*** ERROR: Cannot read %s" % adp_pfile
		return

	folder = '/'.join(adp_pfile.split('/')[:-1])
	adp_file = "/".join(adp_pfile.split("/")[-1:])

	print "WORKING FOLDER:\t%s/" % folder
	print "READING FILE:\t%s" % adp_file

	tel = adp[0].header['TELESCOP'] # ESO-3P6 for HARPS, TNG for HARPS-N
	instr = adp[0].header['INSTRUME'] # instrument used

	print "TELESCOPE:\t%s" % tel
	print "INSTRUMENT:\t%s" % instr

	if instr == 'HARPS': obs = 'ESO'
	elif instr == 'HARPN': obs = 'TNG'
	else:
		print "*** ERROR: Instrument not recognized"

	try: obj = adp[0].header['OBJECT']
	except:
		try: obj = adp[0].header['%s OBS TARG NAME' % obs]
		except:
			print "*** ERROR: Cannot identify object"
			return

	print "OBJECT:\t\t%s" % obj

	# Check if target is science or calibration file
	if obj in ('WAVE,WAVE,THAR1','WAVE,WAVE,THAR2'):
		print '*** ERROR: File is ThAr flux'
		return

	# Override object name with name given in obj_name option
	if obj_name is not None:
		if type(obj_name) is list and len(obj_name) == 1:
			obj = obj_name[0]
		elif type(obj_name) is list and len(obj_name) > 1:
			print "*** ERROR: obj_name requires only one name, more than one given"
			return
		else: obj = obj_name

	date = adp[0].header['DATE-OBS']

	adp.close()

	ccf_pfile = glob.glob('%s/%s.%s*_ccf_*.fits' % (folder,instr,date[:-2]))[0]

	print "CCF FILE:\t%s" % ccf_pfile.split('/')[-1]

	# Test CCF file
	try:
		ccf_fits = pyfits.open(ccf_pfile)
		ccf_fits.close()
	except:
		print "*** ERROR: Cannot read ccf file, but required for wavelength calibration"
		return

	files = {}
	files['adp'] = adp_pfile
	files['ccf'] = ccf_pfile
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
		rv 			float : Radial velocity corrected for BERV [m/s].
		rv_err 		float : Error on the radial velocity [m/s].
		instr		str : Instrument identification.
		data_flg 	None : Flag, not implemented.
		==========  ========================================================
	"""

	print "\nREADING DATA FROM ADP FITS FILE"
	print "-------------------------------"

	if files is None:
		print "*** ERROR: files dictionary is empty"
		return

	flg = None

	obs = files['obs']

	adp = pyfits.open(files['adp'])

	flux = adp[1].data[0][1]
	print "Flux data read success"

	wave_orig = adp[1].data[0][0]
	print "Wave data read success"

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

	ccf_fits = pyfits.open(files['ccf'])

	rv = ccf_fits[0].header['HIERARCH %s DRS CCF RVC' % obs] # [km/s], drift corrected
	rv_err = ccf_fits[0].header['HIERARCH %s DRS DVRMS' % obs] # [m/s]

	ccf_fits.close()

	# RV already corrected for BERV (from CCF file)
	rv = rv * 1000 # convert to m/s

	# Wavelength Doppler shift correction (BERV already corrected)
	delta_wave = rv * wave_orig / lightspeed
	wave = wave_orig - delta_wave
	print "Wavelength Doppler shift (to rest frame) corrected"

	data = {}
	data['flux'] = flux
	data['wave'] = wave
	data['date'] = date
	data['obj'] = obj
	data['bjd'] = bjd
	data['median_snr'] = median_snr
	data['snr'] = snr
	data['rv'] = rv
	data['rv_err'] = rv_err
	data['instr'] = files['instr']
	data['data_flg'] = flg

	return data
