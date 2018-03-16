#!/usr/bin/env python

import glob
import numpy as np
import astropy.io.fits as pyfits

# Constants
lightspeed = 299792458.0 # [m/s]


def read_file_s1d(s1d_pfile, obj_name=None):
	"""
	Reads s1d and s1d_*_rv fits files and recovers all necessary files
	to compute spectrum (including ccf).

		Functionalities:
		- Recognises if instrument is HARPS or HARPS-N.
		- Ignores file if not science (ThAr flux).

	Parameters:
	-----------
	s1d_pfile : str
		s1d or s1d_*_rv fits filename with path.
	obj_name : str (optional)
		Name given to object that overrides the name from fits file. None is default.

	Returns:
	--------
	files : dict
		Dictionary with the file names required to compute the spectrum, object name and date.

		The returned keys are:

		==========  ========================================================
		keys		Description
		----------  --------------------------------------------------------
		s1d 		str : s1d fits filename with path.
		ccf 		{str, None} : CCF fits filename with path, None if not
					found.
		instr		str : Instrument identification.
		obs			str : Code related to instrument to be used in fits
					headers.
		obj 		str : Object identification.
		date 		str : Date of observation in the fits file format.
		==========  ========================================================
	"""

	print "\nREADING FILES FROM s1d FITS FILE"
	print "--------------------------------"

	if s1d_pfile is None:
		print "*** ERROR: s1d file with path required"
		return

	try: s1d = pyfits.open(s1d_pfile)
	except:
		print "*** ERROR: Cannot read %s" % s1d_pfile
		return

	folder = '/'.join(s1d_pfile.split('/')[:-1])
	s1d_file = "/".join(s1d_pfile.split("/")[-1:])
	s1d_file_info = s1d_file.split('_')[0]

	print "WORKING FOLDER:\t%s/" % folder
	print "READING FILE:\t%s" % s1d_file

	tel = s1d[0].header['TELESCOP'] # ESO-3P6 for HARPS, TNG for HARPS-N
	instr = s1d[0].header['INSTRUME'] # instrument used

	print "TELESCOPE:\t%s" % tel[:3]
	print "INSTRUMENT:\t%s" % instr

	if instr == 'HARPS': obs = 'ESO'
	elif instr == 'HARPN': obs = 'TNG'
	else:
		print "*** ERROR: Instrument not recognized"

	try: obj = s1d[0].header['OBJECT']
	except:
		try: obj = s1d[0].header['%s OBS TARG NAME' % obs]
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

	date = s1d_file_info[6:]

	if not folder:
		ccf_pfile = glob.glob('%s_ccf_*.fits' % (s1d_file_info))[0]
	if folder:
		ccf_pfile = glob.glob('%s/%s_ccf_*.fits' % (folder,s1d_file_info))[0]

	print "CCF FILE:\t%s" % ccf_pfile.split('/')[-1]

	# Test ccf file
	try:
		ccf_fits = pyfits.open(ccf_pfile)
		ccf_fits.close()
	except:
		print "*** ERROR: Cannot read ccf file"
		ccf_pfile = None

	files = {}
	files['s1d'] = s1d_pfile
	files['ccf'] = ccf_pfile
	files['instr'] = instr
	files['obs'] = obs
	files['date'] = date
	files['obj'] = obj

	return files


def load_data_s1d(files):
	"""
	Loads data from the fits files returned from read_file_s1d.

		Functionalities:
		- Recognises if instrument is HARPS or HARPS-N.
		- Calculates wavelength using the s1d coefficients.
		- Corrects wavelength for RV (rest frame) if file is s1d (s1d_rv
		  have wavelength already at rest frame).

	Parameters:
	-----------
	files : dict
		Dictionary with the file names required to compute the spectrum.

		The used keys are:

		==========  ========================================================
		keys		Description
		----------  --------------------------------------------------------
		s1d 		str : s1d fits file name with path.
		ccf 		{str, None} : CCF fits file name with path, None if not
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
		flux 		list : Deblazed flux per pixel.
		wave 		list : Wavelength calibrated for BERV and RV (at rest
					frame) per pixel [angstroms].
		snr 		list : SNR at each spectral order.
		median_snr  float : Median SNR of spectrum.
		obj 		str : Object identification.
		date 		str : Date of observation in the fits file format.
		bjd 		float : Barycentric Julian date of observation [days].
		rv 			float : Radial velocity corrected for BERV [m/s].
		rv_err 		float : Error on radial velocity [m/s].
		instr		str : Instrument identification.
		data_flg 	None : Flag, not implemented.
		==========  ========================================================
	"""

	print "\nREADING DATA FROM s1d FITS FILE"
	print "-------------------------------"

	if files is None:
		print "*** ERROR: files dictionary is empty"
		return

	flg = None

	obs = files['obs']

	s1d = pyfits.open(files['s1d'])

	flux = s1d[0].data
	print "Flux data read success"

	header = s1d[0].header
	wave_orig = header['CRVAL1'] + header['CDELT1']*np.arange(header['NAXIS1'])

	obj = files['obj']
	date = files['date']

	bjd = header['HIERARCH %s DRS BJD' % obs]

	# SNR in all orders for HARPS/N
	if obs == 'ESO': orders = 72
	if obs == 'TNG': orders = 69
	snr = [float("%0.1f" % header['HIERARCH %s DRS SPE EXT SN%s' % (obs,k)]) for k in range(orders)]

	s1d.close()

	# Median snr in all orders
	median_snr = np.median(snr)

	# To know if file is s1d or s1d_A_rv
	rest_frame = files['s1d'].split('_')[-1].split('.')[0]

	# If using s1d_*_rv the wavelength is already at rest frame
	if files['ccf'] is None and rest_frame == 'rv':
		rv, rv_err = None
		wave = wave_orig

	# If using s1d_*_rv with ccf file
	elif files['ccf'] and rest_frame == 'rv':
		ccf_fits = pyfits.open(files['ccf'])
		print "CCF data read success"

		rv = ccf_fits[0].header['HIERARCH %s DRS CCF RVC' % obs] # [km/s], drift corrected
		rv_err = ccf_fits[0].header['HIERARCH %s DRS DVRMS' % obs] # [m/s]

		ccf_fits.close()

		# RV already corrected for BERV (from CCF file)
		rv = rv * 1000 # convert to m/s

		wave = wave_orig

	# If using s1d need rv to calibrate wavelength
	elif files['ccf'] is None and rest_frame != 'rv':
		print "*** ERROR: No ccf file found but required to calibrate wavelength"
		return

	# If using s1d use rv to calibrate wavelength
	elif files['ccf'] and rest_frame != 'rv':
		ccf_fits = pyfits.open(files['ccf'])
		print "CCF data read success"

		rv = ccf_fits[0].header['HIERARCH %s DRS CCF RVC' % obs] # [km/s], drift corrected
		rv_err = ccf_fits[0].header['HIERARCH %s DRS DVRMS' % obs] # [m/s]

		ccf_fits.close()

		# RV already corrected for BERV (from CCF file)
		rv = rv * 1000 # convert to m/s

		wave = wave_orig - rv * wave_orig / lightspeed

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
