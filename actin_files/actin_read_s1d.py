#!/usr/bin/env python

import glob
import numpy as np
import astropy.io.fits as pyfits

# Constants
lightspeed = 299792458.0 # [m/s]




def read_file_s1d(s1d_pfile):
	"""
	Reads s1d and s1d_*_rv fits files and recovers all necessary files
	to compute spectrum (including ccf). To be used with load_data_s1d
	function.

		Some functionalities:
		- Automatically recognises if instrument is HARPS or HARPS-N
		- Automatically detects if s1d file is science or calibration

	Parameters:
	-----------
	s1d_pfile : str
		s1d or s1d_*_rv fits file name with path.

	Returns:
	--------
	files : dict
		Dictionary with the file names required to compute the spectrum.

		The returned keys are:

		==========  ========================================================
		keys		Description
		----------  --------------------------------------------------------
		s1d 		str : s1d fits file name with path.
		ccf 		{str, None} : CCF fits file name with path, None if not
					found.
		obj 		str : Object (target) identification.
		date 		str : Date of observation in the fits file format.
		==========  ========================================================

	"""

	if s1d_pfile == None: return

	print "\nREADING FILES FROM s1d FITS FILE"
	print "--------------------------------"

	files = {}

	folder = '/'.join(s1d_pfile.split('/')[:-1])
	print "WORKING FOLDER:\t%s/" % folder

	s1d_file = "/".join(s1d_pfile.split("/")[-1:])
	print "READING FILE:\t%s" % s1d_file

	s1d_file_info = s1d_file.split('_')[0]

	try:
		s1d = pyfits.open(s1d_pfile)
	except:
		print "*** ERROR: Cannot read %s" % s1d_file
		return

	date = s1d[0].header['DATE-OBS']

	tel = s1d[0].header['TELESCOP'] # ESO-3P6 for HARPS, TNG for HARPS-N
	print "TELESCOPE:\t%s" % tel[:3]

	try: obj = s1d[0].header['OBJECT']
	except:
		try: obj = s1d[0].header['%s OBS TARG NAME' % tel[:3]]
		except: return
	print "OBJECT:\t%s" % obj

	if obj in ('WAVE,WAVE,THAR1','WAVE,WAVE,THAR2'):
		print '*** ERROR: File is ThAr flux'
		return

	ccf_pfile = glob.glob('%s/%s_ccf_*.fits' % (folder,s1d_file_info))[0]
	print "CCF FILE:\t%s" % ccf_pfile.split('/')[-1]
	# Test ccf file
	try:
		ccf_fits = pyfits.open(ccf_pfile)
		ccf_fits.close()
	except:
		print "*** ERROR: Cannot read ccf file"
		ccf_pfile = None

	files['s1d'] = s1d_pfile
	files['ccf'] = ccf_pfile
	files['date'] = date # for check_duplicate
	files['obj'] = obj # for check_duplicate

	return files


def load_data_s1d(files):
	"""
	Loads data from the fits files returned from read_file_s1d.

		Some functionalities:
		- Automatically recognises if instrument is HARPS or HARPS-N.
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
		flux 		float : Normalized flux per pixel.
		wave 		float : Wavelength calibrated for BERV and RV (at rest
					frame) per pixel [angstrom].
		snr 		list of floats : SNR at each spectral order.
		median_snr  float : Median SNR of spectrum.
		obj 		str : Object (target) identification.
		date 		str : Date of observation in the fits file format.
		bjd 		float : Barycentric Julian date of observation [days].
		rv 			float : Object radial velocity corrected for BERV [m/s].
		rv_err 		float : Error on the radial velocity (photon noise)
					[m/s].
		b_v 		float : B-V color of object as given in the 'e2ds' file.
		data_flg 	None : Flag, not implemented.
		==========  ========================================================

	"""

	if files == None: return

	print "\nREADING DATA FROM s1d FITS FILE"
	print "-------------------------------"

	data = {}
	flg = None

	try:
		s1d = pyfits.open(files['s1d'])
		flux = s1d[0].data
		print "Flux data read success"
	except:
		print "ERROR: Cannot read s1d file"
		return

	hdr = s1d[0].header

	# Calculate wavelength from s1d header (BERV corrected)
	wave_ref_pixel = hdr['CRVAL1'] # dispersion at ref. pixel (wavelength)
	ref_pixel = hdr['CRPIX1'] # reference pixel
	wave_step = hdr['CDELT1'] # dispersion by pixel (wavelength step)
	wave_orig = []
	for k in range(len(flux)):
		wave_ref_pixel += wave_step
		wave_orig.append(wave_ref_pixel)
	print "Wavelength calculated"

	tel = hdr['TELESCOP'] # ESO-3P6 for HARPS, TNG for HARPS-N
	date = hdr['DATE-OBS'] # Observation date in file format
	obj = hdr['%s OBS TARG NAME' % (tel[:3])] # Target Id
	bjd = hdr['HIERARCH %s DRS BJD' % (tel[:3])] # Barycentric Julian Day
	try: b_v = hdr['HIERARCH %s DRS CAII B-V' % tel[:3]]
	except: b_v = None

	# SNR in all orders for HARPS/N
	if tel[:3] == 'ESO': orders = 72
	if tel[:3] == 'TNG': orders = 69
	snr = [float("%0.2f" % hdr['HIERARCH %s DRS SPE EXT SN%s' % (tel[:3],k)]) for k in range(orders)]

	# median snr in all orders
	median_snr = np.median(snr)

	# to know if file is s1d or s1d_A_rv
	rest_frame = files['s1d'].split('_')[-1].split('.')[0]

	# Test ccf file
	if files['ccf'] == None and rest_frame == 'rv':
		rv, rv_err = None
	elif files['ccf'] == None and rest_frame != 'rv':
		print "*** ERROR: No ccf file found but required to calibrate wavelength"
		return
	else:
		ccf_fits = pyfits.open(files['ccf'])
		print "CCF data read success"

		rv = ccf_fits[0].header['HIERARCH %s DRS CCF RVC' % tel[:3]] # [km/s], drift corrected
		rv_err = ccf_fits[0].header['HIERARCH %s DRS DVRMS' % tel[:3]] # [m/s]
		ccf_fits.close()

		# RV already corrected for BERV (from CCF file)
		rv = rv * 1000 # convert to m/s

	# to know if file is s1d or s1d_*_rv
	rest_frame = files['s1d'].split('_')[-1].split('.')[0]

	if rest_frame == 'rv':
		wave = wave_orig # s1d_*_rv already at rest frame
	else:
		wave = wave_orig - rv * np.asarray(wave_orig) / lightspeed # s1d not at rest frame

	data['flux'] = flux
	data['wave'] = wave
	data['date'] = date
	data['obj'] = obj
	data['bjd'] = bjd
	data['median_snr'] = median_snr
	data['snr'] = snr
	data['rv'] = rv
	data['rv_err'] = rv_err
	data['b-v'] = b_v
	data['data_flg'] = flg

	return data
