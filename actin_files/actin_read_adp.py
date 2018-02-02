#!/usr/bin/env python

import glob
import numpy as np
import astropy.io.fits as pyfits

# Constants
lightspeed = 299792458.0 # [m/s]



def read_file_adp(adp_pfile):
	"""
	Reads ADP fits file and recovers all necessary files to compute
	spectrum (including CCF). To be used with load_data_adp function.

		Some functionalities:
		- Automatically recognises if instrument is HARPS or HARPS-N
		- Automatically detects if ADP file is science or calibration

	Parameters:
	-----------
	adp_pfile : str
		ADP fits file name with path.

	Returns:
	--------
	files : dict
		Dictionary with the file names required to compute the spectrum.

		The returned keys are:

		==========  ========================================================
		keys		Description
		----------  --------------------------------------------------------
		adp 		str : ADP fits file name with path.
		ccf 		str : CCF fits file name with path.
		obj 		str : Object (target) identification.
		date 		str : Date of observation in the fits file format.
		==========  ========================================================

	"""

	if adp_pfile == None: return

	print "\nREADING FILES FROM ADP FITS FILE"
	print "--------------------------------"

	files = {}

	folder = '/'.join(adp_pfile.split('/')[:-1])
	print "WORKING FOLDER:\t%s/" % folder

	adp_file = "/".join(adp_pfile.split("/")[-1:])
	print "READING FILE:\t%s" % adp_file

	try:
		adp = pyfits.open(adp_pfile)
	except:
		print "*** ERROR: Cannot read %s" % adp_file
		return

	date = adp[0].header['DATE-OBS']
	instr = adp[0].header['INSTRUME'] # instrument used

	tel = adp[0].header['TELESCOP'] # ESO-3P6 for HARPS, TNG for HARPS-N
	print "TELESCOPE:\t%s" % tel

	try: obj = adp[0].header['OBJECT']
	except:
		try: obj = adp[0].header['%s OBS TARG NAME' % tel[:3]]
		except: return
	print "OBJECT:\t\t%s" % obj

	if obj in ('WAVE,WAVE,THAR1','WAVE,WAVE,THAR2'):
		print '*** ERROR: File is ThAr flux'
		return

	# using date[:-1] because the time of the production of ccf file can be
	# some miliseconds after the date-obs
	try:
		ccf_pfile = glob.glob('%s/%s.%s*_ccf_*.fits' % (folder,instr,date[:-1]))[0]

		ccf_fits = pyfits.open(ccf_pfile)
		ccf_fits.close()
	except:
		print "*** ERROR: Cannot read ccf file"
		return


	# Test ccf file
	#try:
	#	ccf_fits = pyfits.open(ccf_pfile)
	#	ccf_fits.close()
	#except:
	#	print "*** ERROR: Cannot read ccf file"
	#	return

	files['adp'] = adp_pfile
	files['ccf'] = ccf_pfile
	files['date'] = date # for check_duplicate
	files['obj'] = obj # for check_duplicate

	return files


def load_data_adp(files):
	"""
	Loads data from the files returned from read_file_adp.

		Some functionalities:
		- Automatically recognises if instrument is HARPS or HARPS-N.
		- Corrects wavelength for RV (rest frame).

	Parameters:
	-----------
	files : dict
		Dictionary with the file names required to compute the spectrum.

		The required keys are:

		==========  ========================================================
		keys		Description
		----------  --------------------------------------------------------
		adp 		str : ADP fits file name with path.
		ccf 		str : CCF fits file name with path.
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
		wave 		float : Wavelength calibrated for BERV and RV (at rest
					frame) per pixel [angstrom].
		flux 		float : Flux per pixel.
		median_snr 	float : Median SNR of spectrum.
		obj 		str : Object (target) identification.
		date 		str : Date of observation in the fits file format.
		bjd 		float : Barycentric Julian date of observation [days].
		rv 			float : Object radial velocity corrected for BERV [m/s].
		rv_err 		float : Error on the radial velocity (photon noise)
					[m/s].
		b_v 		float, None : B-V color of object as given in the ADP
					file, None if not available.
		data_flg 	None : Flag, not implemented.
		==========  ========================================================

	"""

	if files == None: return

	print "\nREADING DATA FROM ADP FITS FILE"
	print "-------------------------------"

	data = {}
	flg = None

	try:
		adp = pyfits.open(files['adp'])
	except:
		print "ERROR: Cannot read adp file"
		return

	flux = adp[1].data[0][1]
	print "Flux data read success"

	wave_orig = adp[1].data[0][0]
	print "Wave data read success"

	hdr = adp[0].header

	tel = hdr['TELESCOP'] # ESO-3P6 for HARPS, TNG for HARPS-N
	date = hdr['DATE-OBS'] # Observation date in file format
	obj = hdr['%s OBS TARG NAME' % (tel[:3])] # Target Id
	median_snr = hdr['SNR'] # Median SNR
	bjd = hdr['HIERARCH %s DRS BJD' % (tel[:3])] # Barycentric Julian Day
	try: b_v = hdr[0].header['HIERARCH %s DRS CAII B-V' % tel[:3]]
	except: b_v = None

	# SNR in all orders for HARPS
	if tel[:3] == 'ESO': orders = 71
	if tel[:3] == 'TNG': orders = 68
	snr = ["%0.1f" % hdr['HIERARCH %s DRS SPE EXT SN%s' % (tel[:3],orders)]]

	# convert str to float to calculate median
	median_snr_calc = np.median([float(snr[k]) for k in range(len(snr))]) # not used

	#try:
	ccf_fits = pyfits.open(files['ccf'])
		#print "CCF file read success"
	#except:
		#print "*** ERROR: No ccf file present"
		#return

	rv = ccf_fits[0].header['HIERARCH %s DRS CCF RVC' % tel[:3]] # [km/s], drift corrected
	rv_err = ccf_fits[0].header['HIERARCH %s DRS DVRMS' % tel[:3]] # [m/s]
	berv = ccf_fits[0].header['HIERARCH %s DRS BERV' % tel[:3]] # [km/s]
	ccf_fits.close()

	# RV already corrected for BERV (from CCF file)
	rv = rv * 1000 # convert to m/s
	berv = berv * 1000 # convert to m/s

	# Wavelength Doppler shift (to rest frame) correction (BERV already corrected)
	delta_wave = rv * wave_orig / lightspeed
	wave = wave_orig - delta_wave
	print "Wavelength Doppler shift (to rest frame) corrected"

	data['flux'] = flux
	data['wave'] = wave
	data['date'] = date
	data['obj'] = obj
	data['bjd'] = bjd
	data['median_snr'] = median_snr
	data['snr'] = snr
	data['rv'] = rv
	data['rv_err'] = rv_err
	data['data_flg'] = flg
	data['b-v'] = b_v

	return data
