#!/usr/bin/env python

import glob
import numpy as np
import astropy.io.fits as pyfits

# Constants
lightspeed = 299792458.0 # [m/s]


def check_for_calib_files(e2ds_header, file_type, folder, dif_time_max=1.0):
	"""
	Check for calibration files (wave or blaze) in the working directory
	in case of predetermined files are not present and choses the ones with
	less than 1 day difference from the original e2ds file.

	Parameters:
	-----------
	e2ds_header : fits header
		Fits header of e2ds fits file.
	file_type : str
		File type to be searched, either 'wave' or 'blaze'.
	folder : str
		Path to the fits files.
	dif_time_max : float
		Maximum difference in time between e2ds and
		the chosen calib file [days].

	Returns:
	--------
	calib_pfile : {str, None}
		Selected calibration file with path included, None if
		file not found.
	"""

	mjd_e2ds = e2ds_header['MJD-OBS']
	calib_pfiles = glob.glob('%s/*%s_A.fits' % (folder,file_type))
	if len(calib_pfiles) > 0:
		time_dif = []
		for k in range(len(calib_pfiles)):
			w = pyfits.open(calib_pfiles[k])
			mjd = w[0].header['MJD-OBS']
			w.close()
			time_dif.append(abs(mjd_e2ds-mjd))
		if min(time_dif) < dif_time_max:
			index = time_dif.index(min(time_dif))
			calib_pfile = calib_pfiles[index]
			calib = pyfits.open(calib_pfile) # test new file
			print "New %s file used:\t%s" % (file_type,calib_pfile.split('/')[-1])
			print "%s file time difference to e2ds = %.2f days" % (file_type,min(time_dif))
			calib.close()
			return calib_pfile
		else:
			print "*** WARNING: Closest %s file was produced longer than 1 day" % file_type
			return None
	else:
		print "*** WARNING: No more %s files in folder" % file_type
		return None


def calc_wave(e2ds_pfile, obs):
	"""
	Compute wavelength from e2ds headers.

	Parameters:
	-----------
	e2ds_pfile : str
		e2ds file name with path to its location.
	obs : str
		Code related to instrument to be used in fits headers.

	Returns:
	--------
	wave : list of lists
		2d wavelength list where len(wave) is number of orders [angstroms].
	"""

	e2ds = pyfits.open(e2ds_pfile)

	deg = e2ds[0].header['HIERARCH %s DRS CAL TH DEG LL' % obs]

	ll_coeff = np.zeros((len(e2ds[0].data), deg + 1))

	# Read coefficients
	for i in range(len(e2ds[0].data)):
		for j in range(deg + 1):
			ll_coeff[i, j] = e2ds[0].header['HIERARCH {} DRS CAL TH COEFF '
                                    'LL{}'.format(obs, (j + (deg + 1)*i))]

	# Evaluate polynomials
	x = np.arange(e2ds[0].data.shape[1])  # Pixel array
	wave = np.zeros(e2ds[0].data.shape)  # Wavelength 2D array
	for i in range(len(wave)):
		wave[i] = np.poly1d(ll_coeff[i][::-1])(x)

	e2ds.close()

	return wave


def read_file_e2ds(e2ds_pfile, obj_name=None):
	"""
	Reads e2ds fits file and recovers all necessary files to compute spectrum
	(including wave, blaze and ccf).

		Functionalities:
		- Recognises if instrument is HARPS or HARPS-N.
		- Calculates wavelength if wave file is not present.
		- Searchs for blaze files in the folder if designated file not found.
		- Ignores file if not science (ThAr flux).

	Parameters:
	-----------
	e2ds_pfile : str
		e2ds fits filename with path.
	obj_name : str (optional)
		Name given to object that overrides the name from fits file. None is default.

	Returns:
	--------
	files : dict
		Dictionary with the filenames required to compute the spectrum, object name and date.

		The returned keys are:

		==========  ========================================================
		keys		Description
		----------  --------------------------------------------------------
		e2ds 		str : e2ds fits filename with path.
		wave 		{str, None} : Wave fits filename with path if found,
					None otherwise.
		blaze 		{str, None} : Blaze fits filename with path if found,
					None otherwise.
		ccf 		str : CCF fits filename with path.
		instr		str : Instrument identification.
		obs			str : Code related to instrument to be used in fits
					headers.
		obj 		str : Object identification.
		date 		str : Date of observation in the fits file format.
		==========  ========================================================
	"""

	print "\nREADING FILES FROM e2ds FITS FILE"
	print "--------------------------------"

	if e2ds_pfile is None:
		print "*** ERROR: e2ds file with path required"
		return

	try: e2ds = pyfits.open(e2ds_pfile)
	except:
		print "*** ERROR: Cannot read %s" % e2ds_pfile
		return

	folder = '/'.join(e2ds_pfile.split('/')[:-1])

	e2ds_file = '/'.join(e2ds_pfile.split('/')[-1:])
	e2ds_file_info = e2ds_file.split('_')[0]

	print "WORKING FOLDER:\t%s/" % folder
	print "READING FILE:\t%s" % e2ds_file

	tel = e2ds[0].header['TELESCOP'] # ESO-3P6 for HARPS, TNG for HARPS-N
	instr = e2ds[0].header['INSTRUME'] # instrument used

	print "TELESCOPE:\t%s" % tel
	print "INSTRUMENT:\t%s" % instr

	if instr == 'HARPS': obs = 'ESO'
	elif instr == 'HARPN': obs = 'TNG'
	else:
		print "*** ERROR: Instrument not recognized"

	try: obj = e2ds[0].header['OBJECT']
	except:
		try: obj = e2ds[0].header['%s OBS TARG NAME' % obs]
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

	date = e2ds_file_info[6:] # This is the file date

	wave_file = e2ds[0].header['HIERARCH %s DRS CAL TH FILE' % obs]

	print "WAVE FILE:\t%s" % wave_file

	wave_pfile = "%s/%s" % (folder, wave_file)

	try:
		wave_fits = pyfits.open(wave_pfile)
		wave_fits.close()
	except:
		print "*** WARNING: The wave file associated with this e2ds is not present"
		wave_pfile = None

	blaze_file = e2ds[0].header['HIERARCH %s DRS BLAZE FILE' % obs]

	e2ds.close()

	print "BLAZE FILE:\t%s" % blaze_file

	blaze_pfile = '%s/%s' % (folder,blaze_file)

	# Test blaze file, if not present searchs for other with a close date
	try:
		blaze_fits = pyfits.open(blaze_pfile)
		blaze_fits.close()
	except:
		print "*** WARNING: The blaze file associated with this e2ds is not present"
		print "Looking for other blaze files in the folder..."
		blaze_pfile = check_for_calib_files(e2ds[0].header,'blaze',folder,dif_time_max=1.0)

	if not folder:
		ccf_pfile = glob.glob("%s_ccf_*_A.fits" % (e2ds_file_info))[0]
	if folder:
		ccf_pfile = glob.glob("%s/%s_ccf_*_A.fits" % (folder,e2ds_file_info))[0]

	print "CCF FILE:\t%s" % ccf_pfile.split('/')[-1]

	# Test ccf file
	try:
		ccf_fits = pyfits.open(ccf_pfile)
		ccf_fits.close()
	except:
		print "*** ERROR: Cannot read ccf file"
		return

	files = {}
	files['e2ds'] = e2ds_pfile
	files['wave'] = wave_pfile
	files['blaze'] = blaze_pfile
	files['ccf'] = ccf_pfile
	files['instr'] = instr
	files['obs'] = obs
	files['date'] = date
	files['obj'] = obj

	return files


def load_data_e2ds(files):
	"""
	Loads data from the fits files returned from read_file_e2ds.

		Functionalities:
		- Recognises if instrument is HARPS or HARPS-N.
		- Calculates wavelength using the e2ds coefficients if wave
		  file is not given as input.
		- Corrects wavelength for BERV and RV (rest frame).

	Parameters:
	-----------
	files : dict
		Dictionary with the file names required to compute the spectrum.

		The used keys are:

		==========  ========================================================
		keys		Description
		----------  --------------------------------------------------------
		e2ds 		str : e2ds fits file name with path.
		wave 		{str, None} : Wave fits file name with path, None if not
					present.
		blaze 		{str, None} : Blaze fits file name with path, None if
					not present.
		ccf 		str : CCF fits file name with path.
		instr		str : Instrument identification.
		obs			str : Code related to instrument to be used in fits
					headers.
		obj			str : Object identification.
		date		str : Date of observation in the fits file format.
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
		flux 		list of lists : Flux per pixel per order.
		wave 		list of lists : Wavelength calibrated for BERV and RV
					(at rest frame) per pixel per order [angstroms].
		blaze		list of lists : Blaze function.
		snr 		list : SNR at each spectral order.
		median_snr  float : Median SNR of spectrum.
		obj 		str : Object identification.
		date 		str : Date of observation in the fits file format.
		bjd 		float : Barycentric Julian date of observation [days].
		rv 			float : Radial velocity corrected for BERV [m/s].
		rv_err 		float : Error on radial velocity [m/s].
		instr		str : Instrument identification.
		data_flg 	str : Flag with value 'noDeblazed' when the blaze file
					is not found, None otherwise.
		==========  ========================================================
	"""

	print "\nREADING DATA FROM e2ds FITS FILE"
	print "--------------------------------"

	if files is None:
		print "*** ERROR: files dictionary is empty."
		return

	flg = None

	obs = files['obs']

	e2ds = pyfits.open(files['e2ds'])

	flux = np.asarray(e2ds[0].data)
	print "Flux data read success"

	if files['wave'] is not None:
		wave_fits = pyfits.open(files['wave'])
		wave_orig = wave_fits[0].data
		wave_fits.close()
		print "Wave data read success"
	elif files['wave'] is None:
		print "*** INFO: No wave file present, computing wave from e2ds"
		wave_orig = calc_wave(files['e2ds'],obs)
		wave_orig = np.asarray(wave_orig)

	if files['blaze'] is not None:
		blaze_fits = pyfits.open(files['blaze'])
		blaze = np.asarray(blaze_fits[0].data)
		blaze_fits.close()
		print "Blaze data read success"
	elif files['blaze'] is None:
		print "*** WARNING: No blaze file present, flux not deblazed"
		blaze_norm = np.ones([len(flux),len(flux[0])])
		blaze = np.ones([len(flux),len(flux[0])])
		flg = 'noDeblazed'

	ccf_fits = pyfits.open(files['ccf'])
	print "CCF data read success"

	obj = files['obj']
	date = files['date'] # date as in the fits filename

	bjd = e2ds[0].header['HIERARCH %s DRS BJD' % (obs)]

	e2ds_file = "/".join(files['e2ds'].split("/")[-1:])
	e2ds_file_info = e2ds_file.split('_')[0]

	snr = [float("%0.1f" % e2ds[0].header['HIERARCH %s DRS SPE EXT SN%s' % (obs,k)]) for k in range(len(flux))] # SNR in all orders

	e2ds.close()

	rv = ccf_fits[0].header['HIERARCH %s DRS CCF RVC' % obs] # [km/s], drift corrected
	rv_err = ccf_fits[0].header['HIERARCH %s DRS DVRMS' % obs] # [m/s]
	berv = ccf_fits[0].header['HIERARCH %s DRS BERV' % obs] # [km/s]

	ccf_fits.close()

	# Median SNR of all orders
	median_snr = np.median(snr)

	# RV already corrected for BERV (from CCF file)
	rv = rv * 1000 # convert to m/s
	berv = berv * 1000 # convert to m/s

	# Wavelength Doppler shift (to rest frame) and BERV correction
	delta_wave = (rv - berv) * wave_orig / lightspeed
	wave = wave_orig - delta_wave
	print "Wavelength corrected for RV and BERV (at rest frame)"

	data = {}
	data['flux'] = flux
	data['wave'] = wave
	data['blaze'] = blaze
	data['snr'] = snr
	data['median_snr'] = median_snr
	data['obj'] = obj
	data['date'] = date
	data['bjd'] = bjd # [days]
	data['rv'] = rv # [m/s]
	data['rv_err'] = rv_err # [m/s]
	data['instr'] = files['instr']
	data['data_flg'] = flg

	return data