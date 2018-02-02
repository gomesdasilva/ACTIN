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
			return
	else:
		print "*** WARNING: No more %s files in folder" % file_type
		return


def calc_wave(e2ds_pfile, telescope):
	"""
	Compute wavelength from e2ds headers.

	Parameters:
	-----------
	e2ds_pfile : str
		e2ds file name with path to its location.
	telescope : str
		Telescope used, either 'ESO' for HARPS or, 'TNG' for HARPS-N.

	Returns:
	--------
	wave : list of lists
		2d wavelength list where len(wave) is number of orders [angstroms].
	"""

	e2ds = pyfits.open(e2ds_pfile)

	deg = e2ds[0].header['HIERARCH %s DRS CAL TH DEG LL' % telescope]

	ll_coeff = np.zeros((len(e2ds[0].data), deg + 1))

	# Read coefficients
	for i in range(len(e2ds[0].data)):
		for j in range(deg + 1):
			ll_coeff[i, j] = e2ds[0].header['HIERARCH {} DRS CAL TH COEFF '
                                    'LL{}'.format(telescope, (j + (deg + 1)*i))]

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
	(including wave, blaze and ccf). To be used with load_data_e2ds function.

		Some functionalities:
		- Automatically recognises if instrument is HARPS or HARPS-N.
		- Automatically searchs for wave and blaze files in the folder.
		- Automatically detects if e2ds file is science or calibration.

	Parameters:
	-----------
	e2ds_pfile : str
		e2ds fits file name with path.

	Returns:
	--------
	files : dict
		Dictionary with the file names required to compute the spectrum.

		The returned keys are:

		==========  ========================================================
		keys		Description
		----------  --------------------------------------------------------
		e2ds 		str : e2ds fits file name with path.
		wave 		{str, None} : Wave fits file name with path if found,
					None otherwise.
		blaze 		{str, None} : Blaze fits file name with path if found,
					None otherwise.
		ccf 		str : CCF fits file name with path.
		obj 		str : Object (target) identification.
		date 		str : Date of observation in the fits file format.
		==========  ========================================================

	"""

	print "\nREADING FILES FROM e2ds FITS FILE"
	print "--------------------------------"

	if e2ds_pfile == None:
		print "e2ds_file is 'None'."
		return

	if obj_name != None:
		if type(obj_name) is list and len(obj_name) == 1:
			obj_name = obj_name[0]

	files = {}

	folder = '/'.join(e2ds_pfile.split('/')[:-1])
	print "WORKING FOLDER:\t%s/" % folder

	e2ds_file = "/".join(e2ds_pfile.split("/")[-1:])
	print "READING FILE:\t%s" % e2ds_file

	e2ds_file_info = e2ds_file.split('_')[0]

	try:
		e2ds = pyfits.open(e2ds_pfile)
	except:
		print "*** ERROR: Cannot read %s" % e2ds_file
		return

	#date = e2ds[0].header['DATE-OBS'] # not the same as file date
	date = e2ds_file_info[6:] # This is the file date

	#instr = e2ds[0].header['INSTRUME'] ## only exists for HARPS dates > 2004

	tel = e2ds[0].header['TELESCOP'] # ESO-3P6 for HARPS, TNG for HARPS-N
	print "TELESCOPE:\t%s" % tel

	if obj_name == None:
		try: obj = e2ds[0].header['OBJECT']
		except:
			try: obj = e2ds[0].header['%s OBS TARG NAME' % tel[:3]]
			except:
				print "Cannot identify object"
				return
	else:
		obj = obj_name
	print "OBJECT:\t\t%s" % obj

	# Check if target is science or calibration file
	if obj in ('WAVE,WAVE,THAR1','WAVE,WAVE,THAR2'):
		print '*** ERROR: File is ThAr flux'
		return

	wave_file = e2ds[0].header['HIERARCH %s DRS CAL TH FILE' % tel[:3]]
	print "WAVE FILE:\t%s" % wave_file
	wave_pfile = "%s/%s" % (folder,wave_file)

	# NOT USING SERACH FILE
	# Test wave file, if not present search for others with a close date
	try:
		wave_fits = pyfits.open(wave_pfile)
		wave_fits.close()
	except:
		print "*** WARNING: The wave file associated with this e2ds is not present"
		#print "Looking for other wave files in the folder..."
		#wave_pfile = check_for_calib_files(e2ds[0].header,'wave',
		#									folder,dif_time_max=1.0)
		wave_pfile = None

	blaze_file = e2ds[0].header['HIERARCH %s DRS BLAZE FILE' % tel[:3]]
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

	ccf_pfile = glob.glob("%s/%s_ccf_*_A.fits" % (folder,e2ds_file_info))[0]
	print "CCF FILE:\t%s" % ccf_pfile.split('/')[-1]

	# Test ccf file
	try:
		ccf_fits = pyfits.open(ccf_pfile)
		ccf_fits.close()
	except:
		print "*** ERROR: Cannot read ccf file"
		return

	files['e2ds'] = e2ds_pfile
	files['wave'] = wave_pfile
	files['blaze'] = blaze_pfile
	files['ccf'] = ccf_pfile
	files['date'] = date # used for check_duplicate
	files['obj'] = obj # used for check_duplicate

	return files


def load_data_e2ds(files, obj_name=None):
	"""
	Loads data from the fits files returned from read_file_e2ds.

		Some functionalities:
		- Automatically recognises if instrument is HARPS or HARPS-N.
		- Automatically calculates wavelength using the e2ds coefficients if wave
		file is not given as input.
		- Normalizes flux if blaze file is given as input, otherwise ignores.
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
		flux 		list of floats : Flux per pixel per order.
		flux_deb 	{list of floats, None} : Flux per pixel per order
					normalized by the blaze function, None if CCF file name
					not given as an input parameter.
		wave 		list of floats : Wavelength calibrated for BERV and RV
					(at rest frame) per pixel per order [angstrom].
		snr 		list of floats : SNR at each spectral order.
		median_snr  float : Median SNR of spectrum.
		obj 		str : Object (target) identification.
		date 		str : Date of observation in the fits file format.
		bjd 		float : Barycentric Julian date of observation [days].
		rv 			float : Object radial velocity corrected for BERV [m/s].
		rv_err 		float : Error on the radial velocity (photon noise)
					[m/s].
		b_v 		float : B-V color of object as given in the 'e2ds' file.
		data_flg 	str : Flag with value 'noDeblazed' when the blaze file
					was not found (and flux_deb is real flux), None
					otherwise.
		==========  ========================================================

	"""

	print "\nREADING DATA FROM e2ds FITS FILE"
	print "--------------------------------"

	if files == None: return

	if obj_name != None:
		if type(obj_name) is list and len(obj_name) == 1:
			obj_name = obj_name[0]

	data = {}
	flg = None

	try:
		e2ds = pyfits.open(files['e2ds'])
		flux = np.asarray(e2ds[0].data)
		print "Flux data read success"
	except:
		print "*** ERROR: Cannot read e2ds file"
		return

	tel = e2ds[0].header['TELESCOP'] # ESO-3P6 for HARPS, TNG for HARPS-N

	try:
		wave_fits = pyfits.open(files['wave'])
		wave_orig = wave_fits[0].data
		wave_fits.close()
		print "Wave data read success"
	except:
		print "*** INFO: No wave file present, computing wave from e2ds"
		wave_orig = calc_wave(files['e2ds'],tel[:3]) ####
		wave_orig = np.asarray(wave_orig)

	try:
		blaze_fits = pyfits.open(files['blaze'])
		blaze = np.asarray(blaze_fits[0].data)

		# normalize blaze function
		blaze_norm = [blaze[k]/max(blaze[k]) for k in range(len(blaze))]
		blaze_norm = np.asarray(blaze_norm)

		blaze_fits.close()
		print "Blaze data read success"
	except:
		print "*** WARNING: No blaze file present, flux not deblazed"
		blaze_norm = np.ones([len(flux),len(flux[0])])
		blaze = np.ones([len(flux),len(flux[0])])
		#blaze = None
		flg = 'noDeblazed'

	try:
		ccf_fits = pyfits.open(files['ccf'])
		print "CCF data read success"
	except:
		print "*** ERROR: No ccf file present"
		return

	if obj_name == None:
		obj = e2ds[0].header['%s OBS TARG NAME' % (tel[:3])]
	else: obj = obj_name
	bjd = e2ds[0].header['HIERARCH %s DRS BJD' % (tel[:3])]

	#date = e2ds[0].header['DATE-OBS'] # not the same as file date
	e2ds_file = "/".join(files['e2ds'].split("/")[-1:])
	e2ds_file_info = e2ds_file.split('_')[0]
	date = e2ds_file_info[6:] # This is the file date

	snr = [float("%0.1f" % e2ds[0].header['HIERARCH %s DRS SPE EXT SN%s' % (tel[:3],k)]) for k in range(len(flux))] # SNR in all orders
	try: b_v = e2ds[0].header['HIERARCH %s DRS CAII B-V' % tel[:3]]
	except: b_v = None

	e2ds.close()

	rv = ccf_fits[0].header['HIERARCH %s DRS CCF RVC' % tel[:3]] # [km/s], drift corrected
	rv_err = ccf_fits[0].header['HIERARCH %s DRS DVRMS' % tel[:3]] # [m/s]
	berv = ccf_fits[0].header['HIERARCH %s DRS BERV' % tel[:3]] # [km/s]

	ccf_fits.close()

	# Median SNR of all orders
	median_snr = np.median(snr)

	# Deblazing flux
	#if blaze_norm != None:
	flux_deb = flux/blaze#_norm
		#print "Flux deblazed"
	#else: flux_deb = None

	# RV already corrected for BERV (from CCF file)
	rv = rv * 1000 # convert to m/s
	berv = berv * 1000 # convert to m/s

	# Malavolta:
	#naxis1 = e2ds[0].header['NAXIS1']
	#naxis2 = e2ds[0].header['NAXIS2']
	#wave_rest = wave_orig * ((1.+berv/lightspeed)/(1.+rv/lightspeed))
	#dwave_rest = np.zeros([naxis2,naxis1],dtype=np.double)
	#dwave_rest[:,1:] = wave_rest[:,1:]-wave_rest[:,:-1]
	#dwave_rest[:,0] = dwave_rest[:,1]
	#wave = wave_rest

	# Wavelength Doppler shift (to rest frame) and BERV correction
	delta_wave = (rv - berv) * wave_orig / lightspeed
	wave = wave_orig - delta_wave
	print "Wavelength Doppler shift (to rest frame) and BERV correction"

	data['flux'] = flux
	data['flux_deb'] = flux_deb
	data['wave'] = wave
	data['blaze'] = blaze #####
	data['snr'] = snr
	data['median_snr'] = median_snr
	data['obj'] = obj
	data['date'] = date
	data['bjd'] = bjd # [days]
	data['rv'] = rv # [m/s] ## not required by specha
	data['rv_err'] = rv_err # [m/s] ## not required by specha
	data['data_flg'] = flg
	data['b-v'] = b_v ## not required by specha

	return data
