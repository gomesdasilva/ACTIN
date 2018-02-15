#!/usr/bin/env python

import sys,os
import glob
import time
import numpy as np
import astropy.io.fits as pyfits

path = os.path.dirname(os.path.realpath(__file__)) # directory of actin.py

# location of ACTIN files:
sys.path.append("%s/actin_files/" % path)
import actin_config
import actin_read_e2ds as read_e2ds
import actin_read_s1d as read_s1d
import actin_read_adp as read_adp
import actin_get_win as get_win
import actin_calc_ind as calc_ind
import actin_save_data as save_data
import actin_plot_time as plot
import actin_functions as func


# Configuration file:
config_file = '%s/config_lines.txt' % path


def actin_file(file, calc_index, config_file=config_file, save_output=False, line_plots=False, obj_name=None, targ_list=None, del_out=False, weight=None, norm='npixels'):
    """
    Runs ACTIN for one fits file.
    Accepts files of types: 'e2ds', 's1d', 's1d_*_rv', 'ADP', and 'rdb'.
    Recognizes fits files from HARPS and HARPS-N instruments.

    Parameters:
    -----------
    The same as the actin function below but for one file.
    """
    
    print "\n--------------------"
    print "EXECUTING ACTIN_FILE"
    print "--------------------"

    if type(file) is list: file = file[0]

    # Checking if object in targ_list is the same as the object in fits file
    if targ_list is not None:
        targ_test = func.check_targ(file, targets=targ_list)
        if targ_test is False:
            print "*** ERROR: actin: %s does not belong to any star in the targets list." % file
            return

    data = {}

	# Read config file and retrieve lines information
    if calc_index is not None:
        sel_lines = actin_config.read(config_file, calc_index)
        if sel_lines is None: return

	# Get necessary files
    if "e2ds" in file and ".fits" in file:
        files = read_e2ds.read_file_e2ds(file, obj_name=obj_name)
        file_type = 'e2ds'
    elif "s1d" in file and ".fits" in file:
        files = read_s1d.read_file_s1d(file, obj_name=obj_name)
        file_type = 's1d'
    elif "ADP" in file and ".fits" in file:
        files = read_adp.read_file_adp(file, obj_name=obj_name)
        file_type = 'ADP'
    elif ".rdb" in file:
        file_type = 'rdb'
    else:
        print "*** ERROR: actin: Cannot recognise file type"
        return

    if not files: return

	# Check output file for duplicates
    if save_output is not False and file_type in ('e2ds', 's1d', 'ADP'):
        dupl = save_data.check_duplicate(files['obj'], files['date'], file_type, out_dir=save_output)
        if dupl is True: return

	# Load spectral data
    if file_type == 'e2ds':
        data = read_e2ds.load_data_e2ds(files)

    elif file_type == 's1d':
        data = read_s1d.load_data_s1d(files)

    elif file_type == 'ADP':
        data = read_adp.load_data_adp(files)

    elif file_type == 'rdb':
		# the table must have as headers: 'obj','date','bjd','wave','flux','error_pixel' (optional)
        print "\nLOADING DATA FROM .rdb TABLE"
        print "----------------------------"
        print "READING FILE:\t%s" % file
        data = func.read_rdb(file)[0]
        #else: data['file_orig'] = file.split('/')[-1]
        print "OBJECT:\t\t%s" % data['obj'][0]
        if save_output is True:
            dupl = save_data.check_duplicate(data['obj'][0], data['date'][0], file_type, out_dir=save_output)
            if dupl is True: return

    else: pass

    if not data: return

    data['file_type'] = file_type

    if calc_index is not None:
        # Check selected lines for spectral range and orders
        calc_ind.check_lines(data['wave'], sel_lines)

        # Calculate flux in the required lines
        sel_lines = calc_ind.calc_flux_lines(data, sel_lines, save_plots=line_plots, weight=weight, norm=norm)

        # Calculate chosen indices
        index = calc_ind.calc_ind(sel_lines)

    else:
        print "\nNo indices selected. Insert an index with the '-i' option. Index ids are in the column 'ind_id' of the config_lines.txt file."
        index = None

	# Write output to rdb file in "out_dir"/"obj"
    if save_output is not False:
        save_name = save_data.save_data(data, index, out_dir=save_output)
    elif save_output is False:
        save_name = None

    return data, index, save_name


def actin(files, calc_index, config_file=config_file, save_output=False, line_plots=False, obj_name=None, targ_list=None, del_out=False, weight=None, norm='npixels'):
    """
    Runs ACTIN for one or multiple fits files, for one or multiple stars.
    Accepts files of types: 'e2ds', 's1d', 's1d_*_rv', 'ADP', and 'rdb'.
    Recognizes fits files from HARPS and HARPS-N instruments.

    Parameters:
    -----------
    files : list
        File(s) to be read by ACTIN. Can be of type 'e2ds', 's1d', 's1d_*_rv', 'ADP', or 'rdb' table.
    calc_index : list
        Indices to be calculated. Indices identification must be consistent
        with the ones given in the configuration file.
    config_file : str (optional)
        Configuration file with path. Default is 'config_lines.txt'. If
        this file is moved from the currect directory, the new path should
        be given by this option.
    save_output : {str, False} (optional)
        Directory to save the output. If False, ACTIN is run but data is
        not saved (default).
    line_plots : {str, False} (optional)
        Directory to save the spectral line plots. If 'same', uses the
        directory specified in 'save_output'. If False, plots are not saved
        (default).
    obj_name : {str, None} (optional)
        Object name to overide the one from the fits files. Useful in the
        case where different object names are given in the fits files (ex:
        Gl551, Proxima, ProximaCen). ACTIN will use this name as the save
        directory for the output. If None, uses identification from fits
        (default). NOTE: use this carefully, only if certain that the files
        belong to the same object.
    targ_list : {list, None} (optional)
        List of objects to search the 'files' given as input. ACTIN will run only when the fits files 'obj' id match names in the list. If None,
        ACTIN runs all 'files' given as input. NOTE: if the 'files' list is
        very large, it can take some time to run this option.
    del_out : bool (optional)
        If True, searchs and deletes all output rdb files belonging to the
        objects identified in 'files' and previously saved by ACTIN. Usefull
        when running ACTIN on the same fits files but with different options.
        If False, ACTIN will search for duplicates in the output rdb file and ignore them if found or add data if there are no duplicates (default).
    weight : {str, None} (optional)
        Function to weight the integrated flux. If 'blaze' the flux is multiplied by the blaze function, if None (default) the flux is not weighted (default).
    norm : str (optional)
    	Normalization of the flux: if 'band' the sum is normalized by the bandpass wavelength value in angstroms, if 'npixels' by the number of pixels in the bandpass (default), if 'weight' by the sum of the weight function inside the bandpass, if None the integrated flux is not normalized.
    """

    print "\n#----------------#"
    print "# STARTING ACTIN #"
    print "#----------------#"

    start_time = time.time()

    if files is None:
        print "*** ERROR: actin: No file(s) specified. File name(s) should be the first argument of actin.py."
        return

    if weight == 'None': weight = None # str gives problems in comput_flux
    if norm == 'None': norm = None

    # Remove output file
    if del_out is True:
        print "\nSearching output files to delete..."
        err_msg = func.remove_output(files, save_output=save_output, targ_list=targ_list)

    # Option to make line plots directory the same as the data output dir
    if line_plots == 'same': line_plots = save_output

    # Case of only one file in files list [if called from terminal]
    if type(files) is list and len(files) == 1:
        files = files[0]

    # Case of files is one files (str) [if called from script]
    if type(files) is str:
        # run ACTIN for one file
        actin_file(files, calc_index, config_file=config_file, save_output=save_output, line_plots=line_plots, obj_name=obj_name, targ_list=targ_list, del_out=del_out, weight=weight, norm=norm)

    elif type(files) is list and len(files) > 1:
        # run ACTIN for list of files
        output_rdb = []
        for k in range(len(files)):
            output = actin_file(files[k], calc_index, config_file=config_file, save_output=save_output, line_plots=line_plots, obj_name=obj_name, targ_list=targ_list, del_out=del_out, weight=weight, norm=norm)[2]
            output_rdb.append(output)

        # save plot timeseries
        if output_rdb and calc_index is not None:
            output_rdb = list(set(output_rdb)) # remove duplicates in list
            for k in range(len(output_rdb)):
                plot.plt_time(output_rdb[k], rmv_flgs=False, save_plt=True)
                plot.plt_time_mlty(output_rdb[k], rmv_flgs=False, save_plt=True)

    elapsed_time = time.time() - start_time

    print "\n---------------------------"
    if type(files) is str: nfiles = 1
    else: nfiles = len(files)
    print "Weight used in flux:\t%s" % weight
    print "Normalization of flux:\t%s" % norm
    print "Files analysed:\t\t%i" % nfiles
    print "Save output:\t\t%s" % save_output
    print "Elapsed time:\t\t%.4f min" % (float(elapsed_time)/60.)

    return


if __name__ == "__main__":

    if len(sys.argv) < 2:
        print "You haven't specified any arguments. Use -h to get more details on how to use this command."
        sys.exit(1)

    import argparse

	# initiate the parser
    parser = argparse.ArgumentParser()

	# add short and long argument
    parser.add_argument('files', help='Read file(s)', nargs='+')

    parser.add_argument('--index', '-i', help="Index id to calculate as designated by 'ind_id' in config_index.txt.", nargs='+', default=None)

    parser.add_argument('--save_data', '-s', help='Path to output directory of data table, or False (default).', default=False)#, type=lambda x: (str(x).lower() == 'true'))

    parser.add_argument('--save_plots', '-p', help='Path to output directory of line plots, or False (default)', default=False)#, type=lambda x: (str(x).lower() == 'true'))

    parser.add_argument('--obj_name', '-obj', help='Give target a name that overides the one from the fits files.', default=None)

    parser.add_argument('--targ_list', '-tl', help='Give a list of stars to select from fits files.', nargs='+', default=None)

    parser.add_argument('--del_out', '-del', help='Delete output data file if True.', default=False, type=lambda x: (str(x).lower() == 'true'))

    parser.add_argument('--weight', '-w', help='Weight for index calculation. Allowed values: None (default), blaze.', default=None)

    parser.add_argument('--norm', '-n', help='Normalization of flux. Allowed values: band, npixels (default), weight, None.', default='npixels')

	# read arguments from the command lines
    args = parser.parse_args()

    actin(files=args.files, calc_index=args.index, save_output=args.save_data, line_plots=args.save_plots, obj_name=args.obj_name, targ_list=args.targ_list, del_out=args.del_out, weight=args.weight, norm=args.norm)
