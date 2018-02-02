    #!/usr/bin/env python

import sys, os
import numpy as np

def read(config_file, calc_index):
    """
    Reads data from config file and selects the lines needed for the
    selected indices.

    Parameters:
    -----------
    config_file : string
    	Name of the config file with path.

    calc_index : list of strings
        Index id as provided in the configuration file.

    Returns:
    --------
    sel_lines : dict
    	Dictionary containing the identification of the indices selected and
    	the parameters of the spectral lines required for each index.

    	Each key entry is a list of parameters where the list indices form the
    	rows related to the same spectral line identified with the key 'ln_id'
    	which is related to the spectral index identified by 'ind_id'.

    	The returned keys are:

    	==========  ========================================================
    	keys		Description
    	----------  --------------------------------------------------------
    	ind_id		str : Index identification.
    	ind_c		float : Constant to be multilied in the index equation.
    	ind_var		str : Variable assigned to a given line to be used in
    				the index equation. Ex: 'L1', 'L2', etc, for spectral
    				lines; 'R1', 'R2', etc, for reference lines.
    	ln_id		str : Spectral line identification.
    	ln_ctr 		float : Wavelength of the line centre [angstrom].
    	ln_win 		float : Bandwidth around the line centre to be used in
    				the flux integration [angstrom].
    	==========  ========================================================

    """

    print "\nLOADING DATA FROM CONFIG FILE"
    print "-----------------------------"

    if config_file == None:
    	print "*** ERROR: No config file provided."
    	return

    try: f = open(config_file, 'r')
    except:
    	print "*** ERROR: Cannot read config file."
    	return

    header1 = f.readline()
    header2 = f.readline()
    header3 = f.readline()
    header4 = f.readline()
    header5 = f.readline()
    header6 = f.readline()
    header7 = f.readline()
    header8 = f.readline()
    header9 = f.readline()
    header10 = f.readline() # table headers
    header11 = f.readline() # last "-------" before data

    # convert the lines in the config table into list of strings
    columns = []
    for line in f:
    	if not line.strip(): pass # ignore empty new lines
    	else:
    		line = line.strip() # removes whites spaces at the start and end of line
    		column = line.split() # splits the line into a list of strings
    		columns.append(column)

    f.close()

    # associate each key (header) in table to a column
    lines = {}
    keys = header10.split() # keys as provided in the table in config file
    for k in range(len(keys)):
    	lines[keys[k]] = []
    	for i in range(len(columns)):
    		lines[keys[k]].append(columns[i][k])

    # converts the values in the list to floats and leave the strings
    for k in range(len(lines)):
    	for i, x in enumerate(lines[keys[k]]):
    		try: lines[keys[k]][i] = float(x)
    		except ValueError: pass

    # get indices id from the arguments
    sel_ind = calc_index
    print sel_ind

    # check weather the selected index has lines in config table
    for k in range(len(sel_ind)):
    	if sel_ind[k] not in lines['ind_id']:
    		print "*** ERROR: Index %s is not in the config table." % sel_ind[k]
    	else: pass

    # initiate a dictionary with keys as the headers in table and make empty lists
    sel_lines = {}
    for i in range(len(keys)): sel_lines[keys[i]] = []

    # select the rows that belong to the selected index
    rows = len(lines['ind_id'])
    for k in range(rows):
    	if lines['ind_id'][k] in sel_ind:
    		print lines['ln_id'][k]
    		for i in range(len(keys)):
    			sel_lines[keys[i]].append(lines[keys[i]][k])
    	else: pass

    return sel_lines
