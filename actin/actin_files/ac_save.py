#!/usr/bin/env python

# compatibility with python 2/3:
from __future__ import print_function
from __future__ import division

import sys, os
import pylab as plt
import numpy as np

import datetime

# ACTIN FILES:
import ac_tools
import ac_settings as ac_set


def check_duplicate(obj, date, instr, file_type, out_dir):
    """
    Check if measurement is a duplicate in the output file.

    Parameters:
    -----------
    obj : str
        Object identification.
    date : str
        Date of observation in the fits file format.
    file_type : str
        Type of file used: 'S2D', 'S1D', 'e2ds', 's1d', 'ADP', or 'rdb'.
    output_dir : str
        Directory of output file.

    Returns:
    --------
    bool
        True if duplicate, False otherwise.
    """

    print()
    print("Executing check_duplicate:")

    if obj is None or date is None or out_dir is None: return

    file_name = "{}_{}_{}_{}".format(obj, instr, file_type, ac_set.fnames['data'])
    pfile_name = os.path.join(out_dir, obj, file_name)

    if os.path.isfile(pfile_name):
        try:
            rdb_data = ac_tools.read_rdb(pfile_name)[0]
        except:
            print("*** ERROR: Cannot read rdb file:")
            print(pfile_name)
            sys.exit()
        if date in rdb_data['obs_date']:
            print("Date {} already saved in:".format(date))
            print(pfile_name)
            print("*** ACTION: Ignoring measurement")
            return True
        else:
            print("Date {} not present in {}".format(date, pfile_name))
    else: print("No data saved for {}, {}.".format(obj, date))

    return False


def save_data(data, index, out_dir):
    """
    Saves output data to rdb file.
    For a description of the returned headers see README.md file.
    """

    print()
    print("SAVING DATA")
    print("-----------")

    print("Executing save_data")

    if data is None or out_dir is None:
        print("*** ERROR: 'data' is None or 'out_dir' is None")
        return

    # reading from rdb file
    if data['file_type'] == 'rdb':
        data_keys = ['obj', 'obs_date', 'bjd']
    else:
        data_keys = ac_set.outkeys

    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    if not os.path.isdir(os.path.join(out_dir, data['obj'])):
        os.mkdir(os.path.join(out_dir, data['obj']))

    # Include index information
    if index is not None:
        # convert index data from table list to dictionary
        indices = [index['index'][k] for k in range(len(index['index']))]
        ind = {}
        for k in range(len(indices)):
            ind[indices[k]] = index['value'][k]
            ind['{}_err'.format(indices[k])] = index['error'][k]
            ind['{}_snr'.format(indices[k])] = index['snr'][k]
            ind['{}_flg'.format(indices[k])] = index['flg'][k]
            ind['{}_mfrac_neg'.format(indices[k])] = index['mfrac_neg'][k]

        # For npixels in each line
        index_keys = list(index)
        for k in range(len(index_keys)):
            if index_keys[k].split('_')[-1] == 'npixels':
                ind['{}_npixels'.format(index_keys[k].split('_')[0])] = index['{}_npixels'.format(index_keys[k].split('_')[0])]

        #index_keys = list(ind.keys())
        index_keys = list(ind)
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

    file_data = "{}_{}_{}_{}".format(all_data['obj'], data['instr'], data['file_type'], ac_set.fnames['data'])
    pfile_data = os.path.join(out_dir, all_data['obj'], file_data)

    if not os.path.isfile(pfile_data):
        ac_tools.save_rdb(data_save, keys, pfile_data)
        print("Data saved to:\n{}".format(pfile_data))
    else:
        ac_tools.add_rdb(data_save, keys, pfile_data)
        print("Data added to:\n{}".format(pfile_data))

    return pfile_data


def save_log(info, options, n_files, out_dir):

    if out_dir is False: return

    print()
    print("Executing save_log:")

    file_type = info['file_type']
    instr = info['instr']
    files = info['source_path']
    obj = info['obj']

    file_log = "{}_{}_{}_{}".format(obj, instr, file_type, ac_set.fnames['log_data'])
    pfile_log = os.path.join(out_dir, obj, file_log)

    date_now = str(datetime.datetime.now()).split(".")[0]

    if type(files) is list and len(files) > 1:
        source_path = os.path.commonprefix(files)
    else: source_path = os.path.split(files)[0]

    log = {}
    log['version'] = info['version']
    log['run_date'] = date_now
    log['source_path'] = info['source_path']
    log['tel'] = info['tel']
    log['instr'] = info['instr']
    log['file_type'] = info['file_type']
    log['n_files'] = n_files
    log['frac'] = options['frac']

    with open(pfile_log,'w') as log_file:
        for k in range(len(list(log))):
            print("{}={}".format(list(log)[k],log[list(log)[k]]), file=log_file)

    print("Log saved to:\n{}".format(pfile_log))

    return pfile_log


def save_line_info(info, sel_lines, out_dir):

    if out_dir is False: return
    if not sel_lines: return

    print()
    print("Executing save_line_info:")

    obj = info['obj']
    file_type = info['file_type']
    instr = info['instr']

    keys = ["ind_id", "ind_var", "ln_id", "ln_c", "ln_ctr", "ln_win", "bandtype"]

    file_lines = "{}_{}_{}_{}".format(obj, instr, file_type, ac_set.fnames['lines_data'])
    pfile_lines = os.path.join(out_dir, obj, file_lines)

    ac_tools.save_rdb(sel_lines, keys, pfile_lines)
    print("Lines information saved to:\n{}".format(pfile_lines))

    return pfile_lines
