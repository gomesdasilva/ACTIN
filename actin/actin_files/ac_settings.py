#!/usr/bin/env python

# compatibility with python 2/3:
from __future__ import print_function
from __future__ import division

import numpy as np


def init():

    # Accepted file types
    global ftypes
    ftypes = {}
    ftypes['1d'] = ['S1D',
                    's1d',
                    'ADP',
                    'rdb']

    ftypes['2d'] = ['S2D',
                    'e2ds']

    ftypes['all'] = ftypes['1d'] + ftypes['2d']

    # Accepted instruments
    global instr
    instr =    ['HARPS',
                'HARPN',
                'ESPRESSO']

    # File names
    global fnames
    fnames = {}
    fnames['data'] = "data.rdb"
    fnames['log_data'] = "log.txt"
    fnames['lines_data'] = "lines.txt"
    fnames['ln_plt'] = ".pdf"
    fnames['time_plt'] = "time.pdf"
    fnames['time_mlty_plt'] = "time_mlty.pdf"

    # output keys (used in ac_seve.py)
    global outkeys
    outkeys = {}
    outkeys = ['obj', 'instr', 'obs_date', 'bjd', 'rv', 'rv_err', 'fwhm', 'fwhm_err', 'cont', 'cont_err', 'bis', 'bis_err', 'noise', 'median_snr', 'data_flg', 'bv', 'airmass_end']

    # err messages
    global err
    err = {}
    err['ERROR'] = "*** ERROR:"


def preamble(version_file):

    __author__ = "Joao Gomes da Silva"

    try:
        with open(version_file, 'r') as file:
            version = file.read().splitlines()[0]
        print("\nACTIN {}".format(version))
    except:
        print("*** WARNING | Unable to read 'VERSION' file.")
        version = "Unknown"

    print("Instituto de Astrofisica e Ciencias do Espaco")
    print("Centro de Astrofisica da Universidade do Porto")
    print("Author:", __author__+ ",", "Joao.Silva@astro.up.pt")

    return version
