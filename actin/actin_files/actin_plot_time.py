#!/usr/bin/env python

# compatibility with python 2/3:
from __future__ import print_function
from __future__ import division

import sys, os
import numpy as np
import pylab as plt

import actin_functions as func



# IMPORTANT: all indices names must start with "I_" in config file
def plt_time(rdb_file, save_plt=False, rmv_flgs=False):
    """
    Saves timeseries plots of the indices identified in the rdb file by
    starting with 'I_'.

    Parameters:
    -----------
    rdb_file : str
        Output rdb file from ACTIN with path.
    save_plt : bool (optional)
        If True the plots are saved in the same directory of 'rdb_file',
        if False the plots are not saved (default).
    rmv_flgs : bool (optional)
        If True, removes values with negative flux flags, if False the
        these values are included with red colour (default).

    Returns:
    --------
    Timeseries plots of the indices selected to be calculated by ACTIN.
    Output directory is the same as the one the rdb file is located.
    """

    print("\nPLOTTING TIMESERIES")
    print("-------------------")

    if rdb_file is None:
        print("*** ERROR: No rdb file given to plot timeseries.")
        return

    if type(rdb_file) is list: rdb_file = rdb_file[0]

    folder, file = os.path.split(rdb_file)
    star = file.split('_')[0]
    file_type = file.split('_')[1]

    width, height = func.plot_params()

    data = func.read_rdb(rdb_file)[0]

    bjd = np.asarray(data['bjd']) - 2450000

    ind_ids = []
    for k in range(len(data.keys())):
        if 'I_' in data.keys()[k]:
            ind_ids.append(('_').join(data.keys()[k].split('_')[:2]))

    if not ind_ids:
        print("No indices detected in %s" % rdb_file)
        return

    # remove duplicates of ind_id and gives a list of selected indices
    ind_ids = list(set(ind_ids))

    # make new dictionary just for plots
    ind = {}
    for k in range(len(ind_ids)):
        ind[ind_ids[k]] = data[ind_ids[k]]
        ind['%s_err' % ind_ids[k]] = data['%s_err' % ind_ids[k]]
        ind['%s_flg' % ind_ids[k]] = data['%s_flg' % ind_ids[k]]

    for k in range(len(ind_ids)):
        plt.figure(figsize=(width,height))

        N = 0
        for i in range(len(ind[ind_ids[k]])):
            if '%s_flg' % ind_ids[k] in ind.keys() and ind['%s_flg' % ind_ids[k]][i] == 'None':
                N += 1
                plt.errorbar(bjd[i], ind[ind_ids[k]][i], ind['%s_err' % ind_ids[k]][i], c='k', marker='.',ls='')
            elif ind['%s_flg' % ind_ids[k]][i] != 'None':
                if rmv_flgs is False:
                    N += 1
                    plt.errorbar(bjd[i], ind[ind_ids[k]][i], ind['%s_err' % ind_ids[k]][i], c='r', marker='.',ls='')
                else: pass

        plt.annotate(star, xy=(0.05,0.9), xycoords='axes fraction', textcoords='axes fraction')

        plt.annotate("N = %i" % N, xy=(0.05,0.85), xycoords='axes fraction', textcoords='axes fraction')

        plt.xlabel('bjd - 2450000')
        plt.ylabel(ind_ids[k])

        if save_plt is True:
            save_name = '%s_%s_%s_time.pdf' % (star, file_type, ind_ids[k])
            plt.savefig(os.path.join(folder, save_name))
            print("%s timeseries saved to %s" % (ind_ids[k], os.path.join(folder, save_name)))

        plt.close()

    return


def plt_time_mlty(rdb_file, save_plt=False, rmv_flgs=False, hdrs=['I_CaII', 'I_Ha', 'I_NaI', 'I_HeI']):
    """
    Saves timeseries plots of the indices identified in the rdb file by
    starting with 'I_' in a 'multi-plot' format.

    Parameters:
    -----------
    rdb_file : str
        Output rdb file from ACTIN with path.
    save_plt : bool (optional)
        If True the plots are saved in the same directory of 'rdb_file',
        if False the plots are not saved (default).
    rmv_flgs : bool (optional)
        If True, removes values with negative flux flags, if False the
        these values are included with red colour (default).
    hdrs : list
        List of indices to be plotted. If the rdb does not include any
        index in this list, plot is not saved. If only one index is found
        in the rdb file, plot is not saved.

    Returns:
    --------
    Timeseries multi-plot of the indices selected to be calculated by ACTIN.
    Output directory is the same as the one the rdb file is located.
    """

    if rdb_file is None: return

    if type(rdb_file) is list: rdb_file = rdb_file[0]

    folder, file = os.path.split(rdb_file)
    star = file.split('_')[0]
    file_type = file.split('_')[1]

    width, height = func.plot_params(7, 6.5)

    data = func.read_rdb(rdb_file)[0]

    bjd = np.asarray(data['bjd']) - 2450000

    if len(hdrs) == 0:
        print("No indices selected in plt_time_mlty")
        return

    if len(hdrs) == 1:
        print("Only one index selected, no need for multiplot")
        return

    ind = {} # make new dictionary just for plots
    ind_ids = [] # only the indices in data
    for k in range(len(hdrs)):
        if hdrs[k] in data.keys():
            ind[hdrs[k]] = data[hdrs[k]]
            ind['%s_err' % hdrs[k]] = data['%s_err' % hdrs[k]]
            ind['%s_flg' % hdrs[k]] = data['%s_flg' % hdrs[k]]
            ind_ids.append(hdrs[k])

    if len(ind_ids) == 0:
        print("No index matches the rdb file and hdr option of plt_time_mlty, must have data for one of the following indices: I_CaII, I_Ha, I_NaI or I_HeI.")
        return
    if len(ind_ids) == 1:
        print("Only one index in rdb file, no need for multiplot.")
        return

    plot_n = 0
    n_subplots = len(ind_ids)
    plt.figure(figsize=(width, height))
    for k in range(len(ind_ids)):
        plot_n += 1
        ax1 = plt.subplot(n_subplots, 1, plot_n)

        N = 0
        for i in range(len(ind[ind_ids[k]])):
            if ind['%s_flg' % ind_ids[k]][i] == 'None':
                N += 1

                plt.errorbar(bjd[i], ind[ind_ids[k]][i], ind['%s_err' % ind_ids[k]][i], c='k', marker='.',ls='')
            elif ind['%s_flg' % ind_ids[k]][i] != 'None':
                if rmv_flgs == False:
                    N += 1

                    plt.errorbar(bjd[i], ind[ind_ids[k]][i], ind['%s_err' % ind_ids[k]][i], c='r', marker='.',ls='')
                else: pass

        plt.xlabel('bjd - 2450000')
        plt.ylabel(ind_ids[k])

        from matplotlib.ticker import MaxNLocator
        ax1.yaxis.set_major_locator(MaxNLocator(prune='both'))

        if plot_n != len(ind_ids):
            frame1 = plt.gca()
            frame1.axes.get_xaxis().set_ticklabels([])

    n = len(ind_ids)
    plt.annotate(star, xy=(0.05,1.*n+0.05), xycoords='axes fraction', textcoords='axes fraction')

    plt.annotate("N = %i" % N, xy=(0.4,1.*n+0.05), xycoords='axes fraction', textcoords='axes fraction')

    plt.subplots_adjust(hspace=0.000)

    save_name =  '%s_%s_%s_time_mlty.pdf' % (star, file_type, ('_').join(ind_ids))
    plt.savefig(os.path.join(folder, save_name))
    print("%s multiplot saved to %s" % ((', ').join(ind_ids), os.path.join(folder, save_name)))

    plt.close()

    return


if __name__ == '__main__':

    if len(sys.argv) < 2:
        print("You haven't specified any arguments. Use -h to get more details on how to use this command.")
        sys.exit(1)

    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('rdb_file', help='read .rdb file(s)', nargs='+')
    parser.add_argument('--rmv_flgs', '-flg', help='remove flagged data if true', default=False, type=lambda x: (str(x).lower() == 'true'))
    parser.add_argument('--save_plt', '-s', help='save plot if True', default=False, type=lambda x: (str(x).lower() == 'true'))

    args = parser.parse_args()

    plt_time(rdb_file=args.rdb_file, rmv_flgs=args.rmv_flgs, save_plt=args.save_plt)
    plt_time_mlty(rdb_file=args.rdb_file, rmv_flgs=args.rmv_flgs, save_plt=args.save_plt)
