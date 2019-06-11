#!/usr/bin/env python

# compatibility with python 2/3:
from __future__ import print_function
from __future__ import division

import sys, os
import numpy as np
import pylab as plt

import ac_settings as ac_set
import ac_tools


# IMPORTANT: all indices names must start with "I_" in config file
def plt_time(info, out_dir, save_plt=False, rmv_flgs=False):
    """
    Saves timeseries plots of the indices identified in the rdb file by
    starting with 'I_'.
    """

    if out_dir is False: return

    print()
    print("Executing plt_time:")

    obj = info['obj']
    file_type = info['file_type']
    instr = info['instr']

    rdb_name = "{}_{}_{}_{}".format(obj, instr, file_type, ac_set.fnames['data'])
    rdb_file = os.path.join(out_dir, obj, rdb_name)

    width, height = ac_tools.plot_params()

    data = ac_tools.read_rdb(rdb_file)[0]

    bjd_raw = np.asarray(data['bjd'])
    bjd = bjd_raw - int(min(bjd_raw))

    if len(bjd) == 1:
        print("Only one data point, no need to plot.")
        return

    data_keys = list(data)

    ind_ids = []
    for k in range(len(data_keys)):
        if any(map(data_keys[k].startswith, 'I_')):
            ind_ids.append(('_').join(data_keys[k].split('_')[:2]))

    if not ind_ids:
        print("*** ERROR: No indices detected in file:")
        print("***", rdb_file)
        return

    # remove duplicates of ind_id and gives a list of selected indices
    ind_ids = list(set(ind_ids))

    # make new dictionary just for plots
    ind = {}
    for k in range(len(ind_ids)):
        ind[ind_ids[k]] = data[ind_ids[k]]
        ind_err_key = '{}_err'.format(ind_ids[k])
        ind_flg_key = '{}_flg'.format(ind_ids[k])
        ind[ind_err_key] = data[ind_err_key]
        ind[ind_flg_key] = data[ind_flg_key]

    for k in range(len(ind_ids)):
        plt.figure(figsize=(width,height))

        N = 0
        for i in range(len(ind[ind_ids[k]])):
            ind_flg = ind['{}_flg'.format(ind_ids[k])][i]
            ind_err = ind['{}_err'.format(ind_ids[k])][i]
            ind_val = ind[ind_ids[k]][i]
            if ind_flg == 'None':
                N += 1
                plt.errorbar(bjd[i], ind_val, ind_err, c='k', marker='.', ls='')
            if ind_flg == 'negFlux':
                if rmv_flgs is False:
                    N += 1
                    plt.errorbar(bjd[i], ind_val, ind_err, c='r', marker='.', ls='')
                else: pass

        plt.annotate(obj, xy=(0.05,0.9), xycoords='axes fraction', textcoords='axes fraction')

        plt.annotate("N = {}".format(N), xy=(0.05,0.85), xycoords='axes fraction', textcoords='axes fraction')

        plt.xlabel('BJD - {} [days]'.format(int(min(bjd_raw))))
        plt.ylabel(ind_ids[k])

        if save_plt is True:
            save_name = '{}_{}_{}_{}_{}'.format(obj, instr, file_type, ind_ids[k], ac_set.fnames['time_plt'])
            save_file = os.path.join(out_dir, obj, save_name)
            plt.savefig(save_file)
            print("{} timeseries saved to:\n{}".format(ind_ids[k], save_file))

        plt.close()

    return


def plt_time_mlty(info, out_dir, save_plt=False, rmv_flgs=False, hdrs=['I_CaII', 'I_Ha', 'I_NaI', 'I_HeI']):
    """
    Saves timeseries plots of the indices identified in the rdb file by
    starting with 'I_' in a 'multi-plot' format.
    """

    if not out_dir: return

    print()
    print("Executing plt_time_mlty:")

    if not info:
        print("*** ERROR: There is no 'info' dictionary.")
        return
    if not hdrs:
        print("*** ERROR: No input indices. Use '-i' and insert index names that are listed in the config file. To know location of config file call ACTIN without any arguments. To read a config file from a different location use '-cf' and add the path/file.")
        return
    if len(hdrs) == 1:
        print("Only one index selected, no need for multiplot")
        return

    obj = info['obj']
    file_type = info['file_type']
    instr = info['instr']

    rdb_name = "{}_{}_{}_{}".format(obj, instr, file_type, ac_set.fnames['data'])
    rdb_file = os.path.join(out_dir, obj, rdb_name)

    data = ac_tools.read_rdb(rdb_file)[0]

    bjd_raw = np.asarray(data['bjd'])
    bjd = bjd_raw - int(min(bjd_raw))

    if len(bjd) == 1:
        print("Only one data point, no need to plot.")
        return

    ind = {} # make new dictionary just for plots
    ind_ids = [] # only the indices in data
    for k in range(len(hdrs)):
        if hdrs[k] in list(data):
            ind[hdrs[k]] = data[hdrs[k]]
            ind_err_key = '{}_err'.format(hdrs[k])
            ind_flg_key = '{}_flg'.format(hdrs[k])
            ind[ind_err_key] = data[ind_err_key]
            ind[ind_flg_key] = data[ind_flg_key]
            ind_ids.append(hdrs[k])

    if len(ind_ids) == 1:
        print("Only one index in rdb file, no need for multiplot.")
        return

    if len(ind_ids) > 5:
        print("Too many indices for multiplot.")
        return

    width, height = ac_tools.plot_params(7, 4.5)
    height = 2*len(ind_ids)

    plot_n = 0
    n_subplots = len(ind_ids)
    plt.figure(figsize=(width, height))
    for k in range(len(ind_ids)):
        plot_n += 1
        ax1 = plt.subplot(n_subplots, 1, plot_n)

        N = 0
        for i in range(len(ind[ind_ids[k]])):
            ind_flg = ind['{}_flg'.format(ind_ids[k])][i]
            ind_val = ind[ind_ids[k]][i]
            ind_err = ind['{}_err'.format(ind_ids[k])][i]
            if ind_flg == 'None':
                N += 1
                plt.errorbar(bjd[i], ind_val, ind_err, c='k', marker='.',ls='')
            if ind_flg == 'negFlux':
                if rmv_flgs == False:
                    N += 1
                    plt.errorbar(bjd[i], ind_val, ind_err, c='r', marker='.',ls='')
                else: pass

        plt.xlabel('BJD - {} [days]'.format(int(min(bjd_raw))))
        plt.ylabel(ind_ids[k])

        from matplotlib.ticker import MaxNLocator
        ax1.yaxis.set_major_locator(MaxNLocator(prune='both'))

        if plot_n != len(ind_ids):
            frame1 = plt.gca()
            frame1.axes.get_xaxis().set_ticklabels([])

    n = len(ind_ids)
    plt.annotate(obj, xy=(0.05,1.*n+0.05), xycoords='axes fraction', textcoords='axes fraction')

    plt.annotate("N = {}".format(N), xy=(0.85,1.*n+0.05), xycoords='axes fraction', textcoords='axes fraction')

    plt.subplots_adjust(hspace=0.000)

    if save_plt == True:
        save_name =  '{}_{}_{}_{}_{}'.format(obj, instr, file_type, ('_').join(ind_ids), ac_set.fnames['time_mlty_plt'])
        save_file = os.path.join(out_dir, obj, save_name)
        plt.savefig(save_file)
        print("{} multiplot saved to:\n{}".format((', ').join(ind_ids), save_file))

    plt.close()

    return
