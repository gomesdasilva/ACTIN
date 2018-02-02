#!/usr/bin/env python

import sys, os
import numpy as np
import pylab as plt

# Location of SPECHA files:
sys.path.append("specha_files/")

# SPECHA files
import pyrdb
import functions as func



# IMPORTANT: all indices names must start with "I_" in config file
def plt_time(rdb_file, save_plt=False, rmv_flgs=False):
    print "\nPLOTTING TIMESERIES"
    print "-------------------"

    if rdb_file == None: return

    if type(rdb_file) is list: rdb_file = rdb_file[0]

    width, height = func.plot_params()

    star = rdb_file.split('/')[-1].split('_')[0]
    dir = ('/').join(rdb_file.split('/')[:-1])

    data = pyrdb.read_rdb(rdb_file)[0]

    bjd = np.asarray(data['bjd']) - 2450000

    ind_ids = []
    for k in range(len(data.keys())):
        if 'I_' in data.keys()[k]:
            ind_ids.append(('_').join(data.keys()[k].split('_')[:2]))

    if len(ind_ids) == 0:
        print "No indices detected in %s" % rdb_file
        return

    # remove duplicates of ind_id and gives a list of selected indices
    ind_ids = list(set(ind_ids))

    # make new dictionary just for plots
    ind = {}
    for k in range(len(ind_ids)):
        ind[ind_ids[k]] = data[ind_ids[k]]
        ind['%s_err' % ind_ids[k]] = data['%s_err' % ind_ids[k]]
        ind['%s_flg' % ind_ids[k]] = data['%s_flg' % ind_ids[k]]

    # Check s_mw
    if 's_mw' in data.keys():
        ind_ids.append('s_mw')
        ind['s_mw'] = data['s_mw']
        ind['s_mw_err'] = data['s_mw_err']
        ind['s_mw_flg'] = data['s_mw_flg']
    else: print "No s_mw detected in %s" % rdb_file

    for k in range(len(ind_ids)):
        plt.figure(figsize=(width,height))

        N = 0
        for i in range(len(ind[ind_ids[k]])):
            if '%s_flg' % ind_ids[k] in ind.keys() and ind['%s_flg' % ind_ids[k]][i] == 'None':
                N += 1
                plt.errorbar(bjd[i], ind[ind_ids[k]][i], ind['%s_err' % ind_ids[k]][i], c='k', marker='.',ls='')
            elif ind['%s_flg' % ind_ids[k]][i] != 'None':
                if rmv_flgs == False:
                    N += 1
                    plt.errorbar(bjd[i], ind[ind_ids[k]][i], ind['%s_err' % ind_ids[k]][i], c='r', marker='.',ls='')
                else: pass

        plt.annotate(star, xy=(0.05,0.9), xycoords='axes fraction', textcoords='axes fraction')

        plt.annotate("N = %i" % N, xy=(0.05,0.85), xycoords='axes fraction', textcoords='axes fraction')

        plt.xlabel('bjd - 2450000')
        plt.ylabel(ind_ids[k])

        if save_plt == True:
            save_name = '%s_%s_time.pdf' % (star,ind_ids[k])
            plt.savefig(('%s/%s' % (dir,save_name)))
            print "%s timeseries saved to %s/%s" % (ind_ids[k],dir,save_name)

        plt.close()

    return


def plt_time_mlty(rdb_file, save_plt=False, rmv_flgs=False, hdrs=['s_mw', 'I_Ha', 'I_NaI', 'I_HeI']):

    if rdb_file == None: return

    if type(rdb_file) is list: rdb_file = rdb_file[0]

    width, height = func.plot_params(7, 6.5)

    star = rdb_file.split('/')[-1].split('_')[0]
    dir = ('/').join(rdb_file.split('/')[:-1])

    data = pyrdb.read_rdb(rdb_file)[0]

    bjd = np.asarray(data['bjd']) - 2450000

    if len(hdrs) == 0:
        print "No indices selected in plt_time_mlty"
        return

    ind = {} # make new dictionary just for plots
    ind_ids = [] # only the indices in data
    for k in range(len(hdrs)):
        if hdrs[k] in data.keys():
            ind[hdrs[k]] = data[hdrs[k]]
            ind['%s_err' % hdrs[k]] = data['%s_err' % hdrs[k]]
            ind['%s_flg' % hdrs[k]] = data['%s_flg' % hdrs[k]]
            ind_ids.append(hdrs[k])

    if len(ind_ids) == 0: return

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

    plt.annotate("N = %i" % N, xy=(0.2,1.*n+0.05), xycoords='axes fraction', textcoords='axes fraction')

    plt.subplots_adjust(hspace=0.000)
    #plt.subplots_adjust(bottom=0.075, left=0.11, top = 0.95, right=0.975)

    save_name =  '%s_%s_time_mlty.pdf' % (star,('_').join(ind_ids))
    plt.savefig(('%s/%s' % (dir,save_name)))
    print "%s timeseries saved to %s/%s" % ((', ').join(ind_ids),dir,save_name)

    plt.close()

    return


if __name__ == '__main__':

    if len(sys.argv) < 2:
        print "You haven't specified any arguments. Use -h to get more details on how to use this command."
        sys.exit(1)

    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('rdb_file', help='read .rdb file(s)', nargs='+')
    parser.add_argument('--rmv_flgs', '-flg', help='remove flagged data if true', default=False, type=lambda x: (str(x).lower() == 'true'))
    parser.add_argument('--save_plt', '-s', help='save plot if True', default=False, type=lambda x: (str(x).lower() == 'true'))

    args = parser.parse_args()

    plt_time(rdb_file=args.rdb_file, rmv_flgs=args.rmv_flgs, save_plt=args.save_plt)
    plt_time_mlty(rdb_file=args.rdb_file, rmv_flgs=args.rmv_flgs, save_plt=args.save_plt)
