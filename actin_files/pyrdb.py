#(@) $Id: pyrdb.py,v 3.0 2007/01/17 13:22:12 vltsccm Exp $
#
# who       when      what
# --------  --------  ----------------------------------------------
# dqueloz	27/10/05	


import os,sys,string,pickle
from numpy import *


def read_rdb(filename):
# use: table = pyrdb.read_rdb(file)[0] for data
# use: table = pyrdb.read_rdb(file)[1] to get the keys by order


	f = open(filename, 'r')
	data = f.readlines()
	f.close()
	
	key = string.split(data[0][:-1],'\t')
	output = {}
	for i in range(len(key)): output[key[i]] = []
	
	for line in data[2:]:
		qq = string.split(line[:-1],'\t')
		for i in range(len(key)):
			try: value = float(qq[i])
			except ValueError: value = qq[i]
			output[key[i]].append(value)
	
	return output,key

# edit by J. Gomes da Silva
def save_rdb(dic,keys,file):
	
	out = open(file,'w')
	n_keys = len(keys)
		
	for k in range(n_keys):
		if k != n_keys-1: out.write('%s\t' % (keys[k]))
		else: out.write('%s' % (keys[k]))
	out.write('\n')
	for k in range(n_keys):
		if k != n_keys-1: out.write('%s\t' % ('-'*len(keys[k])))
		else: out.write('%s' % ('-'*len(keys[k])))
	out.write('\n')
		
	for i in range(len(dic[keys[0]])):
		for k in range(n_keys):
			if k != n_keys-1: out.write('%s\t' % (dic[keys[k]][i]))
			else: out.write('%s' % (dic[keys[k]][i]))
		out.write('\n')
	out.close()

# edit by J. Gomes da Silva
def add_rdb(dic,keys,file_name):
    
    out = open(file_name,'a')
    n_keys = len(keys)
    for i in range(len(dic[keys[0]])):
        for k in range(n_keys):
            if k != n_keys-1: out.write('%s\t' % (dic[keys[k]][i]))
            else: out.write('%s\t' % (dic[keys[k]][i]))	
        out.write('\n')
    out.close()


def write_rdb(filename,data,keys,format):

	f = open(filename, 'w')
	
	head1 = string.join(keys,'\t')
	head2 = ''
	for i in head1:
		if i=='\t': head2 = head2+'\t'
		else: head2 = head2+'-'
	
	f.write(head1+'\n')
	f.write(head2+'\n')
	
	if len(data.values()) > 0:
		for i in range(len(data.values()[0])):
			line = []
			for j in keys: line.append(data[j][i])
			f.write(format % tuple(line))
	
	f.close()


def read_rdb_rows(filename,refcol):
	
	f = open(filename, 'r')
	data = f.readlines()
	f.close()
	
	key = string.split(data[0][:-1],'\t')
	iref = key.index(refcol)
	output = {}
	
	for line in data[2:]:
		qq1 = string.split(line[:-1],'\t')
		qq2 = {}
		for i in range(len(key)): qq2[key[i]] = qq1[i]
		output[qq1[iref]] = qq2
	
	return output
    
def ajustement_lenght_for_write(vecteur,max_len):
    vecteur_write=array([1e30]*max_len)
    for i in arange(0,len(vecteur),1):
        vecteur_write[i]=vecteur[i]
    return vecteur_write
    
    
def read_to_long_vecteur(vecteur):
    for i in arange(0,len(vecteur),1):
        if vecteur[i]==1e30:
            i=i-1
            break
    vecteur=vecteur[:i+1]
    return vecteur

