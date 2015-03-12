# -*- coding: utf-8 -*-
"""
Created on Tue Dec 16 15:51:09 2014

@author: dima
"""

import subprocess
import sqlite3
from ReadCoverage import *
import sys
# import pdb
###############################################################################
import argparse
parser = argparse.ArgumentParser('''
''')
parser.add_argument("inFile",
help="input file (.db)")

parser.add_argument("outFile", nargs='?', type=str, default='',
help="output file (.csv); for stdout type '-' ")

parser.add_argument("-l", "--rangel", type=int, default="75",
help="")

parser.add_argument("-s", "--csvseparator", type=str, default= '\t',
help="separator for the input file")

parser.add_argument("-t", "--tag", type=str, default="HL10_ngm_nm4",
help="") 
"ABD192_bwa_nm3"

parser.add_argument("-o", "--outFile", type=str, default= 'out.tab',
help="file for output")

parser.add_argument("-q", "--type", type=str, default= 'H',
help="blob data type (H-- unsigned short int, i -- int, B -- char)")

args = parser.parse_args()
###############################################################################
global rl
rl = args.rangel
subrange = 50;
CHROMOSOMES = [3]
table_base = args.tag + '__coverage_'

condition = ""# " WHERE pos = 16473265"

nrows = 0
with sqlite3.connect(args.inFile) as conn:
    curs = conn.cursor()
    for cc in CHROMOSOMES:
        try:
            curs.execute("SELECT Count(*) FROM " + table_base + condition)
            nrows += curs.fetchone()[0]
            print('chr %u, number of rows: %u' % (cc, nrows) )
        except:
            curs.execute('SELECT name FROM sqlite_master WHERE type = "table"')
            print('requested table has not been found!')
            print('available tables:')
            print(curs.fetchone())
            sys.exit()
    feature_array = np.zeros((nrows, 4*2*subrange))
    ii = 0
    
    for cc in CHROMOSOMES:
        curs.execute("select * from %s" % chrtable + condition)
        for rr in ResultIter(curs):
            # print(rr)
            cs = CoverageSqlite(cc, rr, args)
            # cs.print_range("snp_cov", subrange, h)
            feature_array[ii , :] = cs.feature_vector(subrange, "snp_cov")
            ii += 1
        
# X = feature_array - np.mean(feature_array, axis = 0)
from sklearn.decomposition import RandomizedPCA # import sklearn 
pca = RandomizedPCA(n_components=3, iterated_power = 7)
pca.fit(feature_array)                 
print(pca.explained_variance_ratio_) 
Y = pca.transform(feature_array)

from sklearn.decomposition import FastICA # import sklearn 
pca = FastICA(n_components=3)
pca.fit(feature_array)                 
# print(pca.explained_variance_ratio_) 
Y = pca.transform(feature_array)


import matplotlib.pyplot as plt

fig = plt.figure(figsize=(5,5))
plt.plot(Y[:,0], Y[:,1], 'b.')
plt.show()

fig = plt.figure(figsize=(5,5))
plt.plot(Y[:,0], Y[:,2], 'r.')
plt.show()


fig = plt.figure(figsize=(5,5))
plt.plot(Y[:,1], Y[:,2], 'g.')
plt.show()


        