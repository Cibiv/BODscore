#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 14 12:52:05 2015

@author: dima
"""

import subprocess
from CoverageSummary import *
import sqlite3
# import pdb
###############################################################################
import argparse
parser = argparse.ArgumentParser('''
''')
parser.add_argument("inFile",
help="input sqlite3 database file (.db)")

parser.add_argument("outFile", nargs='?', type=str, default='',
help="output file (.csv); for stdout type '-' ")

parser.add_argument("-l", "--rangel", type=int, default="75",
help="")

parser.add_argument("-s", "--csvseparator", type=str, default= '\t',
help="separator for the input file")

parser.add_argument("-o", "--out_dir", type=str, default= './plots',
help="directory for output")

args = parser.parse_args()
###############################################################################
rl = args.rangel

import shutil
shutil.rmtree(args.out_dir)

os.mkdir(args.out_dir)
###############################################################################
def ResultIter(cursor, arraysize=5000):
    'An iterator that uses fetchmany to keep memory usage down'
    while True:
        results = cursor.fetchmany(arraysize)
        if not results:
            break
        for result in results:
            yield result
###############################################################################
            
condition = " WHERE (snp_ratio < .67) AND (snp_ratio > .33)"
CHROMOSOMES = [1, 4, 5]
nrows = 0
with sqlite3.connect(args.inFile) as conn:
    curs = conn.cursor()
    for cc in CHROMOSOMES:
        ii = 0
        curs.execute("SELECT Count(*) FROM coverage_%u"% cc + condition)
        nrows_chr = curs.fetchone()[0]
        print('chr %u, number of rows: %u' % (cc, nrows_chr) )
        nrows += nrows_chr
        

loc_covs = np.zeros( (nrows, 1), dtype = int)
snp_rats = np.zeros( (nrows, 1), dtype = int)
tot_covs = np.zeros( (nrows, 4*2*rl), dtype = int)

ii = 0
with sqlite3.connect(args.inFile) as conn:
    curs = conn.cursor()
    for cc in CHROMOSOMES:
        ##
        curs.execute("select * from %s" % "coverage_%u"% cc + condition)
        for rr in ResultIter(curs):
            # print('chr:%u, pos:%u' % (cc, rr[0]) ) 
            cs = CoverageSqlite(cc, rr, args)
            loc_covs[ii] = cs.totCounts
            snp_rats[ii] = cs.snp_ratio
            tot_covs[ii,:] = cs.tot_cov.flatten()
            ii+=1

tot_covs.shape
    
    