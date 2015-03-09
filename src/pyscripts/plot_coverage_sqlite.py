# -*- coding: utf-8 -*-
"""
Created on Wed Jan 14 12:52:05 2015

@author: dima
"""

import subprocess
from ReadCoverage import *
import sqlite3
# import pdb
###############################################################################
import argparse
parser = argparse.ArgumentParser('''
''')
parser.add_argument("inFile",
help="input sqlite3 database file (.db)")
# ../../bin/vshape-1.0.0/vcf.db

parser.add_argument("outFile", nargs='?', type=str, default='',
help="output file (.csv); for stdout type '-' ")

parser.add_argument("-l", "--rangel", type=int, default= 75,
help="")

parser.add_argument("-t", "--tag", type=str, default="HL10_ngm_nm4",
help="") 
"ABD192_bwa_nm3"

parser.add_argument("-s", "--csvseparator", type=str, default= '\t',
help="separator for the input file")

parser.add_argument("-o", "--out_dir", type=str, default= '../../plots',
help="directory for output")

parser.add_argument("-q", "--type", type=str, default= 'H',
help="blob data type (H-- unsigned short int, i -- int, B -- char)")

args = parser.parse_args()
###############################################################################
rl = args.rangel

# num_lines = int(subprocess.check_output(['wc', '-l', args.inFile]).decode().split(' ')[0])
# print("lines: %u" % num_lines)

destination_dir = os.path.join(args.out_dir, args.tag)
import shutil

purge = False

if purge:
    print("purging the `plots` directory...")
    shutil.rmtree(destination_dir, ignore_errors=True)

os.makedirs(destination_dir, exist_ok=True)
###############################################################################

CHROMOSOMES = ["Chr1"]
table_base = args.tag + '__coverage'

condition = ""# " WHERE pos = 16473265"

with sqlite3.connect(args.inFile) as conn:
    curs = conn.cursor()
    for cc in CHROMOSOMES:
        condition = "WHERE contig == " % cc + condition
        curs.execute("SELECT Count(*) FROM " + table_base + condition)
        nrows = curs.fetchone()[0]
        print('contig: %s, number of rows: %u' % (cc, nrows) )
        curs.execute("select * from %s" % chrtable + condition)
        for rr in ResultIter(curs):
            # print(rr)
            cs = CoverageSqlite(cc, rr, args)
            cs.plot()
            cs.print_plot( destination_dir , 'coverage')
            # plt.close("all")
            # cs.plot_scatters('snp_cov')
            # cs.show_plot()
            # cs.print_plot(destination_dir, 'scatters')
            plt.close("all")
        
        
    
    