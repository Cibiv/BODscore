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

# subprocess.call(['wc', '-l', args.inFile])
num_lines = int(subprocess.check_output(['wc', '-l', args.inFile]).decode().split(' ')[0])

print("lines: %u" % num_lines)

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
            
CHROMOSOMES = [1]

with sqlite3.connect(args.inFile) as conn:
    curs = conn.cursor()
    for cc in CHROMOSOMES:
        curs.execute("SELECT Count(*) FROM coverage_%u"% cc )
        nrows = curs.fetchone()[0]
        print('chr %u, number of rows: %u' % (cc, nrows) )
        curs.execute("select * from %s" % "coverage_%u"% cc )
        for rr in ResultIter(curs):
            # print(rr)
            cs = CoverageSqlite(cc, rr, args)
            cs.plot()
            cs.print_plot( args.out_dir, 'coverage')
            plt.close("all")
            cs.plot_scatters('snp_cov')
            # cs.show_plot()
            cs.print_plot( args.out_dir, 'scatters')
            plt.close("all")
        
        
    
    