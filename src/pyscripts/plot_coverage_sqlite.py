# -*- coding: utf-8 -*-
"""
Created on Wed Jan 14 12:52:05 2015

@author: dima
"""
from __future__ import print_function
import subprocess
from ReadCoverage import *
import sqlite3
import sys
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

parser.add_argument("-r", "--range", type=int, default= 75,
help="")

parser.add_argument("-t", "--tag", type=str, default="HL7_ngm_nm4",
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
rl = args.range

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
def dict_factory(cursor, row):
    d = {}
    for idx, col in enumerate(cursor.description):
        d[col[0]] = row[idx]
    return d
###############################################################################


CHROMOSOMES = ["Chr1"]
table_base = args.tag + '__coverage'

condition0 = ""# " WHERE pos = 16473265"

with sqlite3.connect(args.inFile) as conn:
    conn.row_factory =  sqlite3.Row
    curs = conn.cursor()
    try:
        query = "SELECT Count(*) FROM " + table_base + condition0
        curs.execute( query )
        nrows = curs.fetchone()[0]
        print( "query: " , query )
        print('total number of rows: %u' % ( nrows) )
        curs.execute('PRAGMA table_info(%s)'% table_base)
        while 1:
            col_names = curs.fetchone()
            if isinstance(col_names, type(None)):
                break
            else:
                print('\t', col_names[1])
    except:
        curs.execute('SELECT name FROM sqlite_master WHERE type = "table"')
        print('requested table has not been found!')
        print('available tables:')
        tables = curs.fetchall()
        for tt in tables:
            print('\t', tt['name'])
        raise
        sys.exit()
    query_contigs = 'SELECT %(clmn)s, COUNT( %(clmn)s ) AS count FROM "%(tbl)s"  GROUP BY %(clmn)s' % \
                {'clmn': "contig", 'tbl': table_base }
    
    curs.execute(query_contigs)
    chromosomes = curs.fetchall()
    print( list(chromosomes[0]) )
    print(chromosomes[0]["contig"])

    for co in chromosomes:
        cc = co["contig"]
        if len(cc) == 0:
            RuntimeWarning('empty contig. Skipping...')
            continue
        condition = " WHERE contig == '%s'" % cc + condition0
        query = "SELECT Count(*) FROM " + table_base + condition
        print( query, file = sys.stderr )
        curs.execute( query )
        nrows = curs.fetchone()[0]
        # print( co.keys() )
        print('contig: %s, number of rows: %u' % (cc, co["count"] ) )
        curs.execute("select * from %s" % table_base + condition)
        for rr in ResultIter(curs):
            # print(rr)
            cs = CoverageSqlite(rr, args)
            cs.plot()
            cs.print_plot( destination_dir , 'coverage')
            # plt.close("all")
            # cs.plot_scatters('snp_cov')
            # cs.show_plot()
            # cs.print_plot(destination_dir, 'scatters')
            plt.close("all")

