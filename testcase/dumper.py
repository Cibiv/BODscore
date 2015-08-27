import sys, os
##############################################################
import subprocess
def runcmd(x):
    assert (type(x) == str), 'not a string!'
    p = subprocess.Popen(x, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    return (out.decode('utf-8'), err.decode('utf-8'))
##############################################################
import argparse
parser = argparse.ArgumentParser('adds a vcf file to an SQLite database')

parser.add_argument("dbpath",
                    help="a path to the database")

parser.add_argument( "csvpath", nargs = '?', default = '',
                     help="path to the output `csv` file")

parser.add_argument( "-t", "--tablename",  default = '', type=str,
                   help="table name")

parser.add_argument( "-s", "--separator",  default = ';',
                    help="csv separator")

args = parser.parse_args()
##############################################################

if not len( args.tablename ):
    tbl = os.path.basename( sys.argv[1] ).rstrip(".db")
else:
    tbl = args.tablename

if len(args.csvpath):
    outfile = args.csvpath
else:
    outfile = tbl + ".tab"

sqlc = """.mode csv
-- use '.separator SOME_STRING' for something other than a comma.
.separator %(sep)s
.headers on
.out %(outfile)s
select chr, pos, counts0+counts1 as counts from %(tbl)s;
"""

sqlc  = sqlc % {"tbl" : tbl, "outfile" : outfile, "sep": args.separator}

sqlcmdfile = "_dumpsum_%s.sql" % tbl
with open( sqlcmdfile , "w") as text_file:
    text_file.write( sqlc )

print("command: " )
print(sqlc )

out = runcmd('sqlite3 %(infile)s <  %(cmd)s ' % {"infile": args.dbpath, "cmd": sqlcmdfile }  )

print(out)

os.remove( sqlcmdfile )
