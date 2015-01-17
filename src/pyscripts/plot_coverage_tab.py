
import subprocess
from CoverageSummary import *

# import pdb
###############################################################################
import argparse
parser = argparse.ArgumentParser('''
''')
parser.add_argument("inFile",
help="input file (.csv) with the gene IDs given in the 'geneID' column")

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

    
# with open('a', 'w') as a, open('b', 'w') as b:
with open(args.inFile, 'r') as f:
    # header = next(f)
    for line in f:
        cs = CoverageSummary( line, args )
        
        print("len(y) %u" % cs.aln_centres.shape[2] )
        print("plotting chr %u, pos %u" % (cs.chromosome, cs.pos) )
        print("centre of mut reads: %3f" % (cs.expCentreMut) )
        print("centre of all reads: %3f" % (cs.expCentre) )
        cs.plot()
        cs.print_plot( args.out_dir, 'coverage')
        plt.close("all")
        
        cs.plot_scatters('snp_cov')
        # cs.show_plot()
        cs.print_plot( args.out_dir, 'scatters')
        plt.close("all")
        
        
        
