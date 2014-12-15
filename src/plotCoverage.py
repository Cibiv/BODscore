import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import os
import subprocess

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

def blockToNpArray(a):
        return np.array(a.strip(args.csvseparator).split(args.csvseparator)).astype(int)
    
class CoverageSummary():
    
    def __init__( self, line ):
        lineList = line.rstrip("\n").split("|")       
        line_begin = lineList[0].split(args.csvseparator)
        self.chromosome = int(line_begin[0])
        self.pos = int(line_begin[1])
        self.score = float(line_begin[2])
        
        self.parseArrays(lineList)
        self.x = np.arange(- rl, rl ).astype(int)
        
        self.expCentreMut = np.sum(self.x * self.snp_alnCentres) / np.sum(self.snp_alnCentres)
        self.expCentre = np.sum(self.x * self.tot_alnCentres) / np.sum(self.tot_alnCentres)
        
        
        self.snp_ratio = (self.covLoF[rl] + self.covLoR[rl])/ \
        (self.covLoF[rl]+ self.covLoR[rl] + self.covHiF[rl]+ self.covHiR[rl])
        return None    
        
    def parseArrays(self,lineList):
        self.covHiF = blockToNpArray(lineList[1])
        self.covHiR = blockToNpArray(lineList[2])
        self.covLoF = blockToNpArray(lineList[3])
        self.covLoR = blockToNpArray(lineList[4])
        self.tot_alnCentres = blockToNpArray(lineList[5])
        self.snp_alnCentres = blockToNpArray(lineList[6])
        self.snp_covHiF = blockToNpArray(lineList[7])
        self.snp_covHiR = blockToNpArray(lineList[8])
        self.snp_covLoF = blockToNpArray(lineList[9])
        self.snp_covLoR = blockToNpArray(lineList[10])
        
        full_range = len(self.snp_covLoF)
        global rl        
        print(rl)
        assert ( full_range == len(self.covHiF) )
        rl = int(np.fix(full_range/2))
        
    
    def plot(self):
        fig = plt.figure(figsize=(8,10))
        gs1 = gridspec.GridSpec(5, 1)
        # gs1.update(left=0.05, right=0.48, wspace=0.05)        
  
        # plt.title('ratio: %3.3f, centre mut: %u, centre all %u, score = %3.3f,' %(self.snp_ratio, self.expCentreMut, self.expCentre, self.score) )
  
        ax1 = plt.subplot(gs1[:2, :]) # plt.subplot(211)
        plt.plot(self.x, self.covHiF + self.covLoF,  'b.-')
        plt.plot(self.x, -(self.covHiR + self.covLoR),  'b.-')
        plt.plot(self.x, self.covHiF, 'g.-')
        plt.plot(self.x, -self.covHiR, 'g.-')
        plt.xlim([-rl, rl]) 
        plt.xlabel('position relative to the SNP, ratio=%3f' % self.snp_ratio)
        plt.ylabel('coverage')
        # plt.title('A tale of 2 subplots')
        plt.plot(self.x, self.covLoF, 'r.-')
        plt.plot(self.x, -self.covLoR, 'r.-')
        
        ax2 = plt.subplot(gs1[2:4, :]) # plt.subplot(211)
        plt.plot(self.x, self.snp_covLoF + self.snp_covHiF,  'b.-')
        plt.plot(self.x, - (self.snp_covLoR + self.snp_covHiR),  'b.-')
        plt.plot(self.x, self.snp_covHiF, 'g.-')
        plt.plot(self.x, - self.snp_covHiR, 'g.-')
        plt.xlim([-rl, rl]) 
        plt.xlabel('position relative to the SNP, score = %3f' % self.score )
        plt.ylabel('local coverage')
        plt.plot(self.x, self.snp_covLoF , 'r.-')
        plt.plot(self.x, -self.snp_covLoR , 'r.-')        
        
        ax3 = plt.subplot(gs1[-1, :]) # plt.subplot(212)
        markerline, stemlines, baseline = plt.stem(self.x, self.tot_alnCentres, 'g.-')
        plt.setp(markerline, 'markerfacecolor', 'g')
        markerline, stemlines, baseline = plt.stem(self.x, self.snp_alnCentres, 'r.-')
        plt.setp(markerline, 'markerfacecolor', 'r')        
        plt.xlim([-rl, rl]) 
        plt.xlabel('position relative to the SNP, \nratio: %3.3f, centre mut: %u, centre all %u, score = %3.3f,' %(self.snp_ratio, self.expCentreMut, self.expCentre, self.score) )
        plt.ylabel('number of reads')
        # plt.show()
        plt.savefig(os.path.join(args.out_dir, 'coverage%u_%u.eps'% (self.chromosome, self.pos) ) )
        plt.close() # 
        

with open(args.inFile) as f:
    # header = next(f)
    for line in f:
        cs = CoverageSummary( line )
        
        print("len(y) %u" % len(cs.tot_alnCentres) )
        print("plotting chr %u, pos %u" % (cs.chromosome, cs.pos) )
        print("centre of mut reads: %3f" % (cs.expCentreMut) )
        print("centre of all reads: %3f" % (cs.expCentre) )
        cs.plot()        
        
        
