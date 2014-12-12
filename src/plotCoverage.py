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
        
        self.expCentreMut = np.sum(self.x * self.alnmCentresMut) / np.sum(self.alnmCentresMut)
        self.expCentre = np.sum(self.x * self.alnmCentres) / np.sum(self.alnmCentres)
        self.diffQ = -np.subtract(self.highQ, self.lowQ)
        
        self.snp_ratio = self.diffQ[rl]/self.lowQ[rl]
        return None    
        
    def parseArrays(self,lineList):
        self.highQ = blockToNpArray(lineList[1])
        self.lowQ = blockToNpArray(lineList[2])
        self.alnmCentres = blockToNpArray(lineList[3])
        self.alnmCentresMut = blockToNpArray(lineList[4])
        self.locHighQ = blockToNpArray(lineList[5])
        self.locLowQ = blockToNpArray(lineList[6])
        
        full_range = len(self.locLowQ)
        global rl        
        print(rl)
        assert ( full_range == len(self.highQ) )
        rl = int(np.fix(full_range/2))
        
    
    def plot(self):
        fig = plt.figure(figsize=(8,10))
        gs1 = gridspec.GridSpec(5, 1)
        # gs1.update(left=0.05, right=0.48, wspace=0.05)        
  
        # plt.title('ratio: %3.3f, centre mut: %u, centre all %u, score = %3.3f,' %(self.snp_ratio, self.expCentreMut, self.expCentre, self.score) )
  
        ax1 = plt.subplot(gs1[:2, :]) # plt.subplot(211)
        plt.plot(self.x, self.lowQ,  'b.-')
        plt.plot(self.x, self.highQ, 'g.-')
        plt.xlim([-rl, rl]) 
        plt.xlabel('position relative to the SNP, ratio=%3f' % self.snp_ratio)
        plt.ylabel('coverage')
        # plt.title('A tale of 2 subplots')
        plt.plot(self.x, self.diffQ, 'r.-')
        
        ax2 = plt.subplot(gs1[2:4, :]) # plt.subplot(211)
        plt.plot(self.x, self.locLowQ,  'b.-')
        plt.plot(self.x, self.locHighQ, 'g.-')
        plt.xlim([-rl, rl]) 
        plt.xlabel('position relative to the SNP, score = %3f' % self.score )
        plt.ylabel('local coverage')
        plt.plot(self.x, self.locLowQ - self.locHighQ, 'r.-')        
        
        ax3 = plt.subplot(gs1[-1, :]) # plt.subplot(212)
        markerline, stemlines, baseline = plt.stem(self.x, self.alnmCentres, 'g.-')
        plt.setp(markerline, 'markerfacecolor', 'g')
        markerline, stemlines, baseline = plt.stem(self.x, self.alnmCentresMut, 'r.-')
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
        
        print("len(y) %u" % len(cs.alnmCentres) )
        print("plotting chr %u, pos %u" % (cs.chromosome, cs.pos) )
        print("centre of mut reads: %3f" % (cs.expCentreMut) )
        print("centre of all reads: %3f" % (cs.expCentre) )
        cs.plot()        
        
        
