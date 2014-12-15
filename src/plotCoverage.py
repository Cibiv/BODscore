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

Lo = 0
Hi = 1
F = 0
R = 1
    
def blockToNpArray(a):
        return np.array(a.strip(args.csvseparator).split(args.csvseparator)).astype(int)
    
def quadrupleArray(lineList, ii):
    a = np.empty([2,2,args.rangel*2])
    a[Hi,F,:] = blockToNpArray(lineList[1 + 4*ii])
    a[Hi,R,:] = blockToNpArray(lineList[2 + 4*ii])
    a[Lo,F,:] = blockToNpArray(lineList[3 + 4*ii])
    a[Lo,R,:] = blockToNpArray(lineList[4 + 4*ii])
    return a
    
class CoverageSummary():
    
    def __init__( self, line ):
        lineList = line.rstrip("\n").split("|")       
        line_begin = lineList[0].split(args.csvseparator)
        self.chromosome = int(line_begin[0])
        self.pos = int(line_begin[1])
        self.score = float(line_begin[2])
        
        self.parseArrays(lineList)
        self.x = np.arange(- rl, rl ).astype(int)
        
        self.expCentreMut = np.sum(self.x * sum(self.aln_centres[Lo,:,:]) )/ np.sum(sum(self.aln_centres[Lo,:,:]))
        self.expCentre = np.sum(self.x * sum(self.aln_centres[Hi,:,:])) / np.sum( sum(self.aln_centres[Hi,:,:]) )
        if np.isnan(self.expCentreMut): 
            self.expCentreMut = 0        
        if np.isnan(self.expCentre) : 
            self.expCentre = 0
        
        self.snp_ratio = (self.tot_cov[Lo][F][rl] + self.tot_cov[Lo][R][rl])/ \
        (self.tot_cov[Lo][F][rl]+ self.tot_cov[Lo][R][rl] + self.tot_cov[Hi][F][rl]+ self.tot_cov[Hi][R][rl])
        return None    
        
        
    def parseArrays(self,lineList):
        self.tot_cov = quadrupleArray(lineList, 0)
        self.aln_centres = quadrupleArray(lineList, 1)
        self.snp_cov = quadrupleArray(lineList, 2)
        
        full_range = len(self.snp_cov[Lo][F])
        global rl        
        print(rl)
        assert ( full_range == len(self.tot_cov[Hi][F]) )
        rl = int(np.fix(full_range/2))
        
    def plot_cov(self, field ):
        plt.plot(self.x, self.__dict__[field][Hi][F] + self.__dict__[field][Lo][F],  'b.-')
        plt.plot(self.x, -(self.__dict__[field][Hi][R] + self.__dict__[field][Lo][R]),  'b.-')
        plt.plot(self.x, self.__dict__[field][Hi][F], 'g.-')
        plt.plot(self.x, -self.__dict__[field][Hi][R], 'g.-')
        plt.xlim([-rl, rl]) 
        # plt.title('A tale of 2 subplots')
        plt.plot(self.x, self.__dict__[field][Lo][F], 'r.-')
        plt.plot(self.x, -self.__dict__[field][Lo][R], 'r.-')
    
    def plot(self):
        self.fig = plt.figure(figsize=(8,10))
        gs1 = gridspec.GridSpec(5, 1)
        # gs1.update(left=0.05, right=0.48, wspace=0.05)        
  
        # plt.title('ratio: %3.3f, centre mut: %u, centre all %u, score = %3.3f,' %(self.snp_ratio, self.expCentreMut, self.expCentre, self.score) )
  
        ax1 = plt.subplot(gs1[:2, :]) # plt.subplot(211)
        self.plot_cov( "tot_cov" )
        plt.ylabel('coverage')    
        plt.xlabel('position relative to the SNP, ratio=%3f' % self.snp_ratio)
        
        
        ax2 = plt.subplot(gs1[2:4, :]) # plt.subplot(211)
        self.plot_cov( "snp_cov" )
        plt.xlabel('position relative to the SNP, score = %3f' % self.score )
        plt.ylabel('local coverage')
        
        
        ax3 = plt.subplot(gs1[-1, :]) # plt.subplot(212)
        markerline, stemlines, baseline = plt.stem(self.x, self.aln_centres[Hi][F], 'g.-')
        plt.setp(markerline, 'markerfacecolor', 'g')
        markerline, stemlines, baseline = plt.stem(self.x, -self.aln_centres[Hi][R], 'g.-')
        plt.setp(markerline, 'markerfacecolor', 'g')
        markerline, stemlines, baseline = plt.stem(self.x, self.aln_centres[Lo][F], 'r.-')
        plt.setp(markerline, 'markerfacecolor', 'r')
        markerline, stemlines, baseline = plt.stem(self.x, -self.aln_centres[Lo][R], 'r.-')
        plt.setp(markerline, 'markerfacecolor', 'r')
        plt.xlim([-rl, rl]) 
        plt.xlabel('position relative to the SNP, '+ \
        '\nratio: %3.3f, centre mut: %u, centre all %u, score = %3.3f,' % \
        (self.snp_ratio, self.expCentreMut, self.expCentre, self.score) )
        plt.ylabel('number of reads')


    def print_plot(self):
        plt.savefig(os.path.join(args.out_dir, 'coverage%u_%u.eps'% (self.chromosome, self.pos) ) )
        plt.close(self.fig) # 
        
    def show_plot(self):
        plt.show(self.fig)

with open(args.inFile) as f:
    # header = next(f)
    for line in f:
        cs = CoverageSummary( line )
        
        print("len(y) %u" % cs.aln_centres.shape[2] )
        print("plotting chr %u, pos %u" % (cs.chromosome, cs.pos) )
        print("centre of mut reads: %3f" % (cs.expCentreMut) )
        print("centre of all reads: %3f" % (cs.expCentre) )
        cs.plot()
        cs.print_plot()
        plt.close("all")
        # cs.show_plot()
        
        
