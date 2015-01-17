# -*- coding: utf-8 -*-
"""
Created on Tue Dec 16 15:35:34 2014

@author: dima
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import struct

def printBlob(R):
    i = 0
    for c in range(len(R) ):
        print( '%2s' %hex( R[c] ) [2:], end = '  ')
        i += 1
        if (i % 16 == 0):
            print('', end = '\n')
    print('', end = '\n')
        
        
def controlPoints(rr):
   return [rr + rr*ii + ii for ii in range(4)]
   
def checkControlPoints(L, cps):
    for cc in cps:
        if not (L[cc]==0):
            printBlob(L)
            print('byte# %u, hex:  %2s' %(cc, hex(L[cc]) [2:] ) )
        assert(L[cc]==0), 'corrupted byte string!'

def blockToNpArray(a, sep = '\t'):
    # int.from_bytes(b'y\xcc\xa6\xbb', byteorder='big')
    return [ord(s)- int(33) for s in list(a)]
        # return np.array(a.strip(sep).split(sep)).astype(int)

def splitByteIntoNArrays(lineList, N, block_size):
    # assert( len(lineList) == (1 + 2*rangel)*N ), 'wrong blob length: %u' % len(lineList)
    " split into N int numpy sub-arrays within list intList"
    it = struct.iter_unpack('B', lineList)
    intList = [[] for Null in range(N)]
    for i, temp_char in enumerate(it):
        # print('ind %u\t arr %u' % (i, i//(2*rangel+1)) )
        intList[i // block_size].append(temp_char[0])
    return intList
    
def quadrupleArrayByte(lineList, rangel, sep = '\t'):
    block_size = 1 + 2*rangel
    assert( len(lineList) == block_size*4 ), 'wrong blob length: %u' % len(lineList)
    checkControlPoints(lineList, controlPoints(2*rangel))
    intList = splitByteIntoNArrays(lineList, 4, block_size)

    a = np.empty([2,2, rangel*2], dtype=int)
    a[Hi,F,:] = intList[0][:-1]
    a[Hi,R,:] = intList[1][:-1]
    a[Lo,F,:] = intList[2][:-1]
    a[Lo,R,:] = intList[3][:-1]
    return a


def quadrupleArrayChar(lineList, rangel, sep = '\t'):
    a = np.empty([2,2, rangel*2], dtype=int)
    ll = lineList.split('\t')
    a[Hi,F,:] = blockToNpArray(ll[0])
    a[Hi,R,:] = blockToNpArray(ll[1])
    a[Lo,F,:] = blockToNpArray(ll[2])
    a[Lo,R,:] = blockToNpArray(ll[3])
    return a
    
def quadrupleArray(lineList, ii, rangel, sep = '\t'):
    a = np.empty([2,2, rangel*2], dtype=int)
    a[Hi,F,:] = blockToNpArray(lineList[0 + 4*ii], sep)
    a[Hi,R,:] = blockToNpArray(lineList[1 + 4*ii], sep)
    a[Lo,F,:] = blockToNpArray(lineList[2 + 4*ii], sep)
    a[Lo,R,:] = blockToNpArray(lineList[3 + 4*ii], sep)
    return a

Lo = 0
Hi = 1
F = 0
R = 1

class PlotCoverage:
    def calcCentres(self):
        self.expCentreMut = np.sum(self.x * sum(self.aln_centres[Lo,:,:]) )/ np.sum(sum(self.aln_centres[Lo,:,:]))
        self.expCentre = np.sum(self.x * sum(self.aln_centres[Hi,:,:])) / np.sum( sum(self.aln_centres[Hi,:,:]) )
        if np.isnan(self.expCentreMut): 
            self.expCentreMut = 0 
        if np.isnan(self.expCentre) : 
            self.expCentre = 0
        
    def calcCountsRatio(self):
        self.totCounts = (self.tot_cov[Lo][F][rl]+ self.tot_cov[Lo][R][rl] + self.tot_cov[Hi][F][rl]+ self.tot_cov[Hi][R][rl])
        self.snp_ratio = (self.tot_cov[Lo][F][rl] + self.tot_cov[Lo][R][rl])/ \
        (self.tot_cov[Lo][F][rl]+ self.tot_cov[Lo][R][rl] + self.tot_cov[Hi][F][rl]+ self.tot_cov[Hi][R][rl])
        
    def plot_cov(self, field ):
        plt.plot(self.x, self.__dict__[field][Hi][F] + self.__dict__[field][Lo][F],  'b.-')
        plt.plot(self.x, -(self.__dict__[field][Hi][R] + self.__dict__[field][Lo][R]),  'b.-')
        plt.plot(self.x, self.__dict__[field][Hi][F], 'g.-')
        plt.plot(self.x, -self.__dict__[field][Hi][R], 'g.-')
        plt.xlim([-rl, rl]) 
        # plt.title('A tale of 2 subplots')
        plt.plot(self.x, self.__dict__[field][Lo][F], 'r.-')
        plt.plot(self.x, -self.__dict__[field][Lo][R], 'r.-')
        
    def plot_centres(self):
        markerline, stemlines, baseline = plt.stem(self.x, self.aln_centres[Hi][F], 'g.-')
        plt.setp(markerline, 'markerfacecolor', 'g')
        markerline, stemlines, baseline = plt.stem(self.x, -self.aln_centres[Hi][R], 'g.-')
        plt.setp(markerline, 'markerfacecolor', 'g')
        markerline, stemlines, baseline = plt.stem(self.x, self.aln_centres[Lo][F], 'r.-')
        plt.setp(markerline, 'markerfacecolor', 'r')
        markerline, stemlines, baseline = plt.stem(self.x, -self.aln_centres[Lo][R], 'r.-')
        plt.setp(markerline, 'markerfacecolor', 'r')
        plt.xlim([-rl, rl]) 
        
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
        self.plot_centres()
        plt.xlabel('position relative to the SNP, '+ \
        '\nratio: %3.3f, centre mut: %u, centre all %u, score = %3.3f,' % \
        (self.snp_ratio, self.expCentreMut, self.expCentre, self.score) )
        plt.ylabel('number of reads')

    def plot_scatters(self, field):
        self.fig = plt.figure(figsize=(10, 5))
        gs1 = gridspec.GridSpec(1, 2)
        ##########################################
        ax1 = plt.subplot(gs1[:,0], aspect = 'equal')
        xmax = max(max( self.__dict__[field][Hi][F] + self.__dict__[field][Hi][R]), \
        max(self.__dict__[field][Lo][F] + self.__dict__[field][Lo][R]))
        
        plt.plot( self.__dict__[field][Hi][F] + self.__dict__[field][Hi][R], \
        self.__dict__[field][Lo][F] + self.__dict__[field][Lo][R], 'bo', mec = 'b', mfc='none', label = "all")
        
        plt.plot( self.__dict__[field][Hi][F], self.__dict__[field][Lo][F], \
        'gx', label = "F")
        
        plt.plot( self.__dict__[field][Hi][R], self.__dict__[field][Lo][R], \
        'r+', label = "R")
        
        plt.xlabel('Hi (match)')
        plt.ylabel('Lo (mis-match)')
        
        ax1.set_xlim([0, xmax])
        ax1.set_ylim([0, xmax])
        plt.legend(loc = 2)
        ##########################################
        ax2 = plt.subplot(gs1[:,1], aspect = 'equal')
        xmax = max(max( self.__dict__[field][Hi][F] + self.__dict__[field][Lo][F]), \
        max(self.__dict__[field][Hi][R] + self.__dict__[field][Lo][R]))
        
        plt.scatter( self.__dict__[field][Hi][F] + self.__dict__[field][Lo][F], \
        self.__dict__[field][Hi][R] + self.__dict__[field][Lo][R], \
        marker ='o', edgecolor = 'b', c = 'b', facecolors='none',  label = "all")
        
        plt.scatter( self.__dict__[field][Lo][F], self.__dict__[field][Lo][R], \
        marker ='x', c = 'g', label = "mis-")
        
        plt.scatter( self.__dict__[field][Hi][F], self.__dict__[field][Hi][R],  \
        marker ='+', c = 'red', label = "match")

        plt.xlabel('R')
        plt.ylabel('F')
        ax2.set_xlim([0, xmax])
        ax2.set_ylim([0, xmax])
        plt.legend(loc = 4)
        # plt.legend(bbox_to_anchor=(xmin, ymax - 0.1*(ymax -ymin)), bbox_transform=plt.gcf().transFigure)
        
    def print_plot(self, out_dir, name):
        plt.savefig(os.path.join(out_dir, name + '_%u_%u.eps'% (self.chromosome, self.pos) ) )
        plt.close(self.fig) # 
        
    def show_plot(self):
        plt.show(self.fig)
    
    def feature_vector(self, field, subrange):
        fv = np.empty(4*2*subrange)
        fv[0:2*subrange] = self.__dict__[field][Hi][F][rl-subrange:rl+subrange]
        fv[2*subrange: 4*subrange] = self.__dict__[field][Hi][R][rl-subrange:rl+subrange]
        fv[4*subrange: 6*subrange] = self.__dict__[field][Lo][F][rl-subrange:rl+subrange] 
        fv[6*subrange: 8*subrange] = self.__dict__[field][Lo][R][rl-subrange:rl+subrange]
        return fv        
        
    def print_range(self, field, subrange, fid, sep = '\t'):
        m = sep.join( map(str, self.__dict__[field][Hi][F][rl-subrange:rl+subrange]) )
        print( m , end = sep, file = fid)
        m = sep.join( map(str, self.__dict__[field][Hi][R][rl-subrange:rl+subrange]) )
        print( m , end = sep, file = fid)
        m = sep.join( map(str, self.__dict__[field][Lo][F][rl-subrange:rl+subrange]) )
        print( m , end = sep, file = fid)
        m = sep.join( map(str, self.__dict__[field][Lo][R][rl-subrange:rl+subrange]) )
        print( m , end = '\n', file = fid)
        
###############################################################################    

class CoverageSummary( PlotCoverage ):
    
    def __init__( self, line, args ):
        lineList = line.rstrip("\n").split(args.csvseparator) 
        line_begin = lineList[0:3]
        global rl, sep, out_dir
        rl = args.rangel
        sep = args.csvseparator
        self.chromosome = int(line_begin[0])
        self.pos = int(line_begin[1])
        self.score = float(line_begin[2])
        
        self.parseArrays(lineList[4:])
        self.x = np.arange(- rl, rl ).astype(int)
        
        self.calcCentres()
        self.calcCountsRatio()
        return None    
        
    def parseArrays(self,lineList):
        self.tot_cov = quadrupleArray(lineList, 0, rl, sep)
        self.aln_centres = quadrupleArray(lineList, 1, rl, sep)
        self.snp_cov = quadrupleArray(lineList, 2, rl, sep)
        ##
        full_range = len(self.snp_cov[Lo][F])
        assert ( full_range == len(self.tot_cov[Hi][F]) )
        global rl
        rl = int(np.fix(full_range/2))
###############################################################################    

class CoverageSqlite( PlotCoverage ):
    
    def __init__( self, chrs, line, args ):
        global rl, sep, out_dir
        rl = args.rangel
        sep = args.csvseparator
        self.chromosome = chrs
        self.pos = line[0]
        self.score = line[1]
        self.totCounts = line[2]
        self.snp_ratio = line[4]
        
        self.parseArrays(line)
        self.x = np.arange(- rl, rl ).astype(int)
        
        self.calcCentres()
        return None    
        
    def parseArrays(self,lineList):
        firstArrInd = 5
        self.tot_cov =     quadrupleArrayByte(lineList[firstArrInd    ], rl, sep)
        self.snp_cov =     quadrupleArrayByte(lineList[firstArrInd + 1], rl, sep)
        self.aln_centres = quadrupleArrayByte(lineList[firstArrInd + 2], rl, sep)
        ##
        full_range = len(self.snp_cov[Lo][F])
        assert ( full_range == len(self.tot_cov[Hi][F]) )
        global rl
        rl = int(np.fix(full_range/2))
        