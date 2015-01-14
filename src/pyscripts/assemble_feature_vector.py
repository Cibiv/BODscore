# -*- coding: utf-8 -*-
"""
Created on Tue Dec 16 15:51:09 2014

@author: dima
"""

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

parser.add_argument("-o", "--outFile", type=str, default= 'out.tab',
help="file for output")

args = parser.parse_args()
###############################################################################
global rl
rl = args.rangel
subrange = 50;
# subprocess.call(['wc', '-l', args.inFile])
num_lines = int(subprocess.check_output(['wc', '-l', args.inFile]).decode().split(' ')[0])
print("lines: %u" % num_lines)

feature_array = np.zeros((num_lines, 4*2*subrange))

ii = 0   
# with open('a', 'w') as a, open('b', 'w') as b:
with open(args.inFile, 'r') as f, open(args.outFile, 'w') as h:
    # header = next(f)
    for line in f:
        cs = CoverageSummary( line, args )
        # cs.print_range("snp_cov", subrange, h)
        feature_array[ii , :] = cs.feature_vector("snp_cov" , subrange)
        ii += 1
        
# X = feature_array - np.mean(feature_array, axis = 0)
from sklearn.decomposition import RandomizedPCA # import sklearn 
pca = RandomizedPCA(n_components=3, iterated_power = 7)
pca.fit(feature_array)                 
print(pca.explained_variance_ratio_) 
Y = pca.transform(feature_array)

from sklearn.decomposition import FastICA # import sklearn 
pca = FastICA(n_components=3)
pca.fit(feature_array)                 
# print(pca.explained_variance_ratio_) 
Y = pca.transform(feature_array)


import matplotlib.pyplot as plt

fig = plt.figure(figsize=(5,5))
plt.plot(Y[:,0], Y[:,1], 'b.')
plt.show()

fig = plt.figure(figsize=(5,5))
plt.plot(Y[:,0], Y[:,2], 'r.')
plt.show()


fig = plt.figure(figsize=(5,5))
plt.plot(Y[:,1], Y[:,2], 'g.')
plt.show()


        