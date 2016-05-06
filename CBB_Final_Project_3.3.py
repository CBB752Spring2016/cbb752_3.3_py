# -*- coding: utf-8 -*-
"""
Created on Thu May  5 11:12:30 2016

@author: kevin
"""
### This is one way to read in arguments in Python. We need to read input file and score file.
parser = argparse.ArgumentParser(description='Enrichment Score Calculator')
parser.add_argument('-gct', '--gct', help='GCT file name', required=True)
parser.add_argument('-cls', '--cls', help='CLS file name', required=True)
parser.add_argument('-gs', '--geneset' , help='Gene set', required=True)
args = parser.parse_args()

import numpy as np
import scipy.stats
import operator
#---------------Read Data into table-----------------#
GCT = open(parser.gct, "rb")
genes = []
for G in GCT:
    genes.append(G.split("\t"))

np.shape(genes)
for i in range(len(genes)):
    genes[i][39] = genes[i][39].replace("\r\n","")
#cleaning
genes = np.array(genes).T
genes = genes[:,range(2,np.shape(genes)[1])]
genes = np.delete(genes, (1), axis=0)
genes = genes.T

CLS = open(parser.cls, "rb")
for C in CLS:
    target = C.split(" ")
target = map(int, target)
numZero = target.count(0)#ALL
numOne = target.count(1)#AML

AML = []
ALL = []
pVals = {}

for i in range(1,len(genes)):
    ALL = map(int, genes[i][range(1,numZero+1)])
    AML = map(int, genes[i][range(numZero+1,numZero+numOne+1)])
    pVals[genes[i][0]] = scipy.stats.ttest_ind(ALL, AML)[1]

sorted_pVals = sorted(pVals.items(), key=operator.itemgetter(1))

geneset = open(parser.geneset, "r")
GeneSet2 = []
for GS in geneset:
    GS = GS.replace("\n","")
    GeneSet2.append(GS)
GeneSet2 = GeneSet2[2:len(GeneSet2)]


sumOfPVals = 0
listOfGenesInBoth = []
numHits = 0
for gene in GeneSet2:
    for j in range(len(sorted_pVals)):
        if gene in sorted_pVals[i][0]:
            listOfGenesInBoth.append(sorted_pVals[i][0])
            sumOfPVals += sorted_pVals[i][1]
            numHits += 1

Phit = 0
for gene in GeneSet2:
    for j in range(len(sorted_pVals)):
        if gene in sorted_pVals[i][0]:
            Phit += sorted_pVals[i][0]/sumOfPVals
nonHits = abs(len(genes)-numHits)

Pmiss = 0
for i in range(len(genes)):
    Pmiss += 1/nonHits

ES = Phit - Pmiss

print ("The Enrichment Score is: " +  str(ES))
