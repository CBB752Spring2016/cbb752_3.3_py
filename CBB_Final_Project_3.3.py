# -*- coding: utf-8 -*-
"""
Created on Thu May  5 11:12:30 2016

@author: kevin
"""

import argparse
import numpy as np
import scipy.stats
import operator
import re
import pandas as pd
pd.set_option('precision',15)
### This is one way to read in arguments in Python..
parser = argparse.ArgumentParser(description='Enrichment Score Calculator')
parser.add_argument('-gct', '--gct', help='GCT file name', required=True)
parser.add_argument('-cls', '--cls', help='CLS file name', required=True)
parser.add_argument('-gs', '--geneset' , help='Gene set', required=True)
args = parser.parse_args()

#args.gct = "all_aml_train.gct"
#args.cls = "all_aml_train.cls"
#args.geneset = "geneset.txt"
#---------------Read Data into table-----------------#
GCT = open(args.gct, "rb")
genes = []
for G in GCT:
    genes.append(G.split("\t"))

np.shape(genes)
for i in range(len(genes)):
    genes[i][39] = genes[i][39].replace("\r\n","")
    #genes[i][0] = re.sub("_[A-z]*","",genes[i][0])##we should not remove the AT stuff
#cleaning
genes = np.array(genes).T
genes = genes[:,range(2,np.shape(genes)[1])]
genes = np.delete(genes, (1), axis=0)
genes = genes.T

CLS = open(args.cls, "rb")
for C in CLS:
    target = C.split(" ")
target = map(int, target)
numZero = target.count(0)#ALL
numOne = target.count(1)#AML

AML = []
ALL = []
pVals = {}

for i in range(1,len(genes)):
    ALL = map(float, genes[i][range(1,numZero+1)])
    AML = map(float, genes[i][range(numZero+1,numZero+numOne+1)])
    pVals[genes[i][0]] = scipy.stats.ttest_ind(ALL, AML,equal_var=False)[1]

sorted_pVals = pd.DataFrame(pVals.items(),columns=['Name', 'P-Val'])
sorted_pVals = sorted_pVals.sort_values("P-Val")
sorted_pVals = sorted_pVals.reset_index(drop = True)


geneset = open(args.geneset, "r")
GeneSet2 = []
for GS in geneset:
    GS = GS.replace("\n","")
    GeneSet2.append(GS)
GeneSet2 = GeneSet2[2:len(GeneSet2)]

listOfGenesInBoth = list(set(GeneSet2).intersection(sorted_pVals['Name'].tolist()))
sorted_pVals['Flag'] = 0
#sorted_pVals[sorted_pVals.Name.isin(listOfGenesInBoth)]['Flag']
#Set the flag for the genes in both
sorted_pVals.loc[sorted_pVals.Name.isin(listOfGenesInBoth), 'Flag'] = 1

#get total sum of values
sumOfPVals = sorted_pVals.loc[sorted_pVals.Flag == 1, 'P-Val'].sum()
#sorted_pVals[sorted_pVals.Name == "L49219"]

#get number of hits
numHits = sorted_pVals.loc[sorted_pVals.Flag == 1, 'Flag'].sum()
#get non hits
nonHits = abs(len(genes)-numHits)
#get weighted
sorted_pVals['Weighted_P-Val'] = 0.0
sorted_pVals.loc[sorted_pVals.Flag == 1,'Weighted_P-Val'] = sorted_pVals.loc[sorted_pVals.Flag == 1, 'P-Val']/sumOfPVals
sorted_pVals['pHit'] = 0
sorted_pVals['pMiss'] = 0

formulaPMiss = 1/float(len(sorted_pVals) - numHits)
#formulaPMiss = 1/float(2759 - numHits)


sumWeightedPval = 0
WeightedPval = sorted_pVals['Weighted_P-Val'].tolist()
pHit = []
pMis = []
ES = []
for i in range(len(WeightedPval)):
    sumWeightedPval += WeightedPval[i]
    pHit.append(sumWeightedPval)
    pMis.append(formulaPMiss*i)
    #sorted_pVals['pHit'][i] = sumWeightedPval
    #sorted_pVals['pMiss'][i] = formulaPMiss*i

sorted_pVals['pHit'] = pHit
sorted_pVals['pMiss'] = pMis


ES = [abs(i - j) for i, j in zip(pHit, pMis)]
max(ES)

print ("The Enrichment Score is: " +  str(max(ES)))
#listTest = ["L05514","L49173","L49219","L49229","L76670","M19045","M23575","S78693",
#            "X14008","X59244","X97230","D16154"]
#
#for text in genes[:,0]:
#    for sList in listTest:
#        if sList in text:
#            print(text)
