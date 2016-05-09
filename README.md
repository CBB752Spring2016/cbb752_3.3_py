# Calculate Enrichment Score of two gene sets

##Input:

Gene set from microarray (usually a microarray experiment)

Files:

1.a GCT file (Microarray Experiment [The experiment matrix])

2.a CLS file (Microarray Experiment Sample Outcomes [The target vector])

3.a TXT file (the other gene set in txt format [all genes must be seperated by new lines "\n" and must have two lines that will be skipped before the first gene similar to the msigdb gene sets]) 

```{python}
 python CBB_Final_Project_3.3.py -gct all_aml_train.preprocessed.gct -cls all_aml_train.cls -gs geneset.txt 

```

###Input Examples:
####GCT file
```
#1.2																																							
7129	38																																						
Name	Description	ALL_19769_B-cell	ALL_23953_B-cell	ALL_28373_B-cell ...
AFFX-BioB-5_at	AFFX-BioB-5_at (endogenous control)	-214	-135	-106 ...
AFFX-BioB-M_at	AFFX-BioB-M_at (endogenous control)	-153	-114	-125 ...
AFFX-BioB-3_at	AFFX-BioB-3_at (endogenous control)	-58	265	-76 ...
AFFX-BioC-5_at	AFFX-BioC-5_at (endogenous control)	88	12	168 ...
...
...
...
```
####CLS file:
```
38 2 1
# ALL AML
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1
```

####Gene Set
```
HALLMARK_ADIPOGENESIS
> Genes up-regulated during adipocyte differentiation (adipogenesis).
ABCA1
ABCB8
ACAA2
ACADL
ACADM
ACADS
ACLY
ACO2
ACOX1
ADCY6
ADIG
ADIPOQ
...
```
##Output:

Print out of the Enrichment Score in the terminal. Example:
```
The Enrichment Score is: 0.474532816974
```

#Explanation of algorithm

We have a gene set G which we use as an "index" to the genes (essentially map G to the genes in the expression data [to the same genes]). We then perform a t-test based on the Y (target vector) of the samples. This will yield p-values and we increment (or decrement) them based on the correlation of the gene and the phenotype (the Y).

Steps:
1. Take the input matrix, and the target vector and perform a t-test on each gene based on the row (the two popluation are the two classes from the CLS file).
2. Sort the p-values
3. Intersect the gene set (from the txt file) with the sorted list
4. Incremenent the value of the gene based on the correleation between the gene and outcome
 
Algorithm formulation:
![Algorithm]
(http://i.imgur.com/aqYa4SV.png)

References:

1.[Gene set enrichment analysis: A knowledge-based
approach for interpreting genome-wide
expression profiles](http://www.pnas.org/content/102/43/15545.full.pdf): Aravind Subramaniana, Pablo Tamayoa, Vamsi K. Moothaa, Sayan Mukherjeed, Benjamin L. Eberta, Michael A. Gillettea, Amanda Paulovichg, Scott L. Pomeroyh, Todd R. Goluba, Eric S. Landera, and Jill P. Mesirova,

2.[Gene set analysis methods: statistical models and methodological differences](http://bib.oxfordjournals.org/content/15/4/504.full): Henryk Maciejewski, Institute of Computer Engineering, Control and Robotics, Wroclaw University of Technology

###NOTE:
The genes in the gene set and the genes in the GCT file must be in the same format, so if the genes have some extension such as _at or _st then your gene set must have this otherwise it will not match up!
