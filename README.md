# Calculate Enrichment Score of two gene sets

##Input:

Gene set from microarray (usually a microarray experiment)

Files:

1.a GCT file (Microarray Experiment [The experiment matrix])
2.a CLS file (Microarray Experiment Sample Outcomes [The target vector])
3.a TXT file (the other gene set in txt format [all genes must be seperated by new lines "\n" and must have two lines that will be skipped before the first gene similar to the msigdb gene sets]) 

##Output:

Printout of the Enrichment Score in the terminal

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
