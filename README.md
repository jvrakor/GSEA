# GSEA method 
as described in original article available at http://www.pnas.org/content/102/43/15545.abstract

Inputs: leukemia.txt (gene expression profiles for two groups of leukemia patients), 
 pathways.txt (list of gene sets that describe the metabolic pathways), 
 N (number of permutations used for computation of nominal p-values and normalized enrichment scores, default N=1000) 

Outputs: gsea.csv (A table of pathway names, number of genes, enrichment scores, normalized enrichment scores and nominal p-values, sorted by normalized enrichment scores)
 
Note that p=1 is used in analysis.  
 
gsea.py - original version, no outputs, slow and incorrect 

gsea_ver1.py - improved version, with an output and better performance, still incorrect

gsea_ver2.py - number of genes added to the output, method slightly differed
