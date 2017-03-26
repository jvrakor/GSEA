# GSEA method 
as described in original article available at http://www.pnas.org/content/102/43/15545.abstract

Inputs: leukemia.txt (gene expression profiles for two groups of leukemia patients)
         pathways.txt (list of gene sets that describe the metabolic pathways) 
         N (number of permutations used for computation of nominal p-values and normalized enrichment scores, default N=1000)

Outputs: gsea.csv (A table of pathway names, enrichment scores, normalized enrichment scores and nominal p-values, sorted by normalized enrichment scores)
 
 Note that p=1 is used in analysis.  
