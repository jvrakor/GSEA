#this is a Python GSEA implementation. Imput files are leukemia.txt and pathways.txt. The implementation follows the original article available at
#(http://www.pnas.org/content/102/43/15545.abstract). p is set to p=1 and thus ommited. There are no outputs yet.

import numpy as np

PH1=[] #reserved for ALL phemotyes data
PH2=[] #reserved for AML phenotypes data
GENE=[] #reserved for gene names from leukemia.txt
N=100 #number of permutations


with open('leukemia.txt', 'r') as f:
	header = np.array(f.readline().replace("\r\n","").split('\t')) #first line, phenotype names
	lines = f.readlines()

for line in lines:
	PHI=[]
	PHII=[]
	for i, val in enumerate(np.array(line.replace("\r\n","").split('\t'))): 
		if header[i]=='ALL': #phenotype ALL
			PHI.append(float(val))
		elif header[i]=='AML': #phenotype AML
			PHII.append(float(val))	 
		else:
			GENE.append(val)	
	PH1.append(PHI)
	PH2.append(PHII)

f.close()	

with open('pathways.txt', 'r') as f:
	lines = f.readlines()
	
PathwayGenes=[] #reserved for gene names 
PathwayName=[] #reserved for pathway names

for line in lines:
	PathwayName.append(line.split('\t')[0])
	if line.split('\t')[1]=='BLACK':  #name of the pathway is always in the first column, the second column is usually a description of the pathway but sometimes a gene BLACK
		PathwayGenes.append(line.replace("\n","").split('\t')[1:])
	else:
		PathwayGenes.append(line.replace("\n","").split('\t')[2:]) 

f.close()

Lph1=len(PH1[0]) #number of ALL phenotype samples
Lph2=len(PH2[0]) #number of AML phenotype samples

PH=np.hstack((PH1,PH2)) #combining ALL and AML phenotype values into a single list

#COMPUTING ES(S) and p-scores, first row in ESNull corresponds to original data and all the rest are permutations, each column corresponds to a different pathway   
ESNull=[]
ESS=[] #coresponds to original data ES(S) values
for n in range(N+1):
	r=[] #reserved for r values
	ES=[] #reserved for ES(S) scores
	
	#assigning random original phenotype names to samples, but permuting values instead of names
	if n!=0:
		PH=PH[:, np.random.permutation(PH.shape[1])]
		PH1=PH[:,:Lph1] #the first portion of a combined array, up to the number of ALL phenotype samples
		PH2=PH[:,Lph1:] #the second portion, corresponding to AML phenotype

	#construction of correlation values array r, Signal2noise method is used
	for (val1, val2) in zip(PH1,PH2):
		r.append((np.mean(val1)-np.mean(val2))/(np.std(val1)+np.std(val2)))

	#sorting of gene names and r according to r value and resorting them so we start with highest correlation first
	GENE=sorted(zip(np.argsort(r),GENE))
	GENE=np.array(GENE)[:,1][::-1]
	r=np.sort(r)[::-1]
	
	#COMPUTING ES(S) for each pathway S
	for pathway in PathwayGenes:
		phit=0	
		pmiss=0
		Phit=[]
		Pmiss=[]
		pathway=set(pathway)
		for i, gene in enumerate(GENE):
			if gene in pathway:
				phit+=abs(r[i])
			else:
				pmiss+=1.0/(len(GENE)-len(pathway))
			Phit.append(phit)
			Pmiss.append(pmiss)		
		if phit>0:
			Phit=np.array(Phit)/phit
		P=Phit-np.array(Pmiss)	
		if abs(max(P))>abs(min(P)):
			ES.append(max(P))
		else:
			ES.append(min(P))	
	if n==0:
		ESS=ES		
	ESNull.append(ES)

ESNull=np.sort(np.array(ESNull).transpose()) #transposing the ESNull matrix and sorting values in each row

Pvalues=[] #reserved for nominal P-values, each row corresponds to a different pathway

for i, es in enumerate(ESS):
	Pvalues.append(float(list(ESNull[i]).index(es))/N) #Nominal P-values computed from the left portion of the histogram

# proposed OUTPUT: a table (CSV), first column:PathwayNames, second column:ESS, third column:Pvalues

