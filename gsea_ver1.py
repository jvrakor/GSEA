# this is a Python GSEA implementation. Input files are leukemia.txt and pathways.txt. The implementation follows the original article
# available at http://www.pnas.org/content/102/43/15545.abstract. p is set to p=1 and thus omitted. Output is a CSV file gsea.csv. 
# It is a table containing pathway names, ES(S) values, NES(S) values and nominal P values.

import numpy as np
import time

N=1000 #number of permutations

PH1=[] #reserved for ALL phemotyes data
PH2=[] #reserved for AML phenotypes data
GENE=[] #reserved for gene names from leukemia.txt

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

PathwayGenes=[] #reserved for gene names 
PathwayName=[] #reserved for pathway names

with open('pathways.txt', 'r') as f:
	lines = f.readlines()

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
ESS=[] #corresponds to original data ES(S) values

for n in range(N+1):
	ES=[] #reserved for ES(S) scores
	
	#assigning random original phenotype names to samples, but permuting values instead of names
	if n!=0:
		PH=PH[:, np.random.permutation(PH.shape[1])]
		PH1=PH[:,:Lph1] #the first portion of a combined array, up to the number of ALL phenotype samples
		PH2=PH[:,Lph1:] #the second portion, corresponding to AML phenotype
	
	#construction of correlation values array r, Signal2noise method is used
	r=[(np.mean(val1)-np.mean(val2))/(np.std(val1)+np.std(val2)) for (val1,val2) in zip(PH1,PH2)]
	
	#sorting of gene names and r according to r value and resorting them so we start with highest correlation first
	combined=np.array(sorted(zip(r,GENE), reverse=True))

	#COMPUTING ES(S) for each pathway S
	for pathway in PathwayGenes:
		phit=0.	
		pmiss=0.
		Phit=[]
		Pmiss=[]
		
		pathway=np.array(pathway)
		
		#find the position of genes in GENE that are in pathway
		matches=np.where(np.in1d(combined[:,1],pathway))[0]
		
		if len(matches)>0 and matches[0]>0:
			pmiss=float(matches[0])
			Pmiss.append(pmiss)
			Phit.append(phit)		
		for j,i in enumerate(matches):
			Pmiss.append(pmiss)
			phit+=abs(float(combined[i,0]))
			Phit.append(phit)
			try:
				pmiss+=float(matches[j+1]-i-1)
				Pmiss.append(pmiss)
				Phit.append(phit)
			except IndexError:
				if i==len(GENE)-1:
					break
				pmiss+=float(len(GENE)-i-1)	
				Pmiss.append(pmiss)
				Phit.append(phit)
			
		if phit>0:
			Phit=np.array(Phit)/phit
		Pmiss=np.array(Pmiss)/(len(GENE)-len(matches)) #corrected method that adjusts ES(S) to be in [-1,1]
		P=Phit-Pmiss
		try:
			if abs(max(P))>abs(min(P)):
				ES.append(max(P))
			else:
				ES.append(min(P))	
		except ValueError: #indicates no matches
			ES.append(-1.0)								
	if n==0:
		ESS=ES		
	ESNull.append(ES)

ESNull=np.sort(np.array(ESNull).transpose()) #transposing the ESNull matrix and sorting values in each row

Pvalues=[] #reserved for nominal P-values, each value corresponds to a different pathway
NES=[] #reserved for normalized ES(S) scores
for i, es in enumerate(ESS):
	mark=np.searchsorted(ESNull[i],0)
	if es>0:
		Pvalues.append(1.0-float(np.where(ESNull[i]==es)[0][0])/N)
		NES.append(es/np.mean(ESNull[i][mark:]))
	else:
		Pvalues.append(float(np.where(ESNull[i]==es)[0][0])/N)
		NES.append(es/abs(np.mean(ESNull[i][:mark])))	
		
f=open('gsea.csv', 'w')
f.write('Pathway, ES(S), NES(S), Nominal P value'+'\n')
for i in np.argsort(NES)[::-1]:
	f.write(str(PathwayName[i])+','+str(ESS[i])+','+str(NES[i])+','+str(Pvalues[i])+'\n')
f.close()
