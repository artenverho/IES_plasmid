# IES_plasmid

Goal: Identify features in the meta-genome and 16s data that can aid the design of the interspecific expression system.

    Investigate the presence of antibiotic resistance genes
    Assess the microbial species abundance
    Search for estrogen degradation genes in the MAGs

 

Useful information:

Manual for using the BioServers is attached to this experiment

Antibiotic resistance genes (ARGs) identification in WWTP samples:

The Comprehensive Antibiotic Resistance Database (CARD; https://card.mcmaster.ca)

ResFinder.

 

Estrogen degradation pathway:

    oecA: https://www.genome.jp/dbget-bin/www_bget?spkc:KC8_09390
    oecB: https://www.genome.jp/dbget-bin/www_bget?spkc:KC8_16650
    oecC: https://www.genome.jp/dbget-bin/www_bget?spkc:KC8_05325

Proteins of 24000 bacteria genomes in GTDB:

https://data.ace.uq.edu.au/public/misc_downloads/refinem/

 

Tools:

Diamond: http://www.diamondsearch.org/index.php

MAFFT: https://mafft.cbrc.jp/alignment/software/

trimAI: http://trimal.cgenomics.org/

fasttree: http://www.microbesonline.org/fasttree/

iqtree: http://www.iqtree.org/

GeneTreeTk: https://github.com/dparks1134/GeneTreeTk

Install Conda and change package directory:

`module load Miniconda3`
`mkdir software`
`mkdir software/condaenvs`
`conda config --append envs_dirs /srv/PHN/users/lcv/software/condaenvs`

Setup abricate:

`conda init bash
conda activate abricate
conda install -c conda-forge -c bioconda -c defaults abricate`

Use abricate in a screen environment

`screen -S abricate #start a new screen session called abricate
screen -r #reattach screen instance`

Start abricate:

`abricate --db card --threads 12 /*.fna  > output.tab #keep in mind to first create output.txt`
	  	  	  	  	  	  	  	 
Using PANDAS to combine the MAG data from Caitlin with the results from Abricate:

Import of both the MAG data and the output.tab data from abricate:

`import pandas as pd`
`df = pd.read_csv('MAG.statistics.STABLEX.20200213.tsv', delimiter='\t')`
`abr = pd.read_csv('output.tab', delimiter='\t')`
`print(abr)`

Trim the MAG filename in abr and df datasets:

`abr['#FILE'] = abr['#FILE'].map(lambda x: x.lstrip('/srv/PHN/users/cms/for_colleagues/Laura/bins/').rstrip('.fna'))`

`df['MAG'] = df['MAG'].map(lambda x: x.rstrip('.fa'))`

Rename first column of abr to match df:

`abr.rename(columns={'#FILE':'MAG'}, inplace=True)`

Set index on first column (MAG). This allows to use the function join (for information check: https://pandas.pydata.org/pandas-docs/stable/user_guide/merging.html):

`df = df.set_index('MAG')`
`abr = abr.set_index('MAG')`

`dfabr2 = df.join(abr, how='outer')`

Make a summary of the results and filter only Aalborg West WWTP:

`sumar2 = dfabr2[['max_abund','GTDBTax','MiDAS3_6Tax','RESISTANCE']]`

`AAV = sumar2.filter(like='AalW', axis=0)`

Export data:

`sumar2 = dfabr2[['max_abund','GTDBTax','MiDAS3_6Tax','RESISTANCE']]`

`dfabr2.to_excel("outputdfabr2.xlsx") sumar2.to_excel("outputsumar2.xlsx") AAV.to_excel("AAV.xlsx")`

Use the R-script for AmpVis2 to obtain abundance data from WWTPs:

Adjust the folder descriptions: 

`setwd("/Users/Maarten/R")`
usearch_otutable <- amp_import_usearch(otutab = "/Users/Maarten/R/ASVtable.tsv",
sintax = "/Users/Maarten/R/ASVs.R1.midas35.sintax")
#from Kasper: use this way of loading metadata for PCA to work

library(openxlsx) metadata <- openxlsx::read.xlsx("/Users/Maarten/R/metadata_BioBank_2019-11-27.xlsx", detectDates = TRUE)`

 Script to match the species names from the replication origin table with the list the tree of ARB:

`import pandas as pd
import re
from Bio import AlignIO
import numpy as np
alignment = AlignIO.read("gtdbtk.bac120.msa.fasta", "fasta")
edited = alignment[1, :] 
repori = pd.read_csv('rep_origins.csv')
repori['wordcount'] = repori['Species'].str.count(' ') + 1
length = len(alignment)
results = []
for i in range(1 , length):
table = alignment[i, :] 
desc = table.description
results.append(desc)`

`#print(edited)
#dir(edited)
#edited.description
df = pd.DataFrame(results, columns=['description'])
df[['ID','description']] = df["description"].str.split(" ", 1, expand=True)
#split of the species name:
df['species']= df['description'].str.split('s__').str[1]
#make sure there are not NaN values (appeariantly needed to let the script run):
df = df.fillna('maarten')
#df[df['species'].str.contains('Pseudomonas')]
#script to match the species names in repori to the list from the Msa:
import numpy as np
i = 0
df['MaartenID'] = '0'
for i in range(0, len(repori)):
repspec = repori.iloc[i,1] 
rep_ID = repori.iloc[i,0] 
rep_ID
repspec
if repori.iloc[i,2] == 1:
    df['MaartenID'] = np.where(df['species'].str.contains(repspec), rep_ID, df['MaartenID'])
else:
    df['MaartenID'] = np.where(df['species'].str.match(repspec), rep_ID, df['MaartenID'])

#df[df['species'].str.contains('Xanthomonas')]

df = df.replace({'0':None})`

 

 
