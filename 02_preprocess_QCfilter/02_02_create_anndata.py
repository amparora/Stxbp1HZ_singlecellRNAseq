'''
author: Amparo Roig Adam
created: 2021/02/20
description: Add all concerning metadata to loom file containing all scRNA-seq data aggregated

input: all samples aggregated scRNA-seq counts in loom file
output: anndata file with counts and metadata
 '''
 
#Packages needed
import loompy
import anndata
import pandas as pd
import numpy as np
from datetime import date

ds = loompy.connect('data/seq_data/aggregated_stxbp1.loom')
ds.shape
ds.attrs.title = 'STXBP1 scRNA-seq'

ds.ca.keys()

#add sample ID info to ds.ca

    #count number of times each suffix appears and create list with that many repetitions of sample ID
suffix_l = [cellID.split('-')[1] for cellID in ds.ca.CellID]

aggregation = pd.read_csv("data/seq_data/agg_stxbp1_libraries_c5.csv")

SampleID = []
for suffix in range(1,15):
    ncells = suffix_l.count(str(suffix))
    SampleID.append([aggregation['library_id'][suffix - 1]]*ncells)
    print(f"nCells sample {aggregation['library_id'][suffix - 1]} :   ",ncells)

#flatten list
SampleID = [item for sublist in SampleID for item in sublist]

ds.ca['SampleID'] = SampleID

ds.close()

#create anndata file
adataSTXBP1 = anndata.read_loom('data/seq_data/aggregated_stxbp1.loom')
adataSTXBP1.var_names_make_unique()
adataSTXBP1.var

#Add metadata - sex, ageweek, genotype, tissue extraction, live cells percent

adataSTXBP1.obs[['Sex', 'AgeWeek', 'Genotype', 'TissueExtraction', 'LiveCellsPercent']] = 'NA'
adataSTXBP1.obs

metadata = pd.read_csv("data/01_mouse_info.csv", sep = ';')

for sample in adataSTXBP1.obs['SampleID'].unique():
    print('Adding metadata to sample:   ',sample)

        #for each sample, the obs slots are modified with the corresponding information in the metadata dataframe
    adataSTXBP1.obs.loc[adataSTXBP1.obs['SampleID'] == sample,['Sex', 'AgeWeek', 'Genotype', 'TissueExtraction', 'LiveCellsPercent']] = metadata[metadata.Linear == float('1711'+sample[-2:])][['Sex', 'AgeWeek', 'Genotype', 'TissueExtraction', 'LiveCellsPercent']].values

#Save full dataset

#All metadata is stored as categorical, for proper data handling downstream, we need to reconvert the pertinent columns into float and re-save
adataSTXBP1.obs['LiveCellsPercent'] = np.float32(adataSTXBP1.obs['LiveCellsPercent'])
adataSTXBP1.obs['AgeWeek'] = np.float32(adataSTXBP1.obs['AgeWeek'])

#type(adataSTXBP1.obs['AgeWeek'][0])
adataSTXBP1.write('data/seq_data/'+date.today().strftime("%Y%m%d")+'_stxbp1_original'+'.h5ad')

print('Saved file:\n', 'data/seq_data/'+date.today().strftime("%Y%m%d")+'_stxbp1_original'+'.h5ad')
