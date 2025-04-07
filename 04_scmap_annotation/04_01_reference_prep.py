'''
author: Amparo Roig Adam
created: 2022/02/03
description: prepare files for cell type annotation based on filtered Allen Mouse Brain dataset (filter Ssp)
inputs: metadata and expression matrix from Allen Mouse Brain 
returns: loom file with filtered dataset to read in R
'''

import loompy
import pandas as pd
import numpy as np
import h5py

####Dataset
#path to reference
PATH_REF = 'STXBP1/data/ALLEN_RefData/'

#ref data matrix
filename = PATH_REF+'expression_matrix.hdf5'

#collect metadata
metadata = pd.read_csv(PATH_REF+'metadata.csv')
#print #cells per region
metadata.region_label.value_counts()

#filter only SSp
metadata = metadata.loc[metadata.region_label == 'SSp',:]

metadata.loc[metadata.index.tolist()].region_label.value_counts()

with h5py.File(filename, "r") as f:
    # List all groups
    print("Keys: %s" % f['data'].keys())
    a_group_key = list(f.keys())[0]
    print(type(f['data']['samples']))

    # gene names are collected as numpy.bytes type, need to decode it to get the proper str
        #genes = [g.decode('UTF-8') for g in f['data']['gene']]
    genes = f['data']['gene'].asstr()[:]    #asstr makes the process faster and more reliable (not dependent on knowing the type of encoding - UTF-8)

    #find cells filtered in metadata (S1 cells) - need to be decoded as the genes
    samples = f['data']['samples'].asstr()[:]

    orig_indices = samples.argsort()
    ind = orig_indices[np.searchsorted(samples[orig_indices], metadata.sample_name)]

    #print(ind)
    samples = samples[ind[ind.argsort()]]   ##added argsort here to have same order as in count matrix
    samples[0]

    #print(samples[0:5])
    # Get the count data with same index as samples
    counts = f['data']['counts'][:,ind[ind.argsort()]]
    
    #add metadata to colattrs dictionary
    colattrs = {'CellID':samples}

    #order metadata according to current sample name order
    metadata_sort_index = np.array(metadata.sample_name).argsort()
    metadata_samples_order = metadata_sort_index[np.searchsorted(np.array(metadata.sample_name)[metadata_sort_index], samples)]

    #metadata needs to be added as dictionary to loom file
        #dictionary needs to be ordered by cell name with the same order as 'samples' and the count matrix
    metadict = dict(metadata.iloc[metadata_samples_order])
    [metadict.update({k:list(metadict[k])}) for k in metadict.keys()]

    #order is preserved - checkpoint
    (metadict['sample_name'] == metadata.loc[metadata.index.tolist()].sample_name).value_counts()

    colattrs.update(metadict)
    colattrs.keys()

    #create loom with
    loompy.create("STXBP1/data/ALLEN_RefData/20220217_ref_ABM_S1_10x.loom", counts, row_attrs = {'Gene': genes}, col_attrs = colattrs)

print(samples)

with loompy.connect("data/ALLEN_RefData/20220217_ref_ABM_S1_10x.loom") as ds:
    print(ds.ca.keys())

    print(ds.ca.CellID)
    print(ds.ca.sample_name)

    print(pd.DataFrame(ds.ca.CellID == ds.ca.sample_name).value_counts())   ###this checks that cell names and metadata matches
