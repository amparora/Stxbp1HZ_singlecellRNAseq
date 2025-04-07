'''
author:  Amparo Roig Adam with support from Lisa Bast
created: 2022/02/20
description: Select high variable genes from ref_ABM_S1_10x reference dataset for automated cell type annotation
inputs: Allen Brain reference dataset in loom format and train set cell names used in HDO genes calculation (04_03_scmap_HDOgenes.R)
returns: text files for HVGs calculated for scmap parameter screen
'''

import scanpy as sc
import os

os.getcwd()

#read reference Dataset
ref_ABM_S1_10x = sc.read_loom("data/ALLEN_RefData/20220217_ref_ABM_S1_10x.loom")
ref_ABM_S1_10x

    #use the following block for GABA and Glutamatergic specific cell annotation runs
train_ref_ABM_S1_10x = ref_ABM_S1_10x[ref_ABM_S1_10x.obs.class_label == 'Glutamatergic',]

#read training set cell names - training set is randomly selected in a previous step (04_03_scmap_HDOgenes.R). Same training set is taken for consistency
train_cells =  open('results/20220217_ALL_ref_scmap_lognorm/cellID_train_ref_ABM_S1_10x.txt','r').read().splitlines() # NOTE - these should be only GABA or Glutamatergic in final annotation run
train_cells

#filter refSet to trainig set - NOTE - do not run this line to get genes for final GABA/Glutamatergic annotation
train_ref_ABM_S1_10x = ref_ABM_S1_10x[train_cells,]
train_ref_ABM_S1_10x

#NORMALISATION
    #remove 0 expressed genes
sc.pp.filter_genes(train_ref_ABM_S1_10x, min_cells=1, inplace=True)
train_ref_ABM_S1_10x0 = train_ref_ABM_S1_10x

for nHVGs in [500,1000,2000, 3000, 5000]: 
    train_ref_ABM_S1_10x = train_ref_ABM_S1_10x0        #reset training dataset
    for meth in ['seurat_v3', 'cell_ranger']: 
        print(nHVGs, meth)

        #only cell_ranger method expects log normalised counts. seurat method expects raw counts.
        if meth == 'cell_ranger':
            sc.pp.normalize_total(train_ref_ABM_S1_10x)
            sc.pp.log1p(train_ref_ABM_S1_10x)

        sc.pp.highly_variable_genes(train_ref_ABM_S1_10x, n_top_genes = nHVGs, inplace = True, flavor = meth )

        HVGs = list(train_ref_ABM_S1_10x.var.index[train_ref_ABM_S1_10x.var.highly_variable == True])

        #save HVGs
        hvgs_file = 'results/20220217_ALL_ref_scmap_lognorm/'+str(nHVGs)+'_'+meth+'_HVGs_ref_ABM_S1_10x.txt'
        print(hvgs_file)

        with open(hvgs_file, 'w')as f:
            for HVG in HVGs:
                f.write(HVG + '\n')

